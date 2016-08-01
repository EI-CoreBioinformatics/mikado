# coding: utf-8

"""
This module defines the last object to be created during the picking,
i.e. the locus.
"""

import itertools
import operator
from .transcript import Transcript
from ..scales.assigner import Assigner
from .sublocus import Sublocus
from .abstractlocus import Abstractlocus
from ..parsers.GFF import GffLine
import collections
from sys import version_info
if version_info.minor < 5:
    from sortedcontainers import SortedDict
else:
    from collections import OrderedDict as SortedDict


class Locus(Sublocus, Abstractlocus):
    """Class that defines the final loci.
    It is a child of monosublocus, but it also has the possibility of adding
    additional transcripts if they are valid splicing isoforms.
    """

    def __init__(self, transcript: Transcript, logger=None):
        """
        Constructor class. Like all loci, also Locus is defined starting from a transcript.

        :param transcript: the transcript which is used to initialize the Locus
        :type transcript: Transcript

        :param logger: the logger instance.
        :type logger: None | logging.Logger
        """

        self.counter = 0
        transcript.attributes["primary"] = True
        self.counter = 0  # simple tag to avoid collisions
        Abstractlocus.__init__(self)
        # this must be defined straight away
        self.monoexonic = transcript.monoexonic
        Abstractlocus.add_transcript_to_locus(self, transcript)
        self.locus_verified_introns = transcript.verified_introns
        self.metrics_calculated = False
        self.scores_calculated = False
        self.score = transcript.score
        # A set of the transcript we will ignore during printing
        # because they are duplications of the original instance. Done solely to
        # get the metrics right.
        self.__orf_doubles = collections.defaultdict(set)
        self.excluded = None
        self.parent = None
        self.tid = transcript.id
        self.logger = logger
        self.attributes = dict()
        self.logger.debug("Created Locus object with {0}".format(transcript.id))
        self.primary_transcript_id = transcript.id
        self.attributes["is_fragment"] = False
        self.metric_lines_store = []
        self.__id = None

    def __str__(self, print_cds=True) -> str:

        self.feature = self.__name__
        assert self.feature != "Monosublocus"
        # Hacky fix to make sure that the primary transcript has the attribute
        # Set to True in any case.
        self.primary_transcript.attributes["primary"] = True
        # BF, really just a hack.
        for transcript in self.transcripts:
            if transcript == self.primary_transcript_id:
                continue
            self.transcripts[transcript].attributes["primary"] = False

        lines = []

        self_line = GffLine('')
        for attr in ["chrom", 'feature', 'source', 'start', 'end', 'strand']:
            setattr(self_line, attr, getattr(self, attr))
        self_line.phase, self_line.score = None, self.score
        self_line.id = self.id
        self_line.name = self.name
        self_line.attributes["superlocus"] = self.parent
        self_line.attributes.update(self.attributes)
        if "is_fragment" in self.attributes and self.attributes["is_fragment"] is False:
            del self_line.attributes["is_fragment"]
        self_line.attributes["multiexonic"] = (not self.monoexonic)
        lines.append(str(self_line))

        for tid in self.transcripts:
            transcript_instance = self.transcripts[tid]
            transcript_instance.source = self.source
            transcript_instance.parent = self_line.id
            self.logger.debug(self.attributes)
            for attribute in self.attributes:
                if attribute not in transcript_instance.attributes:
                    if attribute == "is_fragment" and self.attributes[attribute] is False:
                        continue
                    transcript_instance.attributes[attribute] = self.attributes[attribute]

            lines.append(transcript_instance.format(
                "gff", with_cds=print_cds,
                all_orfs=self.json_conf["pick"]["output_format"]["report_all_orfs"]
            ).rstrip())

        return "\n".join(lines)

    def add_transcript_to_locus(self, transcript: Transcript, **kwargs):
        """Implementation of the add_transcript_to_locus method.
        Before a transcript is added, the class checks that it is a valid splicing isoform
        and that we have not exceeded already the maximum number of isoforms for the Locus.

        :param transcript: the candidate transcript
        :type transcript: Transcript

        :param kwargs: optional keyword arguments are ignored.
        """

        _ = kwargs
        to_be_added = True
        # Total, 5', 3'
        max_utr_lenghts = {
            "total": self.json_conf["pick"]["alternative_splicing"]["max_utr_length"],
            "five": self.json_conf["pick"]["alternative_splicing"]["max_fiveutr_length"],
            "three": self.json_conf["pick"]["alternative_splicing"]["max_threeutr_length"]}
        max_isoforms = self.json_conf["pick"]["alternative_splicing"]["max_isoforms"]

        if len(self.transcripts) >= max_isoforms:
            self.logger.debug("%s not added because the Locus has already too many transcripts.",
                              transcript.id)
            to_be_added = False

        if self.json_conf["pick"]["alternative_splicing"]["only_confirmed_introns"] is True:
            to_check = transcript.introns - self.primary_transcript.introns
            to_check -= transcript.verified_introns
            if len(to_check) > 0:
                self.logger.debug(
                    "%s not added because it has %d non-confirmed intron%s",
                    transcript.id,
                    len(to_check),
                    "s" * min(1, len(to_check) - 1))
                to_be_added = False

        if (to_be_added and transcript.score <
                self.primary_transcript.score *
                self.json_conf["pick"]["alternative_splicing"]["min_score_perc"]):
            self.logger.debug(
                "%s not added because its score (%.2f) is less \
                than %.2f%% of the primary score (%.2f)",
                transcript.id,
                transcript.score,
                self.json_conf["pick"]["alternative_splicing"]["min_score_perc"] * 100,
                self.primary_transcript.score)
            to_be_added = False

        if to_be_added and transcript.strand != self.strand:
            self.logger.debug("%s not added because it has a different strand from %s (%s vs. %s)",
                              transcript.id, self.id, transcript.strand, self.strand)
            to_be_added = False

        if to_be_added:
            is_alternative, ccode, comparison = self.is_alternative_splicing(transcript)
            if is_alternative is False:
                self.logger.debug("%s not added because it is not a \
                valid splicing isoform. Ccode: %s",
                                  transcript.id, ccode)
                to_be_added = False
            else:
                transcript.attributes["ccode"] = ccode
            if self.json_conf["pick"]["alternative_splicing"]["min_cdna_overlap"] > 0:
                overlap = comparison.n_recall[0]
                if overlap < self.json_conf["pick"]["alternative_splicing"]["min_cdna_overlap"]:
                    self.logger.debug(
                        "%s not added because its CDNA overlap is too low (%f%%).",
                        transcript.id,
                        round(overlap * 100, 2))
                    to_be_added = False

        if to_be_added and transcript.combined_utr_length > max_utr_lenghts["total"]:
            self.logger.debug("%s not added because it has too much UTR (%d).",
                              transcript.id,
                              transcript.combined_utr_length)
            to_be_added = False
        if to_be_added and transcript.five_utr_length > max_utr_lenghts["five"]:
            self.logger.debug("%s not added because it has too much 5'UTR (%d).",
                              transcript.id,
                              transcript.five_utr_length)
            to_be_added = False
        if to_be_added and transcript.three_utr_length > max_utr_lenghts["three"]:
            self.logger.debug("%s not added because it has too much 3'UTR (%d).",
                              transcript.id,
                              transcript.three_utr_length)
            to_be_added = False

        self.find_retained_introns(transcript)
        if (to_be_added and
                self.json_conf["pick"]["alternative_splicing"]["keep_retained_introns"] is False):
            if transcript.retained_intron_num > 0:
                self.logger.debug("%s not added because it has %d retained introns.",
                                  transcript.id,
                                  transcript.retained_intron_num)
                to_be_added = False

        if to_be_added and self.json_conf["pick"]["alternative_splicing"]["min_cds_overlap"] > 0:
            if self.primary_transcript.combined_cds_length > 0:
                tr_nucls = set(itertools.chain(
                    *[range(x[0], x[1] + 1) for x in transcript.combined_cds]))
                primary_nucls = set(
                    itertools.chain(
                        *[range(x[0], x[1] + 1) for x in self.primary_transcript.combined_cds]))
                nucl_overlap = len(set.intersection(primary_nucls, tr_nucls))
                overlap = nucl_overlap / self.primary_transcript.combined_cds_length
                if overlap < self.json_conf["pick"]["alternative_splicing"]["min_cds_overlap"]:
                    self.logger.debug(
                        "%s not added because its CDS overlap is too low (%f%%).",
                        transcript.id,
                        round(overlap * 100, 2))
                    to_be_added = False

        if to_be_added is False:
            return

        self.logger.debug("Keeping %s as a valid alternative isoform for %s",
                          transcript.id, self.id)
        transcript.attributes["primary"] = False

        Abstractlocus.add_transcript_to_locus(self, transcript)
        self.locus_verified_introns.update(transcript.verified_introns)

    def is_intersecting(self, *args):
        """Not implemented: this function makes no sense for a single-transcript container.
        :param args: any argument to this nethod will be ignored.
        """
        raise NotImplementedError("""Loci do not use this method, but rather
        assess whether a transcript is a splicing isoform or not.""")

    def other_is_fragment(self,
                          other,
                          minimal_cds_length=0,
                          minimal_exons=2):
        """
        :param other: another Locus to compare against
        :type other: Locus

        :param minimal_cds_length: Minimal CDS length to consider
        a Locus as non-fragment, no matter the ccode.
        :type minimal_cds_length: int

        :param minimal_exons: Maximum number of exons to consider a
        a Locus as non-fragment, no matter the ccode.
        :type minimal_cds_length: int

        This function checks whether another *monoexonic* Locus
        *on the opposite strand* is a fragment,by checking its classification
        according to Assigner.compare.
        Briefly, a transcript is classified as fragment
        if it follows the following criteria:

            - it is monoexonic
            - it has a combined_cds_length inferior to maximal_cds
            - it is classified as x,i,P
        """

        if not isinstance(self, type(other)):
            raise TypeError("I can compare only loci.")

        self.logger.debug("Comparing %s with %s",
                          self.primary_transcript_id,
                          other.primary_transcript_id)

        if other.primary_transcript.combined_cds_length > minimal_cds_length:
            self.logger.debug("%s has a CDS of %d, not a fragment by definition",
                              other.primary_transcript_id,
                              other.primary_transcript.combined_cds_length)
            return False
        elif other.primary_transcript.exon_num > minimal_exons:
            self.logger.debug("%s has %d exons, not a fragment by definition",
                              other.primary_transcript_id,
                              other.primary_transcript.combined_cds_length)
            return False

        result, _ = Assigner.compare(other.primary_transcript, self.primary_transcript)
        # Exclude anything which is completely contained within an intron,
        # or is a monoexonic fragment overlapping/in the neighborhood
        self.logger.debug("Comparison between {0} (strand {3}) and {1}: class code \"{2}\"".format(
            self.primary_transcript.id,
            other.primary_transcript.id,
            result.ccode[0],
            other.strand))
        if result.ccode[0] in ("i", "P", "p", "x"):
            self.logger.debug("{0} is a fragment (ccode {1})".format(
                other.primary_transcript.id, result.ccode[0]))
            return True
        # Adding c's because fragments might very well be contained!
        elif other.strand is None and (result.n_f1[0] > 0 or result.ccode in ("rI", "ri")):
            self.logger.debug("Unstranded {0} is a fragment (ccode {1})".format(
                other.primary_transcript.id, result.ccode[0]))
            return True

        return False

    def set_json_conf(self, jconf: dict):
        """
        Setter for the configuration dictionary.
        :param jconf:
        :type jconf: dict
        """
        if not isinstance(jconf, dict):
            raise TypeError("Invalid configuration of type {0}".format(type(jconf)))
        self.json_conf = jconf

    def get_metrics(self):

        """Quick wrapper to calculate the metrics for all the transcripts."""

        # TODO: Find an intelligent way ot restoring this check

        if self.metrics_calculated is True:
            return

        # self.logger.info("Calculating the intron tree for %s", self.id)
        assert len(self._cds_introntree) == len(self.combined_cds_introns)

        for tid in sorted(self.transcripts):
            self.calculate_metrics(tid)

        self.logger.debug("Finished to calculate the metrics for %s", self.id)

        self.metrics_calculated = True
        return

    def calculate_metrics(self, tid: str):
        """
        :param tid: the name of the transcript to be analysed
        :type tid: str

        This function will calculate the metrics for a transcript which are relative in nature
        i.e. that depend on the other transcripts in the sublocus. Examples include the fraction
        of introns or exons in the sublocus, or the number/fraction of retained introns.
        """

        self.logger.debug("Calculating metrics for %s", tid)
        self.transcripts[tid].finalize()
        if (self.transcripts[tid].number_internal_orfs <= 1 or
                    self.json_conf["pick"]["output_format"]["report_all_orfs"] is False):
            super().calculate_metrics(tid)
        else:
            transcript = self.transcripts[tid]
            selected = transcript.selected_internal_orf
            new_transcript = transcript.copy()
            new_transcript.id = "{0}.orf1".format(new_transcript.id)
            self.transcripts[new_transcript.id] = new_transcript
            super().calculate_metrics(new_transcript.id)
            self.__orf_doubles[tid].add(new_transcript.id)

            for num, orf in enumerate([_ for _ in transcript.internal_orfs if
                                       _ != selected]):
                new_transcript = transcript.copy()
                assert isinstance(new_transcript, Transcript)
                new_transcript.internal_orfs = [orf]
                new_transcript.internal_orfs.extend([_ for _ in transcript.internal_orfs if
                                                     orf != _])
                new_transcript.id = "{0}.orf{1}".format(new_transcript.id, num + 2)
                self.transcripts[new_transcript.id] = new_transcript
                super().calculate_metrics(new_transcript.id)
                self.__orf_doubles[tid].add(new_transcript.id)

        self.logger.debug("Calculated metrics for {0}".format(tid))

    def calculate_scores(self):
        """
        Function to calculate a score for each transcript, given the metrics derived
        with the calculate_metrics method and the scoring scheme provided in the JSON configuration.
        If any requirements have been specified, all transcripts which do not pass them
        will be assigned a score of 0 and subsequently ignored.
        Scores are rounded to the nearest integer.
        """

        if self.scores_calculated is True:
            return

        self.get_metrics()
        if not hasattr(self, "logger"):
            self.logger = None
            self.logger.setLevel("DEBUG")
        self.logger.debug("Calculating scores for {0}".format(self.id))

        self.scores = dict()
        for tid in self.transcripts:
            self.scores[tid] = dict()
            # Add the score for the transcript source
            self.scores[tid]["source_score"] = self.transcripts[tid].source_score

        if self.regressor is None:
            for param in self.json_conf["scoring"]:
                self._calculate_score(param)

            for tid in self.scores:
                self.transcripts[tid].scores = self.scores[tid].copy()

            for tid in self.transcripts:

                if tid in self.__orf_doubles:
                    del self.scores[tid]
                    continue
                self.transcripts[tid].score = sum(self.scores[tid].values())
                self.scores[tid]["score"] = self.transcripts[tid].score

        else:
            valid_metrics = self.regressor.metrics
            metric_rows = SortedDict()
            for tid, transcript in sorted(self.transcripts.items(), key=operator.itemgetter(0)):
                for param in valid_metrics:
                    self.scores[tid][param] = "NA"
                row = []
                for attr in valid_metrics:
                    val = getattr(transcript, attr)
                    if isinstance(val, bool):
                        if val:
                            val = 1
                        else:
                            val = 0
                    row.append(val)
                metric_rows[tid] = row
            # scores = SortedDict.fromkeys(metric_rows.keys())
            for pos, score in enumerate(self.regressor.predict(list(metric_rows.values()))):
                tid = list(metric_rows.keys())[pos]
                if tid in self.__orf_doubles:
                    del self.scores[tid]
                    continue
                self.scores[tid]["score"] = score
                self.transcripts[tid].score = score

        self.metric_lines_store = []
        for row in self.prepare_metrics():
            if row["tid"] in self.__orf_doubles:
                continue
            else:
                self.metric_lines_store.append(row)

        for doubled in self.__orf_doubles:
            for partial in self.__orf_doubles[doubled]:
                if partial in self.transcripts:
                    del self.transcripts[partial]

        self.scores_calculated = True

    def print_scores(self):
        """This method yields dictionary rows that are given to a csv.DictWriter class."""
        self.calculate_scores()
        if self.regressor is None:
            score_keys = sorted(list(self.json_conf["scoring"].keys()) + ["source_score"])
        else:
            score_keys = sorted(self.regressor.metrics + ["source_score"])
        keys = ["tid", "parent", "score"] + score_keys

        for tid in self.scores:
            row = dict().fromkeys(keys)
            row["tid"] = tid
            row["parent"] = self.id
            row["score"] = round(self.scores[tid]["score"], 2)
            calculate_total = (self.regressor is None)
            for key in score_keys:
                if calculate_total:
                    assert self.scores[tid][key] != "NA" and self.scores[tid][key] is not None
                    row[key] = round(self.scores[tid][key], 2)
            if calculate_total is True:
                score_sum = sum(row[key] for key in score_keys)
                #
                assert round(score_sum, 2) == round(self.scores[tid]["score"], 2), (
                    score_sum,
                    self.transcripts[tid].score,
                    tid)
            yield row

    def is_alternative_splicing(self, other):

        """This function defines whether another transcript could be a
        putative alternative splice variant of the primary Locus
        transcript.
        To do so, it compares the candidate against all transcripts in the Locus, and calculates
        the class code using scales.Assigner.compare.
        If all the matches are "n" or "j", the transcript is considered as an AS event.

        :param other: another transcript to compare against
        :type other: Transcript

        """

        is_valid = True
        # main_ccode = None

        valid_ccodes = self.json_conf["pick"]["alternative_splicing"]["valid_ccodes"]
        redundant_ccodes = self.json_conf["pick"]["alternative_splicing"]["redundant_ccodes"]

        main_result, _ = Assigner.compare(other, self.primary_transcript)
        main_ccode = main_result.ccode[0]

        if main_ccode not in valid_ccodes:
            self.logger.debug("%s is not a valid splicing isoform. Ccode: %s",
                              other.id,
                              main_result.ccode[0])
            is_valid = False
        if is_valid:
            for tid in iter(tid for tid in self.transcripts if
                            tid not in (self.primary_transcript_id, other.id)):
                candidate = self.transcripts[tid]
                result, _ = Assigner.compare(other, candidate)
                if result.ccode[0] in redundant_ccodes:
                    self.logger.debug("%s is a redundant isoform of %s (ccode %s)",
                                      other.id, candidate.id, result.ccode[0])
                    is_valid = False
                    break

        return is_valid, main_ccode, main_result

    @property
    def __name__(self):
        if len(self.transcripts) == 0:
            return "locus"
        elif any(transcript.selected_cds_length > 0 for
                 transcript in self.transcripts.values()):
            return "gene"
        else:
            return "ncRNA_gene"

    # pylint: disable=invalid-name
    @property
    def id(self):
        """
        Override of the abstractlocus method.
        :rtype str
        """
        if self.__id is not None:
            return self.__id
        else:
            myid = Abstractlocus.id.fget(self)  # @UndefinedVariable

            if self.counter > 0:
                myid = "{0}.{1}".format(myid, self.counter)
            return myid

    # pylint: disable=arguments-differ
    @id.setter
    def id(self, string):
        """
        Override of the original method from AbstractLocus. This override allows to
        create proper IDs for the final annotation to be output by Mikado.
        :param string:
        :return:
        """

        self.logger.debug("Setting new ID for %s to %s", self.id, string)
        self.__id = string
        primary_id = "{0}.1".format(string)
        old_primary = self.primary_transcript.id
        self.primary_transcript.attributes["Alias"] = self.primary_transcript.id
        self.primary_transcript.id = primary_id
        self.transcripts[primary_id] = self.primary_transcript
        self.primary_transcript_id = primary_id
        del self.transcripts[old_primary]

        order = sorted([k for k in self.transcripts.keys() if k != primary_id],
                       key=lambda xtid: self.transcripts[xtid])

        for counter, tid in enumerate(order):
            counter += 2
            self.transcripts[tid].attributes["Alias"] = tid
            new_id = "{0}.{1}".format(string, counter)
            self.transcripts[tid].id = new_id
            self.transcripts[new_id] = self.transcripts.pop(tid)

    # pylint: enable=invalid-name,arguments-differ

    @property
    def is_fragment(self):
        """
        :rtype : bool
        Flag. It returns the value of self.attributes["is_fragment"]

        """
        return self.attributes["is_fragment"]

    @is_fragment.setter
    def is_fragment(self, val: bool):
        """
        Setter for is_fragment. Only boolean values are accepted.

        :param val: flag
        :type val: bool
        """
        if not isinstance(val, bool):
            raise ValueError(val)
        self.attributes["is_fragment"] = val

    @property
    def primary_transcript(self):
        """
        This property returns the primary transcript of the Locus
        (i.e. the one which has been used for creation and which has the highest score).
        :rtype : Transcript
        """
        return self.transcripts[self.primary_transcript_id]
