# coding: utf-8

"""
This module defines the last object to be created during the picking,
i.e. the locus.
"""

import collections
import itertools
import operator
from collections import deque
import pyfaidx
from ..transcripts.transcript import Transcript
from ..transcripts.transcriptchecker import TranscriptChecker
from .abstractlocus import Abstractlocus
from .sublocus import Sublocus
from ..parsers.GFF import GffLine
from ..scales.assigner import Assigner
from ..utilities import overlap


class Locus(Abstractlocus):
    """Class that defines the final loci.
    It is a child of monosublocus, but it also has the possibility of adding
    additional transcripts if they are valid splicing isoforms.
    """

    def __init__(self, transcript: Transcript, logger=None, json_conf=None):
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
        self.fai = None
        self.json_conf = json_conf
        self.__finalized = False
        # if verified_introns is not None:
        #     self.locus_verified_introns = verified_introns

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

    def finalize_alternative_splicing(self):

        """"This method ensures that all the transcripts retained in the locus
        are within the score threshold. This is due to the fact that the score
        changes depending on the transcript considered together; so that a transcript
        that might have scored relatively well on its own will score pretty badly when
        brought inside the locus."""

        if self._finalized is True:
            return

        self.metrics_calculated = False
        self.scores_calculated = False
        self.calculate_scores()
        max_isoforms = self.json_conf["pick"]["alternative_splicing"]["max_isoforms"]

        while True:
            to_keep = {self.primary_transcript_id}
            order = sorted([(tid, self.transcripts[tid].score) for tid in self.transcripts
                            if tid != self.primary_transcript_id],
                           key=operator.itemgetter(1), reverse=True)
            threshold = self.json_conf["pick"]["alternative_splicing"]["min_score_perc"] * self.primary_transcript.score

            for tid, score in order:
                if len(to_keep) == max_isoforms:
                    self.logger.debug(
                        "Discarding {} from the locus because we have \
reached the maximum number of isoforms for the locus".format(
                            ", ".join(list(set.difference(set(self.transcripts.keys()),
                                                          to_keep)))
                        ))
                    break
                if score < threshold:
                    self.logger.debug(
                        "Discarding {} from the locus because their scores are below the threshold ({})".format(
                            ", ".join(list(set.difference(set(self.transcripts.keys()),
                                                          to_keep))),
                            round(threshold, 2)))
                    break
                to_keep.add(tid)

            if to_keep == set(self.transcripts.keys()):
                self.logger.debug("Finished to discard superfluous transcripts from {}".format(self.id))
                break
            else:
                for tid in set.difference(set(self.transcripts.keys()), to_keep):
                    self.remove_transcript_from_locus(tid)
                assert len(self.transcripts) > 0, to_keep
                self.metrics_calculated = False
                self.scores_calculated = False
                self.calculate_scores()

        if self.json_conf["pick"]["alternative_splicing"]["pad"] is True:
            self.pad_transcripts()
            self.metrics_calculated = False
            self.scores_calculated = False
            self.calculate_scores()

        self.logger.debug("Now checking the retained introns for %s", self.id)
        while True:
            to_remove = set()
            for tid, transcript in self.transcripts.items():
                if tid == self.primary_transcript_id:
                    continue
                self.find_retained_introns(transcript)
                if transcript.retained_intron_num > 0:
                    to_remove.add(tid)
                else:
                    continue
            if not to_remove:
                break
            elif self.json_conf["pick"]["alternative_splicing"]["keep_retained_introns"] is False:
                self.logger.debug("Removing {} because they contain retained introns".format(
                    ", ".join(list(to_remove))))
                for tid in to_remove:
                    self.remove_transcript_from_locus(tid)
                self.metrics_calculated = False
                self.scores_calculated = False
                self.calculate_scores()
            elif self.json_conf["pick"]["alternative_splicing"]["keep_retained_introns"] is True:
                for tid in to_remove:
                    self.transcripts[tid].attributes["retained_intron"] = True
                break

        self._finalized = True

        return

    def remove_transcript_from_locus(self, tid: str):

        """Overloading of the AbstractLocus class, in order to ensure that the primary transcript will *not*
        be removed."""

        if tid == self.primary_transcript_id:
            raise KeyError("%s is the primary transcript of %s!" % (tid, self.id))

        super().remove_transcript_from_locus(tid)
        self._finalized = False


    def add_transcript_to_locus(self, transcript: Transcript, **kwargs):
        """Implementation of the add_transcript_to_locus method.
        Before a transcript is added, the class checks that it is a valid splicing isoform
        and that we have not exceeded already the maximum number of isoforms for the Locus.

        The checks performed are, in order:

        #. whether the locus already has the maximum number of acceptable
           isoforms ("max_isoforms")
        #. (optional) whether all the introns *specific to the transcript when
           compared with the primary transcript* are confirmed by external
           validation tools (eg Portcullis)
        #. Whether the score of the proposed AS event has a score over the
           minimum percentage of the primary transcript score (eg if the minimum
           percentage is 0.6 and the primary is scored 20, a model with a score
           of 11 would be rejected and one with a score of 12 would be accepted)
        #. Whether the strand of the candidate is the same as the one of the locus
        #. Whether the AS event is classified (ie has a class code) which is acceptable as valid AS
        #. Whether the transcript shares enough cDNA with the primary transcript ("min_cdna_overlap")
        #. Whether the proposed model has too much UTR
        #. (optional) Whether the proposed model has a retained intron compared
           to the primary, ie part of its 3' non-coding regions overlaps
           one intron of the primary model
        #. Whether the proposed model shares enough CDS with the primary model (min_cds_overlap)

        :param transcript: the candidate transcript
        :type transcript: Transcript

        :param kwargs: optional keyword arguments are ignored.
        """

        _ = kwargs
        to_be_added = True

        if to_be_added and transcript.strand != self.strand:
            self.logger.debug("%s not added because it has a different strand from %s (%s vs. %s)",
                              transcript.id, self.id, transcript.strand, self.strand)
            to_be_added = False

        if self.json_conf["pick"]["alternative_splicing"]["only_confirmed_introns"] is True:
            to_check = (transcript.introns - transcript.verified_introns) - self.primary_transcript.introns
            if len(to_check) > 0:
                self.logger.debug(
                    "%s not added because it has %d non-confirmed intron%s",
                    transcript.id,
                    len(to_check),
                    "s" * min(1, len(to_check) - 1))
                to_be_added = False

        # Add a check similar to what we do for the minimum requirements and the fragments
        if to_be_added and "as_requirements" in self.json_conf:
            if ("compiled" not in self.json_conf["as_requirements"] or
                        self.json_conf["as_requirements"]["compiled"] is None):
                self.json_conf["as_requirements"]["compiled"] = compile(
                    self.json_conf["as_requirements"]["expression"], "<json>",
                    "eval")
            evaluated = dict()
            for key in self.json_conf["as_requirements"]["parameters"]:
                value = getattr(transcript,
                                self.json_conf["as_requirements"]["parameters"][key]["name"])
                evaluated[key] = self.evaluate(
                        value,
                        self.json_conf["as_requirements"]["parameters"][key])
                # pylint: disable=eval-used
            if eval(self.json_conf["as_requirements"]["compiled"]) is False:
                self.logger.debug("%s fails the minimum requirements for AS events", transcript.id)
                to_be_added = False

        if to_be_added is True:
            is_alternative, ccode, comparison = self.is_alternative_splicing(transcript)
            if is_alternative is False:
                self.logger.debug("%s not added because it is not a \
                valid splicing isoform. Ccode: %s",
                                  transcript.id, ccode)
                to_be_added = False
            else:
                transcript.attributes["ccode"] = ccode
                self.logger.debug("%s is a valid splicing isoform; Ccode: %s", transcript.id, ccode)

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

    def is_putative_fragment(self):

        """This method will use the expression in the "not_fragmentary" section
        of the configuration to determine whether it is itself a putative fragment."""

        self.json_conf["not_fragmentary"]["compiled"] = compile(
            self.json_conf["not_fragmentary"]["expression"], "<json>",
            "eval")

        evaluated = dict()
        for key in self.json_conf["not_fragmentary"]["parameters"]:
            value = getattr(self.primary_transcript,
                            self.json_conf["not_fragmentary"]["parameters"][key]["name"])
            evaluated[key] = self.evaluate(
                value,
                self.json_conf["not_fragmentary"]["parameters"][key])
        if eval(self.json_conf["not_fragmentary"]["compiled"]) is True:
            self.logger.debug("%s cannot be a fragment according to the definitions, keeping it",
                              self.id)
            return False
        else:
            self.logger.debug(
                "%s could be a fragment according to the definitions, tagging it for analysis",
                self.id)
            return True

    def other_is_fragment(self,
                          other):
        """
        :param other: another Locus to compare against
        :type other: Locus

        This function checks whether another *monoexonic* Locus
        *on the opposite strand* is a fragment,by checking its classification
        according to Assigner.compare.
        """

        if not isinstance(self, type(other)):
            raise TypeError("I can compare only loci.")

        if other.primary_transcript_id == self.primary_transcript_id:
            self.logger.debug("Self-comparisons are not allowed!")
            return False, None

        self.logger.debug("Comparing %s with %s",
                          self.primary_transcript_id,
                          other.primary_transcript_id)

        result, _ = Assigner.compare(other.primary_transcript,
                                     self.primary_transcript,
                                     strict_strandedness=True)
        max_distance = self.json_conf["pick"]["fragments"]["max_distance"]
        self.logger.debug("Comparison between {0} (strand {3}) and {1}: class code \"{2}\"".format(
            self.primary_transcript.id,
            other.primary_transcript.id,
            result.ccode[0],
            other.strand))
        if (result.ccode[0] in self.json_conf["pick"]["fragments"]["valid_class_codes"] and
                    result.distance[0] <= max_distance):
            self.logger.debug("{0} is a fragment (ccode {1})".format(
                other.primary_transcript.id, result.ccode[0]))
            return True, result

        return False, None

    def set_json_conf(self, jconf: dict):
        """
        Setter for the configuration dictionary.
        :param jconf:
        :type jconf: dict
        """
        if not isinstance(jconf, dict):
            raise TypeError("Invalid configuration of type {0}".format(type(jconf)))
        self.json_conf = jconf

    # def get_metrics(self):
    #
    #     """Quick wrapper to calculate the metrics for all the transcripts."""
    #
    #     if self.metrics_calculated is True:
    #         return
    #
    #     assert len(self._cds_introntree) == len(self.combined_cds_introns)
    #
    #     for tid in sorted(self.transcripts):
    #         self.calculate_metrics(tid)
    #
    #     self.logger.debug("Finished to calculate the metrics for %s", self.id)
    #
    #     self.metrics_calculated = True
    #     return

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

        super().calculate_scores()

        for index, item in enumerate(reversed(self.metric_lines_store)):
            if item["tid"] in self.__orf_doubles:
                del self.metric_lines_store[index]
            else:
                continue

        for doubled in self.__orf_doubles:
            for partial in self.__orf_doubles[doubled]:
                if partial in self.transcripts:
                    del self.transcripts[partial]
                if partial in self.scores:
                    del self.scores[partial]

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
            if tid in self._not_passing:
                row["score"] = 0
            else:
                row["score"] = round(self.scores[tid]["score"], 2)
            calculate_total = (self.regressor is None)
            for key in score_keys:
                if calculate_total:
                    assert self.scores[tid][key] != "NA" and self.scores[tid][key] is not None
                    row[key] = round(self.scores[tid][key], 2)

            if calculate_total is True:
                score_sum = sum(row[key] for key in score_keys)
                if tid not in self._not_passing and self.scores[tid]["score"] > 0:
                    assert round(score_sum, 2) == round(self.scores[tid]["score"], 2), (
                        score_sum,
                        self.transcripts[tid].score,
                        tid)
                else:
                    assert self.scores[tid]["score"] == 0

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

        if self.json_conf["pick"]["clustering"]["cds_only"] is True:
            main_result, _ = Assigner.compare(other._selected_orf_transcript,
                                              self.primary_transcript._selected_orf_transcript)
            enough_overlap, overlap_reason = self._evaluate_transcript_overlap(
                other._selected_orf_transcript,
                self.primary_transcript._selected_orf_transcript,
                min_cdna_overlap=self.json_conf["pick"]["alternative_splicing"]["min_cdna_overlap"],
                min_cds_overlap=self.json_conf["pick"]["alternative_splicing"]["min_cds_overlap"],
                comparison=main_result,
                is_internal_orf=True)
        else:
            main_result, _ = Assigner.compare(other,
                                              self.primary_transcript)
            enough_overlap, overlap_reason = self._evaluate_transcript_overlap(
                other,
                self.primary_transcript,
                min_cdna_overlap=self.json_conf["pick"]["alternative_splicing"]["min_cdna_overlap"],
                min_cds_overlap=self.json_conf["pick"]["alternative_splicing"]["min_cds_overlap"],
                comparison=main_result,
                is_internal_orf=False)

        main_ccode = main_result.ccode[0]

        if main_ccode not in valid_ccodes:
            self.logger.debug("%s is not a valid splicing isoform. Ccode: %s",
                              other.id,
                              main_result.ccode[0])
            is_valid = False
        elif not enough_overlap:
            self.logger.debug("%s is not a valid splicing isoform. Reason: %s",
                              other.id, overlap_reason)
            is_valid = False

        if is_valid:
            for tid in iter(tid for tid in self.transcripts if
                            tid not in (self.primary_transcript_id, other.id)):
                candidate = self.transcripts[tid]
                if self.json_conf["pick"]["clustering"]["cds_only"] is True:
                    result, _ = Assigner.compare(
                        other._selected_orf_transcript,
                        candidate._selected_orf_transcript)
                else:
                    result, _ = Assigner.compare(other, candidate)
                if result.ccode[0] in redundant_ccodes:
                    self.logger.debug("%s is a redundant isoform of %s (ccode %s)",
                                      other.id, candidate.id, result.ccode[0])
                    is_valid = False
                    break

        return is_valid, main_ccode, main_result

    def pad_transcripts(self):

        """
        """

        try:
            self.fai = pyfaidx.Fasta(self.json_conf["reference"]["genome"])
        except KeyError:
            raise KeyError(self.json_conf.keys())

        five_graph = self.define_graph(self.transcripts, self.__share_extreme, three_prime=False)
        three_graph = self.define_graph(self.transcripts, self.__share_extreme, three_prime=True)

        five_comm = deque(sorted(self.find_communities(five_graph),
                              key=lambda clique: min(self[_].start for _ in clique)))
        three_comm = deque(sorted(self.find_cliques(three_graph),
                           key=lambda clique: max(self[_].end for _ in clique),
                           reverse=True))

        __to_modify = self._find_communities_boundaries(five_comm, three_comm)

        # Now we can do the proper modification
        for tid in __to_modify:
            new_transcript = expand_transcript(self[tid].deepcopy(),
                                               __to_modify[tid][0],
                                               __to_modify[tid][1],
                                               self.fai,
                                               self.logger)
            self.transcripts[tid] = new_transcript

        self.exons = set()
        for tid in self:
            self.exons.update(self[tid].exons)

    def _find_communities_boundaries(self, five_comm, three_comm):

        five_found = set()

        __to_modify = dict()

        while len(five_comm) > 0:

            comm = five_comm.popleft()
            comm = deque(sorted(list(set.difference(set(comm), five_found)),
                         key=lambda internal_tid: self[internal_tid].start))
            if len(comm) == 1:
                continue
            first = comm.popleft()
            five_found.add(first)
            comm_start = self[first].start

            for tid in comm:
                if ((self[tid].start - comm_start + 1) <
                        self.json_conf["pick"]["alternative_splicing"]["ts_distance"] and
                        len([_ for _ in self.splices if comm_start <= _ <= self[tid].start]) <
                        self.json_conf["pick"]["alternative_splicing"]["ts_max_splices"] and
                        self[tid].start > comm_start):
                    __to_modify[tid] = [comm_start, False]
                    five_found.add(tid)
                else:
                    continue
            comm = deque([_ for _ in comm if _ not in five_found])

            if comm:
                five_comm.appendleft(comm)

        # Then do the 3' end

        three_found = set()

        while len(three_comm) > 0:

            comm = three_comm.popleft()
            comm = deque(sorted(list(set.difference(set(comm), three_found)),
                         key=lambda internal_tid: self[internal_tid].end, reverse=True))
            if len(comm) == 1:
                continue
            first = comm.popleft()
            three_found.add(first)
            comm_end = self[first].end
            for tid in comm:
                if ((self[tid].end - comm_end + 1) <
                        self.json_conf["pick"]["alternative_splicing"]["ts_distance"] and
                        len([_ for _ in self.splices if self[tid].end <= _ <= comm_end]) <
                        self.json_conf["pick"]["alternative_splicing"]["ts_max_splices"] and
                        self[tid].end < comm_end):

                    if tid in __to_modify:
                        __to_modify[tid][1] = comm_end
                    else:
                        __to_modify[tid] = [False, comm_end]

                    three_found.add(tid)
                else:
                    continue
            comm = deque([_ for _ in comm if _ not in three_found])
            if comm:
                three_comm.appendleft(comm)

        return __to_modify

    def __share_extreme(self, first, second, three_prime=False):

        """
        :param first:
        :param second:
        :return:
        """

        if not three_prime:
            return (overlap(first.exons[0], second.exons[0]) > 0 and
                    max(first.start, second.start) + 1 - min(first.start, second.start) < self.json_conf[
                        "pick"]["alternative_splicing"]["ts_distance"])
        else:
            return (overlap(first.exons[-1], second.exons[-1]) > 0 and
                    max(first.end, second.end) + 1 - min(first.end, second.end) < self.json_conf[
                        "pick"]["alternative_splicing"]["ts_distance"])

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

        mapper = {old_primary: primary_id}

        for counter, tid in enumerate(order):
            counter += 2
            self.transcripts[tid].attributes["Alias"] = tid
            new_id = "{0}.{1}".format(string, counter)
            self.transcripts[tid].id = new_id
            self.transcripts[new_id] = self.transcripts.pop(tid)
            mapper[tid] = new_id

        if self.scores_calculated is True:
            for tid in mapper:
                self.scores[mapper[tid]] = self.scores.pop(tid)
        if self.metrics_calculated is True:
            for index in range(len(self.metric_lines_store)):
                self.metric_lines_store[index]["tid"] = mapper[self.metric_lines_store[index]["tid"]]
                self.metric_lines_store[index]["parent"] = self.id

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

    @property
    def purge(self):

        """Overloading of the base property. Loci should never purge."""

        return False

    @property
    def _finalized(self):
        return self.__finalized

    @_finalized.setter
    def _finalized(self, value):
        if not isinstance(value, bool):
            raise ValueError(value)
        self.__finalized = value


def expand_transcript(transcript, new_start, new_end, fai, logger):

    # First get the ORFs
    transcript.logger = logger
    if transcript.combined_cds_length > 0:
        internal_orfs = list(transcript.get_internal_orf_beds())
    else:
        internal_orfs = []
    # Remove the CDS and unfinalize
    transcript.strip_cds()
    transcript.unfinalize()

    upstream = 0
    downstream = 0
    if new_start:
        __new_exon = (new_start, transcript.exons[0][1])
        upstream = transcript.start - new_start
        transcript.start = new_start
        transcript.remove_exon(transcript.exons[0])
        transcript.add_exon(__new_exon)
        transcript.exons = sorted(transcript.exons)
    if new_end:
        __new_exon = (transcript.exons[-1][0], new_end)
        downstream = new_end - transcript.end
        transcript.end = new_end
        transcript.remove_exon(transcript.exons[-1])
        transcript.add_exon(__new_exon)
        transcript.exons = sorted(transcript.exons)

    transcript.finalize()
    # Now for the difficult part
    if internal_orfs and (new_start or new_end):
        logger.debug("Enlarging the ORFs for TID %s (%s)",
                       transcript.id, (new_start, new_end))
        new_orfs = []
        seq = "".join(
            TranscriptChecker(
                transcript,
                fai[transcript.chrom][transcript.start-1:transcript.end].seq
            ).fasta.split("\n")[1:])
        assert len(seq) == transcript.cdna_length, (len(seq), transcript.cdna_length, transcript.exons)
        for orf in internal_orfs:
            logger.debug("Old ORF: %s", str(orf))
            try:
                orf.expand(seq, upstream, downstream)
            except AssertionError as err:
                logger.error(err)
                logger.error("%s, %s, %s, %s, %s, %s",
                             new_start,
                             new_end,
                             upstream,
                             downstream,
                             transcript.exons,
                             transcript.cdna_length)
                raise AssertionError(err)
            logger.debug("New ORF: %s", str(orf))
            new_orfs.append(orf)
        transcript.load_orfs(new_orfs)

    # Now finalize again
    transcript.finalize()
    return transcript
