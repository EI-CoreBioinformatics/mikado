# coding: utf-8

"""
This module defines the last object to be created during the picking,
i.e. the locus.
"""

import itertools
from .transcript import Transcript
from ..scales.assigner import Assigner
from .monosublocus import Monosublocus
from .abstractlocus import Abstractlocus


class Locus(Monosublocus, Abstractlocus):
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
        super().__init__(transcript, logger=logger)
        self.logger.debug("Created Locus object with {0}".format(transcript.id))
        self.primary_transcript_id = transcript.id

        self.attributes["is_fragment"] = False
        self.__id = None

    def __str__(self, print_cds=True, source_in_name=True) -> str:

        self.feature = self.__name__
        # Hacky fix to make sure that the primary transcript has the attribute
        # Set to True in any case.
        self.primary_transcript.attributes["primary"] = True
        # BF, really just a hack.
        for transcript in self.transcripts:
            if transcript == self.primary_transcript_id:
                continue
            self.transcripts[transcript].attributes["primary"] = False

        return super().__str__(print_cds=print_cds, source_in_name=source_in_name)

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
            is_alternative, ccode = self.is_alternative_splicing(transcript)
            if is_alternative is False:
                self.logger.debug("%s not added because it is not a \
                valid splicing isoform. Ccode: %s",
                                  transcript.id, ccode)
                to_be_added = False
            else:
                transcript.attributes["ccode"] = ccode
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
        if to_be_added and self.json_conf["pick"]["alternative_splicing"]["min_cdna_overlap"] > 0:
            tr_nucls = set(itertools.chain(*[range(x[0], x[1] + 1) for x in transcript.exons]))
            primary_nucls = set(itertools.chain(*[range(x[0], x[1] + 1) for x in self.primary_transcript.exons]))
            nucl_overlap = len(set.intersection(primary_nucls, tr_nucls))
            overlap = nucl_overlap / len(self.primary_transcript)
            if overlap < self.json_conf["pick"]["alternative_splicing"]["min_cdna_overlap"]:
                self.logger.debug(
                    "%s not added because its CDNA overlap is too low (%f%%).",
                    transcript.id,
                    round(overlap * 100, 2))
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

    def is_intersecting(self):
        """Not implemented: this function makes no sense for a single-transcript container."""
        raise NotImplementedError("""Loci do not use this method, but rather
        assess whether a transcript is a splicing isoform or not.""")

    def other_is_fragment(self, other, minimal_cds_length=0):
        """
        :param other: another Locus to compare against
        :type other: Locus

        :param minimal_cds_length: Minimal CDS length to consider
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
        elif other.strand is None and result.ccode[0] in ("_", "o", "e"):
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

    def is_alternative_splicing(self, other):

        """This function defines whether another transcript could be a
        putative alternative splice variant of the primary Locus
        transcript.
        To do so, it compares the candidate against all transcripts in the Locus, and calculates
        the class code using scales.Assigner.compare.
        If all the matches are "n" or "j", the transcript is considered as an AS event.

        :param other: another transcript to compare against
        :type other: mikado_lib.loci_objects.transcript.Transcript

        """

        is_valid = True
        # main_ccode = None

        valid_ccodes = self.json_conf["pick"]["alternative_splicing"]["valid_ccodes"]

        ccodes = []
        main_result, _ = Assigner.compare(other, self.primary_transcript)
        main_ccode = main_result.ccode[0]
        results = [main_result]

        if main_ccode not in valid_ccodes:
            self.logger.debug("%s is not a valid splicing isoform. Ccode: %s",
                              other.id,
                              main_result.ccode[0])
            is_valid = False
        if is_valid:
            for tid in iter(tid for tid in self.transcripts if
                            tid != self.primary_transcript_id):
                result, _ = Assigner.compare(other, self.transcripts[tid])
                results.append(result)
                ccodes.append(result.ccode[0])
                self.logger.debug("Comparing secondary transcripts %s vs %s. Ccode: %s",
                                  tid, other.id, result.ccode[0])
            if any(_ not in valid_ccodes for _ in ccodes):
                not_valid = set(_ for _ in ccodes if _ not in valid_ccodes)
                self.logger.debug(
                    "%s not a valid AS event. Non valid ccodes: %s",
                    other.id,
                    ", ".join(list(not_valid)))
                is_valid = False

            # if (("_" in ccodes and "_" not in valid_ccodes) or
            #         ("=" in ccodes and "_" not in valid_ccodes)):
            #     self.logger.debug("%s is a redundant valid splicing isoform. Ccode: %s",
            #                       other.id,
            #                       main_result.ccode[0])
            #     is_valid = False
        return is_valid, main_ccode

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
