# coding: utf-8

"""
This module defines the last object to be created during the picking,
i.e. the locus.
"""

import itertools
from mikado_lib.loci_objects.transcript import Transcript
from mikado_lib.scales.assigner import Assigner
from mikado_lib.loci_objects.monosublocus import Monosublocus
from mikado_lib.loci_objects.abstractlocus import Abstractlocus


class Locus(Monosublocus, Abstractlocus):
    """Class that defines the final loci.
    It is a child of monosublocus, but it also has the possibility of adding
    additional transcripts if they are valid splicing isoforms.
    """

    __name__ = "locus"

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

    def __str__(self, print_cds=True) -> str:

        self.feature = self.__name__
        return super().__str__(print_cds=print_cds)

    def add_transcript_to_locus(self, transcript: Transcript, **kwargs):
        """Implementation of the add_transcript_to_locus method.
        Before a transcript is added, the class checks that it is a valid splicing isoform
        and that we have not exceeded already the maximum number of isoforms for the Locus.

        :param transcript: the candidate transcript
        :type transcript: Transcript

        :param kwargs: optional keyword arguments are ignored.
        """

        if len(self.transcripts) >= self.json_conf["alternative_splicing"]["max_isoforms"]:
            self.logger.debug("{0} not added because the Locus has already too many transcripts.".format(transcript.id))
            return
        if not self.is_alternative_splicing(transcript):
            self.logger.debug("{0} not added because it is not a valid splicing isoform.".format(transcript.id))
            return
        if transcript.combined_utr_length > self.json_conf["alternative_splicing"]["max_utr_length"]:
            self.logger.debug("{0} not added because it has too much UTR ({1)}.".format(transcript.id,
                                                                                        transcript.combined_utr_length))
            return
        if transcript.five_utr_length > self.json_conf["alternative_splicing"]["max_fiveutr_length"]:
            self.logger.debug("{0} not added because it has too much 5'UTR ({1)}.".format(transcript.id,
                                                                                          transcript.five_utr_length))
            return
        if transcript.three_utr_length > self.json_conf["alternative_splicing"]["max_threeutr_length"]:
            self.logger.debug("{0} not added because it has too much 5'UTR ({1)}.".format(transcript.id,
                                                                                          transcript.three_utr_length))
            return

        if self.json_conf["alternative_splicing"]["keep_retained_introns"] is False:
            self.find_retained_introns(transcript)
            if transcript.retained_intron_num > 0:
                l = transcript.retained_intron_num
                self.logger.debug("{0} not added because it has {1} retained introns.".format(transcript.id,
                                                                                              l))
                return
        if self.json_conf["alternative_splicing"][
                "min_cds_overlap"] > 0 and self.primary_transcript.combined_cds_length > 0:
            tr_nucls = set(itertools.chain(*[range(x[0], x[1] + 1) for x in transcript.combined_cds]))
            primary_nucls = set(itertools.chain(*[range(x[0], x[1] + 1) for x in self.primary_transcript.combined_cds]))
            nucl_overlap = len(set.intersection(primary_nucls, tr_nucls))
            ol = nucl_overlap / self.primary_transcript.combined_cds_length
            if ol < self.json_conf["alternative_splicing"]["min_cds_overlap"]:
                self.logger.debug(
                    "{0} not added because its CDS overlap with the primary CDS is too low ({1:.2f}%).".format(
                        transcript.id,
                        ol * 100))
                return

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

        :param minimal_cds_length: Minimal CDS length to consider a Locus as non-fragment no matter the ccode.
        :type minimal_cds_length: int


        This function checks whether another *monoexonic* Locus on the opposite strand* is a fragment,
        by checking its classification according to Assigner.compare.
        Briefly, a transcript is classified as fragment if it follows the following criteria:
        
            - it is monoexonic
            - it has a combined_cds_length inferior to maximal_cds
            - it is classified as x,i,P
        """

        if type(self) != type(other):
            raise TypeError("I can compare only loci.")

        self.logger.debug("Comparing {0} with {1}".format(self.primary_transcript_id, other.primary_transcript_id))

        if other.primary_transcript.combined_cds_length > minimal_cds_length:
            self.logger.debug("{0} has a CDS of {1}, not a fragment by definition".format(
                other.primary_transcript_id, other.primary_transcript.combined_cds_length))
            return False

        result, _ = Assigner.compare(other.primary_transcript, self.primary_transcript)
        # Exclude anything which is completely contained within an intron,
        # or is a monoexonic fragment overlapping/in the neighborhood
        self.logger.debug("Comparison between {0} and {1}: class code \"{2}\"".format(
            self.primary_transcript.id,
            other.primary_transcript.id,
            result.ccode[0]))
        if result.ccode[0] in ("i", "P", "p", "x"):
            return True
        # if result.ccode[0] == "x":
        #     # If the transcript bridges an intron completely, it probably is not a a fragment.
        #     # Otherwise, return True.
        #     ostart, oend = other.primary_transcript.start, other.primary_transcript.end
        #     if sum(1 for x in filter(lambda intron:
        #                              self.overlap(intron, (ostart, oend)) >= (intron[1]-intron[0]+1),
        #                              self.primary_transcript.introns)) >= 1:
        #         return False
        #     else:
        #         return True

        return False

    def set_json_conf(self, jconf: dict):
        """
        Setter for the configuration dictionary.
        :param jconf:
        :type jconf: dict
        """
        if type(jconf) is not dict:
            raise TypeError("Invalid configuration of type {0}".format(type(jconf)))
        self.json_conf = jconf

    def is_alternative_splicing(self, other):

        """This function defines whether another transcript could be a putative alternative splice variant.
        To do so, it compares the candidate against all transcripts in the Locus, and calculates
        the class code using scales.Assigner.compare.
        If all the matches are "n" or "j", the transcript is considered as an AS event.

        :param other: another transcript to compare against
        :type other: mikado_lib.loci_objects.transcript.Transcript

        """

        if other.id == self.primary_transcript_id:
            return False
        if self.overlap((other.start, other.end), (self.start, self.end)) < 0:
            return False

        if other.strand != other.strand:
            return False
        if other.retained_intron_num > 0:
            return False

        for tid in self.transcripts:
            result, _ = Assigner.compare(other, self.transcripts[tid])
            self.logger.debug("{0} vs. {1}: {2}".format(tid, other.id, result.ccode[0]))
            if result.ccode[0] not in ("j", "n") or \
                    (self.transcripts[tid].monoexonic is True and result.ccode[0] in ("o", "O")):
                self.logger.debug("{0} is not a valid splicing isoform. Ccode: {1}".format(other.id, result.ccode[0]))
                return False
            if result.n_f1 == 0 or ((self.transcripts[tid].monoexonic is False and result.j_f1 == 0) or
                                    self.transcripts[tid].monoexonic is True):
                self.logger.debug("{0} is not a valid splicing isoform. N_f1: {1}; J_f1: {2}".format(other.id,
                                                                                                     result.n_f1,
                                                                                                     result.j_f1
                                                                                                     ))
                return False

        return True

    @property
    def id(self):
        """
        Override of the abstractlocus method.
        :rtype str
        """
        myid = Abstractlocus.id.fget(self)  # @UndefinedVariable
        if self.counter > 0:
            myid = "{0}.{1}".format(myid, self.counter)
        return myid

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
        if not type(val) is bool:
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
