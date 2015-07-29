# coding: utf-8

"""
This module defines a child of the Transcript class, which is used
to verify that e.g. the assigned strand is correct.
"""

from mikado_lib.loci_objects.transcript import Transcript
from mikado_lib.parsers.GFF import GffLine
from mikado_lib.parsers.GTF import GtfLine
from mikado_lib.exceptions import IncorrectStrandError


class TranscriptChecker(Transcript):
    """This is a subclass of the generic transcript class. Its purpose is to compare
    the information of the transcript instance with the information contained in a
    genomic FASTA file, to verify some of the information.
    At the moment, the class implements only a check on the strandedness made by extracting
    the FASTA sequence of the splice sites and verifying that they are concordant with the annotated strand.
    Keyword arguments:
        - strand_specific: if set, monoexonic transcripts are not set to "unknown" strand.
        - lenient: boolean. If set to True, a transcript with mixed splices will not throw an exception
        but rather will just report the number of splices supporting each strand.
    """

    def __init__(self, gffline, fasta_index, strand_specific=False, lenient=False):

        """
        Constructor method. It inherits from Transcript, with some modifications.

        :param gffline: annotation line to begin the construction.
        :type gffline: GffLine | GtfLine

        :param fasta_index: a SeqIO indexed FASTA file

        :param strand_specific: flag. If set, transcripts will not have their strand changed.
        :type strand_specific: bool

        :param lenient: boolean flag. If set, incorrect transcripts will be flagged rather than be discarded.
        :type lenient: bool
        """
        self.__strand_specific = False
        self.mixed_attribute = ""

        if fasta_index is None:
            raise ValueError()
        super().__init__(gffline)
        self.original_strand = gffline.strand
        assert self.original_strand == self.strand
        self.parent = gffline.parent
        self.fasta_index = fasta_index
        self.strand_specific = strand_specific
        self.checked = False
        self.lenient = lenient
        self.mixed_splices = False
        self.reversed = False

    @property
    def strand_specific(self):
        """
        Flag, set from the constructor. If True, transcript will not have their strand changed.
        :rtype: bool
        """
        return self.__strand_specific

    @strand_specific.setter
    def strand_specific(self, value):
        """
        Setter for strand_specific. Only boolean values are considered as valid.
        :param value: flag
        :type value: bool

        """
        if type(value) is not bool:
            raise TypeError("Invalid value for boolean property: {0}".format(value))
        self.__strand_specific = value

    def __str__(self, print_cds=True, to_gtf=False):

        self.check_strand()
        if self.mixed_splices is True:
            self.attributes["mixed_splices"] = self.mixed_attribute

        return super().__str__(print_cds=print_cds, to_gtf=to_gtf)

    def check_strand(self):
        """
        This method will check that the transcript instance has all the splice sites on one strand,
        or at most with non-canonical (therefore unknowable) splice junctions.
        If the transcript is monoexonic and strand_specific is set to False, the strand of the transcript
        will be set to None.

        The finalize method is called preliminarly before any operation.
        """
        self.finalize()
        if self.checked is True:
            return

        canonical_splices = [
            ("GT", "AG"),
            ("GC", "AG"),
            ("AT", "AC")
        ]

        if self.strand_specific is False and self.monoexonic is True:
            self.strand = None
            return

        elif self.monoexonic is False:
            canonical_counter = dict()
            for strand in ("+", "-", None):
                canonical_counter[strand] = 0

            for intron in self.introns:
                splice_donor = self.fasta_index[self.chrom][intron[0] - 1:intron[0] + 1]
                splice_acceptor = self.fasta_index[self.chrom][intron[1] - 2:intron[1]]
                if self.strand == "-":
                    rsa, rsd = splice_acceptor.reverse_complement(), splice_donor.reverse_complement()
                    splice_donor, splice_acceptor = rsa, rsd
                if (splice_donor.seq, splice_acceptor.seq) in canonical_splices:
                    canonical_counter["+"] += 1
                else:
                    rsa, rsd = splice_acceptor.reverse_complement(), splice_donor.reverse_complement()
                    splice_donor, splice_acceptor = rsa, rsd
                    if (splice_donor.seq, splice_acceptor.seq) in canonical_splices:
                        canonical_counter["-"] += 1
                    else:
                        canonical_counter[None] += 1

            if canonical_counter[None] == len(self.introns):
                if self.lenient is False:
                    raise IncorrectStrandError("No correct strand found for {0}".format(self.id))

            elif canonical_counter["+"] > 0 and canonical_counter["-"] > 0:
                if self.lenient is False:
                    raise IncorrectStrandError(
                        "Transcript {0} has {1} positive and {2} negative splice junctions. Aborting.".format(
                            self.id,
                            canonical_counter["+"],
                            canonical_counter["-"]
                        )
                    )
                else:
                    self.mixed_splices = True

                if canonical_counter["+"] >= canonical_counter["-"]:
                    self.mixed_attribute = "{0}concordant,{1}discordant".format(canonical_counter["+"],
                                                                                canonical_counter["-"])
                else:
                    self.reverse_strand()
                    self.mixed_attribute = "{0}concordant,{1}discordant".format(canonical_counter["-"],
                                                                                canonical_counter["+"])

            elif canonical_counter["-"] > 0:
                self.reverse_strand()
                self.reversed = True

        self.checked = True
