# coding: utf-8

"""
This module defines a child of the Transcript class, which is used
to verify that e.g. the assigned strand is correct.
"""

from .transcript import Transcript
from ..exceptions import IncorrectStrandError
from collections import Counter
from itertools import zip_longest


# pylint: disable=too-many-instance-attributes
class TranscriptChecker(Transcript):
    """This is a subclass of the generic transcript class. Its purpose is to compare
    the information of the transcript instance with the information contained in a
    genomic FASTA file, to verify some of the information.
    At the moment, the class implements only a check on the strandedness made by extracting
    the FASTA sequence of the splice sites and verifying that they are concordant
    with the annotated strand.
    Keyword arguments:
        - strand_specific: if set, monoexonic transcripts are not set to "unknown" strand.
        - lenient: boolean. If set to True, a transcript with mixed splices will
        not throw an exception but rather will just report the number
        of splices supporting each strand.
    """

    __translation_table = str.maketrans("ACGT", "TGCA")

    # The arguments are all needed here.
    # pylint: disable=too-many-arguments
    def __init__(self, gffline, seq,
                 strand_specific=False, lenient=False,
                 canonical_splices=(("GT", "AG"), ("GC", "AG"), ("AT", "AC")),
                 logger=None):

        """
        Constructor method. It inherits from Transcript, with some modifications.

        :param gffline: annotation line to begin the construction.
        :type gffline: (Mikado.parsers.gfannotation.GFAnnotation|Transcript)

        :param seq: a SeqIO indexed FASTA file

        :param strand_specific: flag. If set, transcripts will not have their strand changed.
        :type strand_specific: bool

        :param lenient: boolean flag. If set, incorrect transcripts will be
        flagged rather than discarded.
        :type lenient: bool
        """
        self.__strand_specific = False
        self.mixed_attribute = ""

        if seq is None:
            raise ValueError("No sequence provided!")
        if isinstance(gffline, (Transcript, TranscriptChecker)):
            self.__dict__.update(gffline.__dict__)
        else:
            super().__init__(gffline)
        self.original_strand = gffline.strand
        assert self.original_strand == self.strand
        self.attributes.update(gffline.attributes)
        self.parent = gffline.parent
        self.fasta_seq = seq
        self.strand_specific = strand_specific
        self.checked = False
        self.lenient = lenient
        self.mixed_splices = False
        self.reversed = False
        self.canonical_splices = []
        if not isinstance(canonical_splices, (tuple, list)):
            raise ValueError("Canonical splices should be provided as lists or tuples")

        for canonical_splice in canonical_splices:
            self.canonical_splices.append((str(canonical_splice[0]),
                                           str(canonical_splice[1])))

        self.canonical_junctions = []
        self.logger = logger
    # pylint: enable=too-many-arguments

    @property
    def translation_table(self):
        """
        Returns the table used to reverse complement FASTA strings.
        """
        return self.__translation_table

    @classmethod
    def get_translation_table(cls):
        """Class method to access the translation table."""

        return cls.__translation_table

    @classmethod
    def rev_complement(cls, string):

        """
        Quick method to perform the reverse complement of a given string,
        using the class translation table.

        :param string: the sequence to be rev-complented
        :type string: str
        """

        return "".join(x for x in reversed(
            string.translate(cls.get_translation_table())
        ))

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
        if not isinstance(value, bool):
            raise TypeError("Invalid value for boolean property: {0}".format(value))
        self.__strand_specific = value

    def __str__(self, print_cds=True, to_gtf=False, with_introns=False):

        self.check_strand()
        if self.mixed_splices is True:
            self.attributes["mixed_splices"] = self.mixed_attribute

        return super().__str__(print_cds=print_cds, to_gtf=to_gtf)

    def format(self, format_name, with_introns=False, with_cds=True):

        self.check_strand()
        if self.mixed_splices is True:
            self.attributes["mixed_splices"] = self.mixed_attribute

        return super().format(format_name, with_cds=with_cds, with_introns=with_introns)

    def check_strand(self):
        """
        This method will check that the transcript instance has all the splice sites
        on one strand,or at most with non-canonical (therefore unknowable) splice junctions.
        If the transcript is monoexonic and strand_specific is set to False,
        the strand of the transcript will be set to None.

        The finalize method is called preliminarly before any operation.
        """
        self.finalize()
        assert self.exons[0][0] - self.start == 0
        assert self.exons[-1][1] - self.start + 1 == len(self.fasta_seq)
        if self.checked is True:
            return

        if self.strand_specific is False and self.monoexonic is True:
            self.strand = None

        elif self.monoexonic is False:
            canonical_counter = Counter()
            canonical_index = dict()

            for pos, intron in enumerate(self.introns):
                canonical_index[pos+1] = self._check_intron(intron)

            canonical_counter.update(canonical_index.values())

            if canonical_counter[None] == len(self.introns):
                self.logger.debug("Transcript %s only has non-canonical splices", self.id)
                if self.lenient is False:
                    raise IncorrectStrandError("No correct strand found for {0}".format(
                        self.id))

            elif canonical_counter["+"] > 0 and canonical_counter["-"] > 0:
                if self.lenient is False:
                    err_messg = """Transcript {0} has {1} positive and {2} negative splice junctions. \
    Aborting.""".format(self.id,
                        canonical_counter["+"],
                        canonical_counter["-"])
                    self.logger.warning(err_messg)
                    raise IncorrectStrandError(err_messg)

                else:
                    self.mixed_splices = True
                    self.logger.warning("Transcript %s has %d positive and %s negative splice junctions",
                                        self.id,
                                        canonical_counter["+"],
                                        canonical_counter["-"])

                    if canonical_counter["+"] >= canonical_counter["-"] or self.strand_specific is True:
                        self.mixed_attribute = "{0}concordant,{1}discordant".format(
                            canonical_counter["+"],
                            canonical_counter["-"])
                    else:
                        self.reverse_strand()
                        self.reversed = True
                        self.mixed_attribute = "{0}concordant,{1}discordant".format(
                            canonical_counter["-"],
                            canonical_counter["+"])

            elif canonical_counter["-"] > 0 and self.strand_specific is False:
                self.reverse_strand()
                self.reversed = True
            elif canonical_counter["-"] > 0 and self.strand_specific is True:
                self.logger.warning(
                    "Transcript %s has been assigned to the wrong strand, tagging it but leaving it on this strand.",
                    self.id,
                    canonical_counter["+"],
                    canonical_counter["-"])
                self.attributes["canonical_on_reverse_strand"] = True
            elif self.strand_specific is False:
                assert canonical_counter["+"] + canonical_counter[None] == len(self.introns)

            # self.attributes["canonical_splices"] = ",".join(
            #     str(k) for k in canonical_index.keys() if
            #     canonical_index[k] in ("+", "-"))

            if canonical_counter["+"] >= canonical_counter["-"] or self.strand_specific is True:
                strand = "+"
            else:
                strand = "-"

            self.attributes["canonical_number"] = canonical_counter[strand]
            self.attributes["canonical_proportion"] = canonical_counter[strand] / len(self.introns)

            self.canonical_junctions = [_ for _ in canonical_index if
                                        canonical_index[_] == strand]
            self.attributes["canonical_junctions"] = ",".join([str(_) for _
                                                               in self.canonical_junctions])

        self.checked = True
        return

    def _check_intron(self, intron):

        """
        Private method that checks whether an intron has canonical splice sites
        or not.
        :param intron: the intron tuple (int,int) in 1-base offset
        [("AG","GT")]
        :return: strand of the intron (None | "+" | "-")
        """

        splice_donor = self.fasta_seq[intron[0] - self.start:intron[0]-self.start + 2]
        assert len(splice_donor) == 2
        splice_acceptor = self.fasta_seq[intron[1] - 1 - self.start:intron[1] - self.start + 1]

        if not isinstance(splice_donor, str):
            splice_donor = str(splice_donor.seq)
            splice_acceptor = str(splice_acceptor.seq)

        assert isinstance(splice_acceptor, str)

        # splice_donor = self.fasta_index[self.chrom][intron[0] - 1:intron[0] + 1]
        # splice_acceptor = self.fasta_index[self.chrom][intron[1] - 2:intron[1]]
        if self.strand == "-":
            splice_donor, splice_acceptor = (self.rev_complement(splice_acceptor),
                                             self.rev_complement(splice_donor))
        if (splice_donor, splice_acceptor) in self.canonical_splices:
            strand = "+"
        else:
            splice_donor, splice_acceptor = (self.rev_complement(splice_acceptor),
                                             self.rev_complement(splice_donor))
            if (splice_donor, splice_acceptor) in self.canonical_splices:
                strand = "-"
            else:
                strand = None
        return strand

    @property
    def fasta(self):
        """
        This property calculates and returns the FASTA sequence associated with
        the transcript instance.
        The FASTA sequence itself will be formatted to be in lines with 60 characters
        (the standard).
        :return:
        """

        self.check_strand()
        fasta = [">{0}".format(self.id)]
        sequence = ''

        for exon in self.exons:
            start = exon[0] - self.start
            end = exon[1] + 1 - self.start
            _ = self.fasta_seq[start:end]
            sequence += _

        # pylint: disable=no-member
        if not isinstance(sequence, str):
            if hasattr(sequence, "seq"):
                sequence = str(sequence.seq)
            else:
                raise TypeError("Invalid object for sequence: {0} ({1})".format(
                    type(sequence),
                    repr(sequence)
                ))
        # pylint: enable=no-member

        if self.strand == "-":
            sequence = self.rev_complement(sequence)
        # assert len(sequence) == sum(exon[1] + 1 - exon[0] for exon in self.exons), (
        #     len(sequence), sum(exon[1] + 1 - exon[0] for exon in self.exons)
        # )

        sequence = self.grouper(sequence, 60)
        fasta.extend(sequence)
        fasta = "\n".join(fasta)
        return fasta

    @staticmethod
    def grouper(iterable, num, fillvalue=""):
        """Collect data into fixed-length chunks or blocks. From core documentation.

        :param iterable: the iterable to be considered for grouping.
        :param num: length of the chunks.
        :type num: int

        :param fillvalue: the default filler for missing positions while grouping.
        :type fillvalue: str
        """
        # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
        args = [iter(iterable)] * num
        return list("".join(x) for x in zip_longest(*args, fillvalue=fillvalue))
