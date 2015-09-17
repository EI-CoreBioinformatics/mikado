# coding: utf-8

"""
Module to parse BED12 objects. Much more basic than what PyBedtools could offer,
but at the same time more pythonic.
"""


import random
import os
from Bio import SeqIO
import Bio.SeqRecord
from mikado_lib.parsers import Parser


# These classes do contain lots of things, it is correct like it is
# pylint: disable=too-many-instance-attributes
class BED12:

    """
    BED12 parsing class.
    """

    def __init__(self, *args: str, fasta_index=None, transcriptomic=False):

        """
        :param args: the BED12 line.
        :type args: str

        :param fasta_index: Optional FAI index

        :param transcriptomic: boolean flag
        :type transcriptomic: bool

        Constructor method.

        Each instance will have:

        :param chrom: chromosome
        :type chrom: str

        :param header: whether the line is a header line or not
        :type header: bool

        :param start: start position
        :type start: int

        :param end: end position
        :type end: int

        :param name: Name of the feature
        :type name: str

        :param score: score assigned to the feature.
        :type score: float
        :type score: None

        :param strand of the feature

        :param thickStart: "internal" start of the feature (e.g. CDS start)
        :type thickStart: int

        :param thickEnd: "internal" end of the feature (e.g. CDS end)
        :type thickEnd: int

        :param rgb: RGB color scheme. Currently not checked
        :type rgb: str

        :param blockCount: number of blocks (e.g. exons)
        :type blockCount: int

        :param blockSizes: sizes of the blocks (e.g. exons). Its length must be equal to blockCount
        :type blockSizes: list(int)

        :param blockStarts: list of start positions for the blocks.
        Its length must be equal to blockCount
        :type blockStarts: list(int)

        Additional parameters calculated inside the class:

        :param fasta_length: length of the feature (see len() method)
        :type fasta_length: int

        :param has_start_codon: flag. For transcriptomic BED12, it indicates
        whether we have a start codon or not.
        :type has_start_codon: bool

        :param has_stop_codon: flag. For transcriptomic BED12, it indicates
        whether we have a stop codon or not.
        :type has_stop_codon: bool

        :param start_codon: string of the start codon, if found
        :type start_codon: None
        :type start_codon: str

        :param stop_codon: string of the stop codon, if found
        :type stop_codon: None
        :type stop_codon: str
        """

        self.__has_start = False
        self.__has_stop = False
        self.__transcriptomic = False
        self.__invalid = None
        self.__strand = None
        self.__internal_stop_codons = 0
        self.chrom = None
        self.start = self.end = self.thick_start = self.thick_end = 0
        self.score = 0
        self.strand = None
        self.rgb = ''
        self.block_sizes = [0]
        self.block_starts = [0]
        self.block_count = 1
        self.invalid_reason = ''
        self.fasta_length = None

        if len(args) == 0:
            self.header = True
            return

        self._line = args[0]
        if isinstance(self._line, str) or self._line is None:
            if self._line is None:
                self._line = ''
            self._line = self._line.rstrip()
            if len(self._line) == 0 or self._line[0] == "#":
                self.header = True
                return
            self._fields = self._line.split("\t")
        elif not (isinstance(self._line, list) or isinstance(self._line, tuple)):
            raise TypeError("I need an ordered array, not {0}".format(type(self._line)))

        if len(self._fields) != 12:
            self.header = True
            return

        self.transcriptomic = transcriptomic
        self.header = False

        self.__set_values_from_fields()
        self.__check_validity(transcriptomic, fasta_index)

    def __set_values_from_fields(self):

        """
        Private method that sets the correct values from the fields derived from the input line.
        :return:
        """
        self.chrom, self.start, self.end, \
            self.name, self.score, self.strand, \
            self.thick_start, self.thick_end, self.rgb, \
            self.block_count, block_sizes, block_starts = self._fields

        self.chrom = self._fields[0]

        self.start = int(self.start) + 1
        self.end = int(self.end)
        self.score = float(self.score)
        self.thick_start = int(self.thick_start) + 1
        self.thick_end = int(self.thick_end)
        self.block_count = int(self.block_count)
        self.block_sizes = [int(x) for x in block_sizes.split(",")]
        self.block_starts = [int(x) for x in block_starts.split(",")]
        self.has_start_codon = None
        self.has_stop_codon = None
        self.start_codon = None
        self.stop_codon = None
        self.fasta_length = len(self)
        return

    def __check_validity(self, transcriptomic, fasta_index):
        """
        Private method that checks that the BED12 object has been instantiated correctly.

        :return:
        """

        if transcriptomic is True:
            self.has_start_codon = False
            self.has_stop_codon = False

            # if self.strand == "-":
            #     self.thick_end -= 3

        if self.invalid is True:
            return

        if transcriptomic is True and fasta_index is not None:
            assert self.id in fasta_index
            self.fasta_length = len(fasta_index[self.id])
            if self.invalid is True:
                return
            sequence = fasta_index[self.id].seq

            orf_sequence = sequence[self.thick_start-1:self.thick_end]
            if self.strand == "-":
                orf_sequence = orf_sequence.reverse_complement()

            self.start_codon = str(orf_sequence)[:3]
            self.stop_codon = str(orf_sequence[-3:])

            if self.start_codon == "ATG":
                self.has_start_codon = True
            else:
                self.has_start_codon = False
            if self.stop_codon in ("TAA", "TGA", "TAG"):
                self.has_stop_codon = True
            else:
                self.has_stop_codon = False

            translated_seq = orf_sequence.translate()
            self.__internal_stop_codons = str(translated_seq).count("*")

            # if self.has_stop_codon is False and \
            #         (self.strand != "-" and self.thick_end < len(self) - 2) or\
            #         (self.strand == "-" and self.thick_start > 2):
            #     self.__recheck_stop_codon(sequence)

    # def __recheck_stop_codon(self, sequence):
    #
    #     """
    #     Due to a bug in TransDecoder, sometimes valid ORFs with valid stop
    #     codons are tuncated. This private method rechecks the consistency
    #      of the ORF against the transcript underlying sequence.
    #     :param fasta_index:
    #     :return:
    #     """
    #
    #     if self.strand != "-":
    #         num = self.thick_end + 3
    #         for num in range(self.thick_end + 3, self.end, 3):
    #             codon = sequence[num:num + 3]
    #             if str(codon) in ("TAA", "TGA", "TAG"):
    #                 self.has_stop_codon = True
    #                 break
    #         self.thick_end = num - 3
    #     else:
    #         num = self.thick_start
    #         # This while loop is necessary b/c range does not function backwards
    #         # and reversed mingles poorly with non-unary steps
    #         # i.e reversed(range(1,10,3)) = [9,6,3,0], not [10,7,4,1] ...
    #         while num > self.start:
    #             num -= 3
    #             codon = sequence[num - 3:num]
    #             # Reversed version, save on reversal, should save time
    #             if str(codon) in ("TTA", "TCA", "CTA"):
    #                 self.has_stop_codon = True
    #                 break
    #         self.thick_start = num
    #     orf_sequence = sequence[self.thick_start-1:self.thick_end]
    #     if self.strand == "-":
    #         orf_sequence = orf_sequence.reverse_complement()
    #
    #     translated_seq = orf_sequence.translate()
    #     self.__internal_stop_codons = str(translated_seq).count("*")
    #
    #     if self.invalid is True:
    #         self.invalid_reason = "Wrong CDS detection"

    def __str__(self):

        line = [self.chrom, self.start - 1, self.end, self.name, self.score]
        if self.strand is None:
            line.append(".")
        else:
            line.append(self.strand)
        line.extend([self.thick_start - 1, self.thick_end, self.rgb, self.block_count])
        line.append(",".join([str(x) for x in self.block_sizes]))
        line.append(",".join([str(x) for x in self.block_starts]))
        return "\t".join([str(x) for x in line])

    def __eq__(self, other):
        for key in ["chrom", "strand", "start",
                    "end", "thick_start", "thick_end",
                    "block_count", "block_sizes",
                    "block_starts"]:
            if getattr(self, key) != getattr(other, key):
                return False
        return True

    def __hash__(self):
        return super().__hash__()

    def __len__(self):
        return self.end - self.start + 1

    @property
    def strand(self):
        """
        Strand of the feature. It must be one of None,+,-
        :rtype None | str
        """
        return self.__strand

    @strand.setter
    def strand(self, strand: str):
        """
        Setter for strand. It verifies that the value is correct.
        :param strand: New strand value
        :type strand: str | None
        """

        if strand in (".", "?", None):
            self.__strand = None
        elif strand in ("+", "-"):
            self.__strand = strand
        else:
            raise ValueError("Erroneous strand provided: {0}".format(self.strand))

    @property
    def cds_len(self):
        """
        Return the length of the internal feature i.e. the CDS:
        thickEnd-thickStart+1

        :rtype int
        """
        return self.thick_end - self.thick_start + 1

    @property
    def has_start_codon(self):
        """
        Property. True if the interval contains a start codon.
        :rtype bool
        :rtype None
        """

        return self.__has_start

    @has_start_codon.setter
    def has_start_codon(self, value: bool):
        """
        Setter for has_stop_codon. Valid values are boolean or None
        :param value: boolean flag
        :type value: bool
        :type value: None
        """

        if value not in (None, True, False):
            raise ValueError()
        self.__has_start = value

    @property
    def has_stop_codon(self):
        """
        Property. True if the interval contains a termination codon.
        :rtype bool
        :rtype None
        """
        return self.__has_stop

    @has_stop_codon.setter
    def has_stop_codon(self, value: bool):
        """
        Setter for has_stop_codon. Valid values are boolean.
        :param value: boolean flag
        :type value: bool
        :type value: None
        """
        if value not in (None, True, False):
            raise ValueError()
        self.__has_stop = value

    @property
    def full_orf(self):
        """
        Property. True if the BED12 is transcriptomic and has
        both start and stop codon, False otherwise.
        :rtype bool
        """
        return self.has_stop_codon and self.has_start_codon

    # pylint: disable=invalid-name
    @property
    def id(self):
        """
        Property. It returns the name of the feature.
        :rtype str
        """

        if self.transcriptomic:
            return self.chrom
        else:
            return self.name
    # pylint: enable=invalid-name

    @property
    def invalid(self):
        """
        Property. It performs basic checks on the BED line to verify its integrity.
        :rtype bool
        """

        if self.__internal_stop_codons > 1:
            return True

        if self.thick_start < self.start or self.thick_end > self.end:
            invalid = "thickStart {0} <start {1}: {2}; end {3} <thickEnd {4} {5}"
            self.invalid_reason = invalid.format(self.thick_start,
                                                 self.start,
                                                 self.thick_start < self.start,
                                                 self.end,
                                                 self.thick_end,
                                                 self.thick_end > self.end)
            self.__invalid = True
            return True

        if self.fasta_length is None:
            self.fasta_length = len(self)

        if len(self) != self.fasta_length:
            self.invalid_reason = "Fasta length != BED length: {0} vs. {1}".format(
                self.fasta_length,
                len(self)
            )
            self.__invalid = True
            return True

        if self.transcriptomic is True and self.cds_len % 3 != 0:
            self.invalid_reason = "Invalid CDS length: {0} % 3 = {1}".format(
                self.cds_len,
                self.cds_len % 3
            )
            self.__invalid = True
            return True

        self.__invalid = False
        return False

    @property
    def transcriptomic(self):
        """
        Flag. If set to True, it indicates the BED contains
        transcriptomic rather than genomic coordinates.
        :rtype bool
        """
        return self.__transcriptomic

    @transcriptomic.setter
    def transcriptomic(self, value):
        """
        Setter for transcriptomic. A valid value must be boolean.
        :type value: bool
        """

        if not isinstance(value, bool):
            raise ValueError("Invalid value: {0}".format(value))
        self.__transcriptomic = value


class Bed12Parser(Parser):
    """Parser class for a Bed12Parser file.
    It accepts optionally a fasta index which is used to
    determine whether an ORF has start/stop codons."""

    def __init__(self, handle, fasta_index=None, transcriptomic=False):
        """
        Constructor method.
        :param handle: the input BED file.

        :param fasta_index: optional FAI file

        :param transcriptomic: flag. If set to True, it indicates
        that the BED file contains ORF information (like that derived
        from transdecoder) and is in transcriptomic rather than genomic coordinates.
        :type transcriptomic: bool
        """

        Parser.__init__(self, handle)
        self.transcriptomic = transcriptomic

        if isinstance(fasta_index, dict):
            # check that this is a bona fide dictionary ...
            assert isinstance(
                fasta_index[random.sample(fasta_index.keys(), 1)],
                Bio.SeqRecord.SeqRecord)
        elif fasta_index is not None:
            if isinstance(fasta_index, str):
                assert os.path.exists(fasta_index)
                fasta_index = SeqIO.index(fasta_index, "fasta")
            else:
                assert "SeqIO" in repr(fasta_index) and "index" in repr(fasta_index)

        self.fasta_index = fasta_index

        self.header = False

    def __iter__(self):
        return self

    def __next__(self):
        line = self._handle.readline()
        if line == '':
            raise StopIteration
        return BED12(line, fasta_index=self.fasta_index, transcriptomic=self.transcriptomic)
