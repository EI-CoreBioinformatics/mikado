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

"""Generic module for parsing Bed12Parser files."""


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

        :param blockStarts: list of start positions for the blocks. Its length must be equal to blockCount
        :type blockStarts: list(int)
        
        
        Additional parameters calculated inside the class:
        
        :param fasta_length: length of the feature (see len() method)
        :type fasta_length: int
        
        :param has_start_codon: flag. For transcriptomic BED12, it indicates whether we have a start codon or not.
        :type has_start_codon: bool

        :param has_stop_codon: flag. For transcriptomic BED12, it indicates whether we have a stop codon or not.
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
        self.__strand = None

        if len(args) == 0:
            self.header = True
            return

        self._line = args[0]
        if type(self._line) is str or self._line is None:
            if self._line is None or len(self._line) == 0 or self._line[0] == "#":
                self.header = True
                return

            self._fields = self._line.rstrip().split("\t")
        elif type(self._line) not in (list, tuple):
            raise TypeError("I need an ordered array, not {0}".format(type(self._line)))
        if len(self._fields) != 12:
            self.header = True
            return
            # raise ValueError("Erroneous number of fields detected")

        self.transcriptomic = transcriptomic
        self.header = False
        self.chrom, self.start, self.end, \
            self.name, self.score, self.strand, \
            self.thickStart, self.thickEnd, self.rgb, \
            self.blockCount, self.blockSizes, self.blockStarts = self._fields

        self.start = int(self.start) + 1
        self.end = int(self.end)
        self.score = float(self.score)
        self.thickStart = int(self.thickStart) + 1
        self.thickEnd = int(self.thickEnd)
        self.blockCount = int(self.blockCount)
        self.blockSizes = [int(x) for x in self.blockSizes.split(",")]
        self.blockStarts = [int(x) for x in self.blockStarts.split(",")]
        self.has_start_codon = None
        self.has_stop_codon = None
        self.start_codon = None
        self.stop_codon = None
        self.fasta_length = len(self)

        if transcriptomic is True:
            self.has_start_codon = False
            self.has_stop_codon = False

            if self.strand == "-":
                self.thickEnd -= 3

        if self.invalid is True:
            return

        if transcriptomic is True and fasta_index is not None:
            assert self.id in fasta_index
            self.fasta_length = len(fasta_index[self.id])
            if self.invalid is True:
                return

            start_codon = fasta_index[self.id][self.thickStart - 1:self.thickStart + 2]
            stop_codon = fasta_index[self.id][self.thickEnd:self.thickEnd + 3]
            if self.strand == "-":
                start_codon = fasta_index[self.id][self.thickStart - 1:self.thickStart + 2]
                start_codon = start_codon.reverse_complement()
                stop_codon = stop_codon.reverse_complement()
                start_codon, stop_codon = stop_codon, start_codon
            self.start_codon = str(start_codon.seq)
            self.stop_codon = str(stop_codon.seq)

            if self.start_codon == "ATG":
                self.has_start_codon = True
            else:
                self.has_start_codon = False
            if self.stop_codon in ("TAA", "TGA", "TAG"):
                self.has_stop_codon = True
            else:
                self.has_stop_codon = False

            #             #This is a bug in TD by which sometimes truncates ORFs even when there would be a read through
            if self.has_stop_codon is False and (self.strand != "-" and self.thickEnd < len(self) - 2) or (
                    self.strand == "-" and self.thickStart > 2):
                sequence = fasta_index[self.id]
                tstart, tend = self.thickStart, self.thickEnd
                if self.strand != "-":
                    num = self.thickEnd + 3
                    for num in range(self.thickEnd + 3, self.end, 3):
                        codon = sequence[num:num + 3]
                        if str(codon.seq) in ("TAA", "TGA", "TAG"):
                            self.has_stop_codon = True
                            break
                    self.thickEnd = num - 3
                else:
                    num = self.thickStart
                    # This while loop is necessary b/c range does not function backwards
                    # and reversed mingles poorly with non-unary steps
                    # i.e reversed(range(1,10,3)) = [9,6,3,0], not [10,7,4,1] ...
                    while num > self.start:
                        num -= 3
                        codon = sequence[num - 3:num]
                        # Reversed version, save on reversal, should save time
                        if str(codon.seq) in ("TTA", "TCA", "CTA"):
                            self.has_stop_codon = True
                            break
                    self.thickStart = num
                assert self.invalid is False, ((tstart, tend), (self.strand, self.thickStart, self.thickEnd))

        assert self.blockCount == len(self.blockStarts) == len(self.blockSizes)

    def __str__(self):

        line = [self.chrom, self.start - 1, self.end, self.name, self.score]
        if self.strand is None:
            line.append(".")
        else:
            line.append(self.strand)
        line.extend([self.thickStart - 1, self.thickEnd, self.rgb, self.blockCount])
        line.append(",".join([str(x) for x in self.blockSizes]))
        line.append(",".join([str(x) for x in self.blockStarts]))
        return "\t".join([str(x) for x in line])

    def __eq__(self, other):
        for key in ["chrom", "strand", "start", "end", "thickStart", "thickEnd", "blockCount", "blockSizes",
                    "blockStarts"]:
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
        :rtype None
        :rtype str
        """
        return self.__strand

    @strand.setter
    def strand(self, strand: str):
        """
        Setter for strand. It verifies that the value is correct.
        :param strand: New strand value
        :type strand: str
        :type strand: None
        """

        if strand in (".", "?"):
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
        return self.thickEnd - self.thickStart + 1

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
        Property. True if the BED12 is transcriptomic and has both start and stop codon, False otherwise.
        :rtype bool
        """
        return self.has_stop_codon and self.has_start_codon

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

    @property
    def invalid(self):
        """
        Property. It performs basic checks on the BED line to verify its integrity.
        :rtype bool
        """
        if self.thickStart < self.start or self.thickEnd > self.end:
            return True

        if "fasta_length" not in self.__dict__:
            self.fasta_length = len(self)

        if len(self) != self.fasta_length:
            return True
        if self.transcriptomic is True and (self.thickEnd - self.thickStart + 1) % 3 != 0:
            return True
        return False

    @property
    def transcriptomic(self):
        """
        Flag. If set to True, it indicates the BED contains transcriptomic rather than genomic coordinates.
        :rtype bool
        """
        return self.__transcriptomic

    @transcriptomic.setter
    def transcriptomic(self, value):
        """
        Setter for transcriptomic. A valid value must be boolean.
        :type value: bool
        """

        if type(value) is not bool:
            raise ValueError("Invalid value: {0}".format(value))
        self.__transcriptomic = value


class Bed12Parser(Parser):
    """Parser class for a Bed12Parser file.
    It accepts optionally a fasta index which is used to determine whether an ORF has start/stop codons."""

    def __init__(self, handle, fasta_index=None, transcriptomic=False):
        """
        Constructor method.
        :param handle: the input BED file.

        :param fasta_index: optional FAI file

        :param transcriptomic: flag. If set to True, it indicates that the BED file contains ORF information
        (like that derived from transdecoder) and is in transcriptomic rather than genomic coordinates.
        :type transcriptomic: bool
        """

        super().__init__(handle)
        self.transcriptomic = transcriptomic

        if type(fasta_index) is dict:
            # check that this is a bona fide dictionary ...
            assert type(fasta_index[random.sample(fasta_index.keys(), 1)]) is Bio.SeqRecord.SeqRecord
        elif fasta_index is not None:
            if type(fasta_index) is str:
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
