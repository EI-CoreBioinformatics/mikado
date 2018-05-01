# coding: utf-8

"""
Module to parse BED12 objects. Much more basic than what PyBedtools could offer,
but at the same time more pythonic.
"""


import random
import os
from Bio import SeqIO
from Bio import Seq
import Bio.SeqRecord
from . import Parser
from sys import intern
import copy
from ..parsers.GFF import GffLine
from typing import Union
import re
# import numpy as np


# These classes do contain lots of things, it is correct like it is
# pylint: disable=too-many-instance-attributes
class BED12:

    """
    BED12 parsing class.
    """

    def __init__(self, *args: Union[str, list, tuple, GffLine],
                 fasta_index=None,
                 phase=None,
                 sequence=None,
                 transcriptomic=False,
                 max_regression=0,
                 start_adjustment=True,
                 coding=True):

        """
        :param args: the BED12 line.
        :type args: (str, list, tuple, GffLine)

        :param fasta_index: Optional FAI index
        :param sequence: Optional sequence string

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

        self._line = None
        self.header = False
        self.__phase = None  # Initialize to None for the non-transcriptomic objects
        self.__has_start = False
        self.__has_stop = False
        self.__transcriptomic = False
        self.transcriptomic = transcriptomic
        if phase:
            self.phase = phase  # This will check the validity of the phase itself
        self.__strand = None
        self.__max_regression = 0
        # This value will be set when checking for the internal sequence
        # If >=1, i.e. at least one internal stop codon, the ORF is invalid
        self.__internal_stop_codons = 0
        self.chrom = None
        self.start = self.end = self.thick_start = self.thick_end = 0
        self.name = ""
        self.score = 0
        self.strand = None
        self.rgb = ''
        self.block_sizes = [0]
        self.block_starts = [0]
        self.block_count = 1
        self.invalid_reason = ''
        self.fasta_length = None
        self.__in_index = True
        self.max_regression = max_regression
        self.start_adjustment = start_adjustment
        self.coding = coding

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
            if len(self._fields) == 12:
                self.__set_values_from_fields()
                self.header = False
            else:
                self.header = True
                return
        elif isinstance(self._line, GffLine):
            if self._line.header is True:
                self.header = True
                return
            elif self.transcriptomic is False:
                raise TypeError(
                    "GFF lines can be used as BED12-equivalents only in a transcriptomic context.")
            else:
                self.header = False
                if sequence:
                    fasta_length = len(sequence)
                elif fasta_index:
                    fasta_length = len(fasta_index[self._line.chrom])
                else:
                    raise ValueError("Either a sequence or a FAI index are needed")
                self.__set_values_from_gff(fasta_length)

        elif not (isinstance(self._line, list) or isinstance(self._line, tuple)):
            raise TypeError("I need an ordered array, not {0}".format(type(self._line)))
        else:
            self._fields = self._line
            self.__set_values_from_fields()

        if "phase=" in self.name and "coding=" in self.name:  # Hack to include the properties
            groups = dict(re.findall("([^(;|=)]*)=([^;]*)", self.name))
            if groups["phase"].isdigit():
                self.phase = int(groups["phase"])
            else:
                self.phase = None
            if groups["coding"] in ("True", "False"):
                self.coding = eval(groups["coding"])
            else:
                raise ValueError(groups["coding"])
            self.name = groups["ID"]
        elif "coding=" in self.name:
            groups = dict(re.findall("([^(;|=)]*)=([^;]*)", self.name))
            # self.phase = int(groups["phase"])
            if groups["coding"] in ("True", "False"):
                self.coding = eval(groups["coding"])
            else:
                raise ValueError(groups["coding"])
            self.name = groups["ID"]

        self.__check_validity(transcriptomic, fasta_index, sequence)
        if self.invalid and self.coding:
            self.coding = False

    def __set_values_from_fields(self):

        """
        Private method that sets the correct values from the fields derived from the input line.
        :return:
        """
        self.chrom, self.start, self.end, \
            self.name, self.score, self.strand, \
            self.thick_start, self.thick_end, self.rgb, \
            self.block_count, block_sizes, block_starts = self._fields

        # Reduce memory usage
        intern(self.chrom)
        self.start = int(self.start) + 1
        self.end = int(self.end)
        self.score = float(self.score)
        self.thick_start = int(self.thick_start) + 1
        self.thick_end = int(self.thick_end)
        self.block_count = int(self.block_count)
        if isinstance(block_sizes, (str, bytes)):
            self.block_sizes = [int(x) for x in block_sizes.split(",") if x]
        else:
            self.block_sizes = [int(x) for x in block_sizes]
        if isinstance(block_starts, (str, bytes)):
            self.block_starts = [int(x) for x in block_starts.split(",") if x]
        else:
            self.block_starts = [int(x) for x in block_starts]
        self.has_start_codon = None
        self.has_stop_codon = None
        self.start_codon = None
        self.stop_codon = None
        self.fasta_length = len(self)
        return

    def __set_values_from_gff(self, fasta_length):
        """
        Private method that sets the correct values from the fields derived from an input GFF line.
        :return:
        """

        (self.chrom, self.thick_start,
         self.thick_end, self.strand, self.name) = (self._line.chrom,
                                                    self._line.start,
                                                    self._line.end, self._line.strand, self._line.id)
        intern(self.chrom)
        assert self.name is not None
        self.start = 1
        self.end = fasta_length
        self.score = self._line.score
        self.rgb = None
        self.block_count = 1
        self.block_sizes = [self.thick_end - self.thick_start +1]
        self.block_starts = [self.thick_start]
        self.has_start_codon = None
        self.has_stop_codon = None
        self.start_codon = None
        self.stop_codon = None
        self.fasta_length = fasta_length
        return

    def __check_validity(self, transcriptomic, fasta_index, sequence):
        """
        Private method that checks that the BED12 object has been instantiated correctly.

        :return:
        """

        if transcriptomic is True:
            self.has_start_codon = False
            self.has_stop_codon = False

            # if self.strand == "-":
            #     self.thick_end -= 3

        if transcriptomic is True and self.coding is True and (fasta_index is not None or sequence is not None):

            if sequence is not None:
                self.fasta_length = len(sequence)
                if isinstance(sequence, str):
                    sequence = Seq.Seq(sequence)
            else:

                if self.id not in fasta_index:
                    self.__in_index = False
                    return

                self.fasta_length = len(fasta_index[self.id])
                sequence = fasta_index[self.id].seq

            # Just double check that the sequence length is the same as what the BED would suggest
            if self.invalid is True:
                return

            if self.strand == "+":
                orf_sequence = sequence[
                               (self.thick_start - 1 if not self.phase else self.start + self.phase - 1):self.thick_end]
            elif self.strand == "-":
                orf_sequence = sequence[(self.thick_start - 1):(
                    self.thick_end if not self.phase else self.end - self.phase)].reverse_complement()
            else:
                pass

            self.start_codon = str(orf_sequence)[:3].upper()
            self.stop_codon = str(orf_sequence[-3:]).upper()

            if self.start_codon == "ATG":
                self.has_start_codon = True
                self.phase = 0
            else:
                # We are assuming that if a methionine can be found it has to be
                # internally, not externally, to the ORF

                self.has_start_codon = False

                if self.start_adjustment is True:
                    self._adjust_start(orf_sequence)

            if self.stop_codon in ("TAA", "TGA", "TAG"):
                self.has_stop_codon = True
            else:
                self.has_stop_codon = False
                # Expand the ORF to include the end of the sequence
                if self.end - self.thick_end <= 2:
                    self.thick_end = self.end

            translated_seq = orf_sequence[:-3].translate()
            self.__internal_stop_codons = str(translated_seq).count("*")

            if self.invalid is True:
                return

    def _adjust_start(self, orf_sequence):

        for pos in range(3,
                         int(len(orf_sequence) * self.max_regression),
                         3):
            if orf_sequence[pos:pos + 3] == "ATG":
                # Now we have to shift the start accordingly
                self.has_start_codon = True
                if self.strand == "+":
                    self.thick_start += pos
                else:
                    # TODO: check that this is right and we do not have to do some other thing
                    self.thick_end -= pos
                break
            else:
                continue

        if self.has_start_codon is False:
            # The validity will be automatically checked
            if self.strand == "+":
                self.phase = max(self.thick_start - self.start, 0)
                self.thick_start = self.start
            else:
                if self.end - self.thick_end <= 2:
                    self.phase = self.end - self.thick_end
                    self.thick_end = self.end
                else:
                    self.phase = 0

    def __str__(self):

        if self.header is True:
            if self._line is not None:
                return self._line
            else:
                return "#"

        line = [self.chrom, self.start - 1, self.end, self.name]
        if not self.score:
            line.append(0)
        else:
            line.append(self.score)
        if self.strand is None:
            line.append(".")
        else:
            line.append(self.strand)
        line.extend([self.thick_start - 1, self.thick_end])
        if not self.rgb:
            line.append(0)
        else:
            line.append(self.rgb)
        line.append(self.block_count)
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

    def copy(self):

        return copy.deepcopy(self)

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

        if self.transcriptomic is True:
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

        if self.__internal_stop_codons >= 1:
            self.invalid_reason = "{} internal stop codons found".format(self.__internal_stop_codons)
            return True

        if self.transcriptomic is True and self.__in_index is False:
            self.invalid_reason = "{} not found in the index!".format(self.chrom)
            return True

        assert isinstance(self.thick_start, int)
        assert isinstance(self.thick_end, int)
        assert isinstance(self.start, int)
        assert isinstance(self.end, int)

        if self.thick_start < self.start or self.thick_end > self.end:
            if self.thick_start == self.thick_end == self.block_sizes[0] == 0:
                pass
            else:
                invalid = "thickStart {0} <start {1}: {2}; end {3} <thickEnd {4} {5}"
                self.invalid_reason = invalid.format(self.thick_start,
                                                     self.start,
                                                     self.thick_start < self.start,
                                                     self.end,
                                                     self.thick_end,
                                                     self.thick_end > self.end)
                return True

        if self.fasta_length is None:
            self.fasta_length = len(self)

        if len(self) != self.fasta_length:
            self.invalid_reason = "Fasta length != BED length: {0} vs. {1}".format(
                self.fasta_length,
                len(self)
            )
            return True

        if self.transcriptomic is True and (self.cds_len - self.phase) % 3 != 0 and self.thick_end != self.end:
            self.invalid_reason = "Invalid CDS length: {0} % 3 = {1}".format(
                self.cds_len - self.phase,
                (self.cds_len - self.phase) % 3
            )
            return True

        self.invalid_reason = ''
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
        if value and self.phase is None:
            self.phase = 0
        elif not value:
            self.phase = None

    @property
    def phase(self):
        """This property is used for transcriptomic BED objects
        and indicates what the phase of the transcript is.
        So a BED object with an open 5'ORF whose first codon
        starts at the 1st base would have frame 1, at the second
        frame 2. In all other cases, the frame has to be 0.
        If the BED object is not transcriptomic, its frame is null
        (None).
        """

        return self.__phase

    @phase.setter
    def phase(self, val):

        if val not in (None, 0, 1, 2):
            raise ValueError("Invalid frame specified for {}: {}. Must be None or 0, 1, 2".format(
                self.name, val))
        elif self.transcriptomic is True and val not in (0, 1, 2):
            raise ValueError("A transcriptomic BED cannot have null frame.")
        self.__phase = val

    @property
    def _max_regression(self):
        """
        This property is used to indicate how far downstream we should go in the
          FASTA sequence to find a valid start codon, in terms of percentage
          of the the cDNA sequence. So eg in a 300 nt cDNA a max_regression value
          of 0.3 would instruct the class to look only for the first 90 bps for
          a Met.
        """

        return self.__max_regression

    @_max_regression.setter
    def _max_regression(self, value):
        if not (isinstance(value, (int, float)) and 0 <= value <= 1):
            raise ValueError(
                "Invalid value specified for _max_regression (must be between 0 and 1): {}".format(value))
        self.__max_regression = value

    def expand(self, sequence, upstream, downstream):

        """This method will expand a """
        # assert len(sequence) >= len(self)
        assert len(sequence) == len(self) + upstream + downstream, (len(sequence),
                                                                    len(self),
                                                                    upstream,
                                                                    downstream,
                                                                    len(self) + upstream + downstream)
        if len(self) == len(sequence):
            return
        if self.transcriptomic is False:
            raise ValueError("I cannot expand a non-transcriptomic BED12!")
        if self.strand == "-":
            raise NotImplementedError("I can only expand ORFs on the sense strand for the moment")

        self.fasta_length = len(sequence)

        # I presume that the sequence is already in the right orientation

        old_sequence = sequence[upstream:len(self) + upstream]

        self.start_codon = str(old_sequence[self.thick_start + self.phase:self.thick_start + self.phase + 3]).upper()
        last_codon_start = self.thick_end + ((self.thick_end - self.thick_start + 1 + self.phase) % 3 - 3)

        self.stop_codon = str(old_sequence[last_codon_start:self.thick_end]).upper()

        assert len(self.stop_codon) <= 3 and len(self.stop_codon) > 0, self.stop_codon

        # Now expand
        self.end = len(sequence)
        self.thick_start += upstream
        self.thick_end += upstream
        last_codon_start += upstream
        if self.start_codon != "ATG":
            for pos in range(self.thick_start - self.phase,
                             0,
                             -3):
                codon = sequence[pos:pos + 3]
                self.thick_start = pos
                if codon == "ATG":
                    # self.thick_start = pos
                    self.start_codon = codon
                    self.__has_start = True
                    break

        if self.start_codon != "ATG":
            self.phase = self.thick_start % 3
            self.thick_start = 1
        else:
            self.phase = 0
            self.__has_start = True

        if self.stop_codon not in ("TAA", "TGA", "TAG"):
            for pos in range(last_codon_start,
                             self.end,
                             3):
                codon = sequence[pos:pos + 3]
                if codon in ("TAA", "TGA", "TAG"):
                    self.thick_end = pos + 3
                    self.stop_codon = codon
                    self.__has_stop = True
                    break
        if self.stop_codon not in ("TAA", "TGA", "TAG"):
            self.thick_end = self.end

        self.block_sizes = [self.thick_end - self.thick_start]
        self.block_starts = [self.thick_start]

        return

    @property
    def blocks(self):

        """This will return the coordinates of the blocks, with a 1-offset (as in GFF3)"""

        # First thing: calculate where each start point will be
        _blocks = []
        starts = [_ + self.start - 1 for _ in self.block_starts]
        for pos in range(self.block_count):
            _blocks.append((starts[pos] + 1, starts[pos] + self.block_sizes[pos]))

        return _blocks

    def to_transcriptomic(self, sequence=None, fasta_index=None, start_adjustment=False):

        """This method will return a transcriptomic version of the BED12. If the object is already transcriptomic,
        it will return itself."""

        if self.transcriptomic is True:
            return self

        # First six fields of a BED object

        # Now we have to calculate the thickStart
        block_count = self.block_count

        # Now we have to calculate the thickStart, thickEnd ..
        tStart, tEnd = None, None
        seen = 0
        # if self.strand == "+":
        #     adder = (0, 1)
        # else:
        #     adder = (0, 1)

        for block in self.blocks:
            if tStart and tEnd:
                # We have calculated them
                break
            if not tStart and block[0] <= self.thick_start <= block[1]:
                tStart = seen + self.thick_start - block[0] + 0  # adder[0]
            if not tEnd and block[0] <= self.thick_end <= block[1]:
                tEnd = seen + self.thick_end - block[0] + 1
            seen += block[1] - block[0] + 1

        assert tStart is not None and tEnd is not None, (tStart, tEnd, self.thick_start, self.thick_end, self.blocks)

        if self.strand == "+":
            bsizes = self.block_sizes[:]
        else:
            bsizes = list(reversed(self.block_sizes[:]))
            tStart, tEnd = sum(self.block_sizes) - tEnd, sum(self.block_sizes) - tStart

        bstarts = [0]
        for bs in bsizes[:-1]:
            bstarts.append(bs + bstarts[-1])
        assert len(bstarts) == len(bsizes) == self.block_count, (bstarts, bsizes, self.block_count)

        if self.coding:
            new_name = "ID={};coding={};phase={}".format(self.name.split(";")[0],
                                                         self.coding,
                                                         self.phase if self.phase is not None else 0)
        else:
            new_name = "ID={};coding={}".format(self.name.split(";")[0], self.coding)

        new = list((self.name.split(";")[0],
                    0,
                    sum(self.block_sizes),
                    new_name,
                    self.score,
                    "+"))

        new.extend(list((
            tStart,
            tEnd,
            self.rgb,
            self.block_count,
            bsizes,
            bstarts
        )))
        new = BED12(new,
                    phase=self.phase,
                    sequence=sequence,
                    coding=self.coding,
                    fasta_index=fasta_index,
                    transcriptomic=True,
                    start_adjustment=start_adjustment)
        # assert new.invalid is False
        assert isinstance(new, type(self)), type(new)
        return new


class Bed12Parser(Parser):
    """Parser class for a Bed12Parser file.
    It accepts optionally a fasta index which is used to
    determine whether an ORF has start/stop codons."""

    def __init__(self, handle,
                 fasta_index=None,
                 transcriptomic=False,
                 max_regression=0,
                 is_gff=False):
        """
        Constructor method.
        :param handle: the input BED file.

        :param fasta_index: optional FAI file

        :param transcriptomic: flag. If set to True, it indicates
        that the BED file contains ORF information (like that derived
        from transdecoder) and is in transcriptomic rather than genomic coordinates.
        :type transcriptomic: bool

        :param max_regression: parameter to pass directly to BED12, indicating how much
        we should backtrack in the sequence to find a valid start codon
        """

        Parser.__init__(self, handle)
        self.transcriptomic = transcriptomic
        self.__max_regression = 0
        self._max_regression = max_regression

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
        self.__closed = False
        self.header = False
        self._is_bed12 = (not is_gff)

    def __iter__(self):
        return self

    def __next__(self):

        if self._is_bed12 is True:
            return self.bed_next()
        else:
            return self.gff_next()

    def bed_next(self):
        """

        :return:
        """

        bed12 = None
        while bed12 is None:
            line = self._handle.readline()
            if line == '':
                raise StopIteration
            bed12 = BED12(line,
                          fasta_index=self.fasta_index,
                          transcriptomic=self.transcriptomic,
                          max_regression=self._max_regression)
        return bed12

    def gff_next(self):
        """

        :return:
        """

        bed12 = None
        while bed12 is None:
            line = self._handle.readline()
            if line == "":
                raise StopIteration
            line = GffLine(line)

            if line.feature != "CDS":
                continue
            # Compatibility with BED12
            bed12 = BED12(line,
                          fasta_index=self.fasta_index,
                          transcriptomic=self.transcriptomic,
                          max_regression=self._max_regression)
        # raise NotImplementedError("Still working on this!")
        return bed12

    def close(self):
        if self.__closed is False:
            self._handle.close()
            self.__closed = True

    @property
    def _max_regression(self):
        """
        This property is used to indicate how far downstream we should go in the
          FASTA sequence to find a valid start codon, in terms of percentage
          of the the cDNA sequence. So eg in a 300 nt cDNA a max_regression value
          of 0.3 would instruct the class to look only for the first 90 bps for
          a Met.
        """

        return self.__max_regression

    @_max_regression.setter
    def _max_regression(self, value):
        if not (isinstance(value, (int, float)) and 0 <= value <= 1):
            raise ValueError(
                "Invalid value specified for _max_regression (must be between 0 and 1): {}".format(value))
        self.__max_regression = value