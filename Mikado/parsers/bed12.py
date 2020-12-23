# coding: utf-8

"""
Module to parse BED12 objects. Much more basic than what PyBedtools could offer,
but at the same time more pythonic.
"""

from time import sleep
import os
from Bio import Seq
import Bio.SeqRecord
from . import Parser
from sys import intern
import copy
from ..parsers.GFF import GffLine
from typing import Union
import re
import pysam
import functools
from Bio import BiopythonWarning
import warnings
from Bio.Data.IUPACData import ambiguous_dna_letters as _ambiguous_dna_letters
from Bio.Data.IUPACData import ambiguous_rna_letters as _ambiguous_rna_letters
from Bio.Data import CodonTable
import multiprocessing as mp
import msgpack
import logging
import logging.handlers as logging_handlers
from ..utilities.log_utils import create_null_logger
import pyfaidx
import zlib
import numpy as np
import random


backup_valid_letters = set(_ambiguous_dna_letters.upper() + _ambiguous_rna_letters.upper())
standard = CodonTable.ambiguous_dna_by_id[1]
standard.start_codons = ["ATG"]


@functools.lru_cache(typed=True, maxsize=2**10)
def get_tables(table, to_stop=False, gap=None):
    forward_table = table.forward_table.forward_table.copy()
    stop_codons = set(table.stop_codons)
    dual_coding = [c for c in stop_codons if c in forward_table]
    if dual_coding:
        c = dual_coding[0]
        if to_stop:
            raise ValueError("You cannot use 'to_stop=True' with this table "
                             "as it contains {} codon(s) which can be both "
                             " STOP and an  amino acid (e.g. '{}' -> '{}' or "
                             "STOP)."
                             .format(len(dual_coding), c, forward_table[c]))
        warnings.warn("This table contains {} codon(s) which code(s) for both "
                      "STOP and an amino acid (e.g. '{}' -> '{}' or STOP). "
                      "Such codons will be translated as amino acid."
                      .format(len(dual_coding), c, forward_table[c]),
                      BiopythonWarning)

    for stop in stop_codons:
        forward_table[stop] = "*"
    if gap is not None:
        forward_table[gap * 3] = "*"

    if table.nucleotide_alphabet is not None:
        valid_letters = set(table.nucleotide_alphabet.upper())
    else:
        # Assume the worst case, ambiguous DNA or RNA:
        valid_letters = backup_valid_letters

    getter = np.vectorize(forward_table.get, otypes=["<U"])

    return forward_table, getter, valid_letters


def _translate_str(sequence, table, stop_symbol="*", to_stop=False, cds=False, pos_stop="X", gap=None):
    """Translate nucleotide string into a protein string (PRIVATE).

    Arguments:
     - sequence - a string
     - table - a CodonTable object (NOT a table name or id number)
     - stop_symbol - a single character string, what to use for terminators.
     - to_stop - boolean, should translation terminate at the first
       in frame stop codon?  If there is no in-frame stop codon
       then translation continues to the end.
     - pos_stop - a single character string for a possible stop codon
       (e.g. TAN or NNN)
     - cds - Boolean, indicates this is a complete CDS.  If True, this
       checks the sequence starts with a valid alternative start
       codon (which will be translated as methionine, M), that the
       sequence length is a multiple of three, and that there is a
       single in frame stop codon at the end (this will be excluded
       from the protein sequence, regardless of the to_stop option).
       If these tests fail, an exception is raised.
     - gap - Single character string to denote symbol used for gaps.
       Defaults to None.

    Returns a string.

    e.g.

    >>> from Bio.Data import CodonTable
    >>> table = CodonTable.ambiguous_dna_by_id[1]
    >>> _translate_str("AAA", table)
    'K'
    >>> _translate_str("TAR", table)
    '*'
    >>> _translate_str("TAN", table)
    'X'
    >>> _translate_str("TAN", table, pos_stop="@")
    '@'
    >>> _translate_str("TA?", table)
    Traceback (most recent call last):
       ...
    Bio.Data.CodonTable.TranslationError: Codon 'TA?' is invalid

    In a change to older versions of Biopython, partial codons are now
    always regarded as an error (previously only checked if cds=True)
    and will trigger a warning (likely to become an exception in a
    future release).

    If **cds=True**, the start and stop codons are checked, and the start
    codon will be translated at methionine. The sequence must be an
    while number of codons.

    >>> _translate_str("ATGCCCTAG", table, cds=True)
    'MP'
    >>> _translate_str("AAACCCTAG", table, cds=True)
    Traceback (most recent call last):
       ...
    Bio.Data.CodonTable.TranslationError: First codon 'AAA' is not a start codon
    >>> _translate_str("ATGCCCTAGCCCTAG", table, cds=True)
    Traceback (most recent call last):
       ...
    Bio.Data.CodonTable.TranslationError: Extra in frame stop codon found.
    """

    if cds and len(sequence) % 3 != 0:
        raise CodonTable.TranslationError("Sequence length {0} is not a multiple of three".format(
            len(sequence)
        ))
    elif gap is not None and (not isinstance(gap, str) or len(gap) > 1):
        raise TypeError("Gap character should be a single character "
                        "string.")

    forward_table, getter, valid_letters = get_tables(table, to_stop=to_stop, gap=gap)

    sequence = sequence.upper()
    if not valid_letters.issuperset(set(sequence)):
        raise CodonTable.TranslationError("Invalid letters in the sequence: {}".format(
            set.difference(set(*sequence), valid_letters)
        ))

    amino_acids = getter(np.array(
        [sequence[start:start + 3] for start in range(0, len(sequence) - len(sequence) % 3, 3)], dtype="<U"))

    if cds and amino_acids[0] != "M":
        raise CodonTable.TranslationError(
            "First codon '{0}' is not a start codon".format(sequence[:3]))

    nones = np.where(amino_acids == None)[0]

    if nones.shape[0] > 0:
        assert pos_stop is not None
        amino_acids[nones] = pos_stop

    _stop_locations = np.where(amino_acids == stop_symbol)[0]
    found_stops = _stop_locations.shape[0]

    if cds and found_stops > 1:
        raise CodonTable.TranslationError("Extra in frame stop codon found.")
    elif cds and found_stops and _stop_locations[0] < len(amino_acids) - 1:
        raise CodonTable.TranslationError("Extra in frame stop codon found.")
    if to_stop and found_stops > 0:
        amino_acids = amino_acids[:_stop_locations[0]]

    return "".join(amino_acids)


# These classes do contain lots of things, it is correct like it is
# pylint: disable=too-many-instance-attributes
class BED12:

    """
    BED12 parsing class.
    """

    __valid_coding = {"True": True, "False": False, True: True, False: False}

    _attribute_pattern = re.compile(r"([^;]*)=([^$=]*)(?:;|$)")

    def __init__(self, *args: Union[str, list, tuple, GffLine],
                 fasta_index=None,
                 phase=None,
                 sequence=None,
                 transcriptomic=False,
                 max_regression=0,
                 start_adjustment=True,
                 coding=True,
                 lenient=False,
                 table=0,
                 logger=create_null_logger()):

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
        self.__parent = None
        self.transcriptomic = transcriptomic
        if phase is not None:
            self.phase = phase  # This will check the validity of the phase itself
        self.__strand = None
        self.__max_regression = 0
        # This value will be set when checking for the internal sequence
        # If >=1, i.e. at least one internal stop codon, the ORF is invalid
        self._internal_stop_codons = 0
        self.chrom = None
        self.__start = self.__end = self.__thick_start = self.__thick_end = 0
        self.name = ""
        self.score = 0
        self.strand = None
        self.rgb = ''
        self.__block_sizes = np.zeros(1, dtype=np.int64)
        self.__block_starts = np.zeros(1, dtype=np.int64)
        self.__block_count = 1
        self.__invalid = None
        self.invalid_reason = None
        self.fasta_length = None
        self.__in_index = True
        self.max_regression = max_regression
        self.start_adjustment = start_adjustment
        self.coding = coding
        if self.coding and self.phase is None:
            self.phase = 0
        self.__table = standard
        self.__table_index = 0
        self.table = table
        self.__lenient = lenient
        self.alias = None
        self.__logger = create_null_logger()
        self.logger = logger
        self.logger.debug("Set the basic properties for %s", self.chrom)

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
            if len(self._fields) in (12, 13):
                self.__set_values_from_fields()
                self.header = False
            else:
                self.header = True
                return
        elif isinstance(self._line, type(self)):  # Re-initialising with another object
            self.__set_values_from_bed12(args[0])
        elif isinstance(self._line, GffLine):
            if self._line.header is True:
                self.header = True
                return
            elif self.transcriptomic is False:
                error = "GFF lines can be used as BED12-equivalents only in a transcriptomic context."
                self.logger.error(error)
                raise TypeError(error)
            else:
                self.header = False
                if sequence:
                    fasta_length = len(sequence)
                elif fasta_index:
                    if isinstance(fasta_index, pysam.FastaFile):
                        fasta_length = fasta_index.get_reference_length(self._line.chrom)
                    elif isinstance(fasta_index, pyfaidx.Fasta):
                        sequence = fasta_index[self._line.chrom]
                        fasta_length = len(fasta_index[self._line.chrom])
                    else:
                        raise TypeError("Invalid FASTA index")
                else:
                    raise ValueError("Either a sequence or a FAI index are needed")
                self.__set_values_from_gff(fasta_length)

        elif not (isinstance(self._line, list) or isinstance(self._line, tuple)):
            raise TypeError("I need an ordered array, not {0}".format(type(self._line)))
        else:
            self._fields = self._line
            self.__set_values_from_fields()

        self.__check_validity(transcriptomic, fasta_index, sequence)

        if self.invalid and self.coding:
            self.logger.debug("%s cannot be coding as it is invalid (reason: %s)", self.chrom, self.invalid_reason)
            self.coding = False

        if self.coding and self.phase is None:
            self.phase = 0

    @property
    def is_transcript(self):
        """BED12 files are always transcripts for Mikado."""

        return True

    @property
    def is_gene(self):
        """BED12 files are never "genes" according to Mikado"""
        return False

    @property
    def source(self):
        return "bed12"

    @property
    def gene(self):
        return self.parent[0] if self.parent else None

    @property
    def parent(self):
        return self.__parent

    @property
    def table(self):
        return self.__table

    @table.setter
    def table(self, table):
        if table is None:
            self.__table = standard
            self.__table_index = 0
        elif isinstance(table, int):
            if table == 0:
                self.__table = standard
            else:
                self.__table = CodonTable.ambiguous_dna_by_id[table]
            self.__table_index = 0
        elif isinstance(table, str):
            self.__table = CodonTable.ambiguous_dna_by_name[table]
            self.__table_index = self.__table._codon_table.id
        elif isinstance(table, bytes):
            self.__table = CodonTable.ambiguous_dna_by_name[table.decode()]
            self.__table_index = self.__table._codon_table.id
        else:
            raise ValueError("Invalid table: {} (type: {})".format(
                    table, type(table)))
        return

    @parent.setter
    def parent(self, parent):
        if parent is not None and not isinstance(parent, str):
            raise TypeError(type(parent))
        self.__parent = [parent]

    def __getstate__(self):

        state = copy.deepcopy(dict((key, val) for key, val in self.__dict__.items()
                                   if key not in ("_BED12_table") and
                                   not isinstance(val, logging.Logger) and
                                   not isinstance(val, CodonTable.CodonTable)))

        return state

    def __setstate__(self, state):
        # del state["table"]
        self.__dict__.update(state)
        self.table = self.__table_index

    def _parse_attributes(self, attributes):

        """
        Private method that parses the last field of the GFF line.
        :return:
        """

        self.attribute_order = []

        infolist = self._attribute_pattern.findall(attributes.rstrip().rstrip(";"))

        for item in infolist:
            key, val = item
            if key.lower() in ("parent", "geneid"):
                self.parent = val
            elif "phase" in key.lower():
                self.phase = int(val)
                self.coding = True
            elif key.lower() == "coding":
                self.coding = self.__valid_coding.get(val, False)
                if self.transcriptomic is True:
                    self.phase = 0
            elif key.lower() == "alias":
                self.alias = val
            elif key.lower() == "id":
                self.name = val
            else:
                continue

    def __set_values_from_fields(self):

        """
        Private method that sets the correct values from the fields derived from the input line.
        :return:
        """
        self.chrom, self.start, self.end, \
            self.name, self.score, self.strand, \
            self.thick_start, self.thick_end, self.rgb, \
            self.block_count, block_sizes, block_starts = self._fields[:12]

        # Reduce memory usage
        intern(self.chrom)
        self.start = int(self.start) + 1
        self.end = int(self.end)
        try:
            self.score = float(self.score)
        except (ValueError, TypeError):
            self.score = None
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
        self._parse_attributes(self.name)
        if len(self._fields) == 13:
            self._parse_attributes(self._fields[-1])
        self.has_start_codon = None
        self.has_stop_codon = None
        self.start_codon = None
        self.stop_codon = None
        self.fasta_length = len(self)
        return

    def __set_values_from_bed12(self, line):

        self.__setstate__(line.__getstate__())
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

        del self.invalid

        if transcriptomic is True and self.coding is True:
            if not (fasta_index is not None or sequence is not None):
                self.logger.debug("No further check on the validity of %s as no sequence has been provided.",
                                  self.chrom)
                return

        if transcriptomic is True:
            self.has_start_codon = False
            self.has_stop_codon = False

        if transcriptomic is True and self.coding is True and (fasta_index is not None or sequence is not None):
            self.logger.debug("Starting to check the validity of %s", self.chrom)
            self.validity_checked = True
            if sequence is not None:
                self.fasta_length = len(sequence)
                if hasattr(sequence, "seq"):
                    sequence = str(sequence.seq)
                if not isinstance(sequence, str):
                    sequence = str(sequence)
            else:
                if self.id not in fasta_index:
                    self.logger.warning("%s not found in the index. Aborting the check, we will trust the ORF as-is.")
                    self.__in_index = False
                    return
                self.fasta_length = len(fasta_index[self.id])
                sequence = fasta_index[self.id]
                if hasattr(sequence, "seq"):
                    sequence = str(sequence.seq)
                if not isinstance(sequence, str):
                    sequence = str(sequence)

            assert isinstance(sequence, str)
            # Just double check that the sequence length is the same as what the BED would suggest
            if self.__is_invalid() is True:
                self.logger.debug("%s is invalid (%s)", self.chrom, self.invalid_reason)
                self.coding = False
                return

            if self.strand != "-":
                orf_sequence = sequence[
                               (self.thick_start - 1 if not self.phase else self.start + self.phase - 1):self.thick_end]
            else:
                orf_sequence = Seq.reverse_complement(
                    sequence[(self.thick_start - 1):(
                        self.thick_end if not self.phase else self.end - (3 - self.phase) % 3)])

            self.start_codon = str(orf_sequence)[:3].upper()
            self.stop_codon = str(orf_sequence[-3:]).upper()

            if self.start_codon in self.table.start_codons and (self.phase is None or self.phase == 0):
                self.logger.debug("Found start codon for %s. Setting phase to 0", self.chrom)
                self.has_start_codon = True
                self.phase = 0
            else:
                self.has_start_codon = False
                if self.start_adjustment is True:
                    self._adjust_start(sequence, orf_sequence)

            if self.stop_codon in self.table.stop_codons:
                self.has_stop_codon = True
            else:
                self.has_stop_codon = False
                if self.strand == "+" and self.end - self.thick_end < 3:
                    self.thick_end = self.end
                elif self.strand == "-" and self.thick_start - self.start < 3:
                    self.thick_start = 1

            self.logger.debug("%s with start codon (%s) and stop codon (%s). Valid: %s",
                              self.chrom, self.has_start_codon, self.has_stop_codon, not self.invalid)

            # Get only a proper multiple of three
            if self.__lenient is False:
                if self.strand != "-":
                    orf_sequence = sequence[
                                   (self.thick_start - 1 if not self.phase
                                    else self.start + self.phase - 1):self.thick_end]
                else:
                    orf_sequence = Seq.reverse_complement(
                        sequence[
                        (self.thick_start - 1):
                        (self.thick_end if not self.phase else self.end - self.phase)])

                last_pos = -3 - ((len(orf_sequence)) % 3)
                translated_seq = _translate_str(orf_sequence[:last_pos],
                                                table=self.table,
                                                gap='N')

                self._internal_stop_codons = str(translated_seq).count("*")
            del self.invalid
            if self.__is_invalid() is True:
                return

    def _adjust_start(self, sequence, orf_sequence):

        # self.logger.debug("Checking %s", self.chrom)
        assert len(orf_sequence) == (self.thick_end - self.thick_start + 1 - self.phase), (
                len(orf_sequence), (self.thick_end - self.thick_start + 1 - self.phase)
        )
        # Let's check UPstream first.
        # This means that we DO NOT have a starting Met and yet we are starting far upstream.
        # self.logger.debug("Starting to adjust the start of %s", self.chrom)
        if self.strand == "+" and self.thick_start > 3:
            for pos in range(self.thick_start, 3, -3):
                self.thick_start -= 3
                codon = sequence[pos - 4:pos - 1]
                is_start, is_stop = ((codon in self.table.start_codons),
                                     (codon in self.table.stop_codons))
                self.logger.debug("Checking pos %s (%s) for %s, start: %s; stop: %s",
                                  pos, codon, self.chrom, is_start, is_stop)
                if is_start:
                    # We have found a valid methionine.
                    break
                elif is_stop:
                    self.stop_codon = codon
                    self._internal_stop_codons = 1
                    assert self.invalid is True
                    self.logger.debug(
                        "Found in-frame stop codon for %s while expanding, stopping here. Invalid: %s (reason %s)",
                                        self.chrom, self.invalid, self.invalid_reason)
                    break
                continue

        elif self.strand == "-" and self.end - self.thick_end > 3:
            for pos in range(self.thick_end, self.end - 3, 3):
                self.thick_end += 3
                codon = Seq.reverse_complement(sequence[pos - 3:pos])
                is_start, is_stop = ((codon in self.table.start_codons),
                                     (codon in self.table.stop_codons))
                # self.logger.debug("Checking pos %s (%s) for %s, start: %s; stop: %s",
                #                   pos, codon, self.chrom, is_start, is_stop)
                if is_start:
                    # We have found a valid methionine.
                    self.logger.debug("Found correct start codon for %s while expanding, stopping here.",
                                      self.chrom)
                    break
                elif is_stop:
                    self.stop_codon = codon
                    self._internal_stop_codons = 1
                    assert self.invalid is True
                    self.logger.debug(
                        "Found in-frame stop codon for %s while expanding, stopping here. Invalid: %s (reason %s)",
                        self.chrom, self.invalid, self.invalid_reason)
                    break
        else:
            self.__regression(orf_sequence)

        if self.has_start_codon is False:
            # The validity will be automatically checked
            # self.logger.debug("Making adjustments in %s for missing start codon", self.chrom)
            if self.strand == "+":
                if self.thick_start - self.start <= 2:
                    new_phase = max(self.thick_start - self.start, 0)
                    self.phase = new_phase
                    self.thick_start = self.start
                else:
                    self.phase = 0
            else:
                if self.end - self.thick_end <= 2:
                    new_phase = max(self.end - self.thick_end, 0)
                    self.phase = new_phase
                    self.thick_end = self.end
                else:
                    self.phase = 0
        else:
            self.logger.debug("Setting phase of %s at 0 (end: %s; thick end: %s; thick start %s)",
                              self.chrom, self.end, self.thick_end, self.thick_start)
            self.phase = 0

        del self.invalid
        if self.invalid:
            self.logger.debug("%s is not coding after checking. Reason: %s", self.chrom, self.invalid_reason)
            self.coding = False

    def __regression(self, orf_sequence):
        self.logger.debug(
            "Starting the regression algorithm to find an internal start for %s (end: %s; thick start/end: %s, %s; phase %s)",
            self.chrom, self.end, self.thick_start, self.thick_end, self.phase)
        if self.strand != "-":
            # self.thick_start = self.phase + 3
            self.logger.debug("Starting to analyse %s; positions %s-%s",
                              self.chrom,
                              self.phase + 3,
                              self.phase + 3 + int(len(orf_sequence) * self.max_regression),
                              )
            for pos in range(self.phase + 3,
                             int(len(orf_sequence) * self.max_regression),
                             3):
                codon = orf_sequence[pos:pos + 3]
                # self.logger.debug("Testing position %s-%s (%s)", pos, pos + 3, codon)
                if codon in self.table.start_codons:
                    # Now we have to shift the start accordingly
                    self.has_start_codon = True
                    self.thick_start += pos
                    self.phase = 0
                    break
                else:
                    continue
            self.logger.debug("Final internal coords for %s: %s-%s", self.chrom, self.thick_start, self.thick_end)
        elif self.strand == "-":
            if self.end - self.thick_end < 3:
                self.phase = (3 - (self.end - self.thick_end) % 3) % 3
            self.logger.debug("Starting to analyse %s (phase %s); positions %s-%s",
                              self.chrom,
                              self.phase,
                              self.phase + 3,
                              self.phase + 3 + int(len(orf_sequence) * self.max_regression),
                              )
            for pos in range(self.phase + 3,
                             int(len(orf_sequence) * self.max_regression),
                             3):
                codon = orf_sequence[pos:pos + 3]
                # self.logger.debug("Testing position %s-%s (%s)", pos, pos + 3, codon)
                if codon in self.table.start_codons:
                    # Now we have to shift the start accordingly
                    self.has_start_codon = True
                    self.thick_end -= pos
                    self.phase = 0
                    break
            self.logger.debug("Final internal coords for %s: %s-%s", self.chrom, self.thick_start, self.thick_end)

    def __str__(self):

        if self.header is True:
            if self._line is not None:
                return self._line
            else:
                return "#"

        line = [self.chrom, self.start - 1, self.end]

        if self.transcriptomic is True:
            name = "ID={};coding={}".format(self.id, self.coding)
            if self.coding:
                name += ";phase={}".format(self.phase)
            if self.alias is not None and self.alias != self.id:
                name += ";alias={}".format(self.alias)

            line.append(name)
        else:
            line.append(self.name)

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

    def __len__(self):
        return self.end - self.start + 1

    def copy(self):

        return copy.deepcopy(self)

    def as_simple_dict(self):

        return {
            "chrom": self.chrom,
            "id": self.id,
            "start": self.start,
            "end": self.end,
            "name": self.name,
            "strand": self.strand,
            "thick_start": self.thick_start,
            "thick_end": self.thick_end,
            "score": self.score,
            "has_start_codon": self.has_start_codon,
            "has_stop_codon": self.has_stop_codon,
            "cds_len": self.cds_len,
            "phase": self.phase,
            "transcriptomic": self.transcriptomic,
        }

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

        if self.__invalid is None:
            self.__invalid = self.__is_invalid()

        return self.__invalid

    @invalid.deleter
    def invalid(self):
        self.__invalid = None


    def __is_invalid(self):

        if self._internal_stop_codons >= 1:
            self.invalid_reason = "{} internal stop codons found".format(self._internal_stop_codons)
            return True

        if self.fasta_length is None:
            self.fasta_length = len(self)

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

        if self.transcriptomic is True:
            if self.__in_index is False:
                self.invalid_reason = "{} not found in the index!".format(self.chrom)
                return True

            if len(self) != self.fasta_length:
                self.invalid_reason = "Fasta length != BED length: {0} vs. {1}".format(
                    self.fasta_length,
                    len(self)
                )
                return True

            if self.__lenient is True:
                pass
            else:
                if (self.cds_len - self.phase) % 3 != 0:
                    if self.strand == "+" and self.thick_end != self.end:
                        self.invalid_reason = "Invalid CDS length: {0} % 3 = {1} ({2}-{3}, {4})".format(
                            self.cds_len - self.phase,
                            (self.cds_len - self.phase) % 3,
                            self.thick_start, self.thick_end, self.phase)
                        return True
                    elif self.strand == "-" and self.thick_start != self.start:
                        self.invalid_reason = "Invalid CDS length: {0} % 3 = {1} ({2}-{3}, {4})".format(
                            self.cds_len - self.phase,
                            (self.cds_len - self.phase) % 3,
                            self.thick_start, self.thick_end, self.phase)
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
    def start(self):
        return self.__start

    @start.setter
    def start(self, value):
        try:
            value = int(value)
        except (ValueError, TypeError):
            raise ValueError("Start must be an integer, not {}! Value: {}".format(type(value), value))
        self.__start = value
        del self.invalid

    @start.deleter
    def start(self):
        self.__start = 0
        del self.invalid

    @property
    def end(self):
        return self.__end

    @end.setter
    def end(self, value):
        try:
            value = int(value)
        except (ValueError, TypeError):
            raise ValueError("End must be an integer, not {}! Value: {}".format(type(value), value))
        self.__end = value
        del self.invalid

    @end.deleter
    def end(self):
        self.__end = 0
        del self.invalid

    @property
    def thick_start(self):
        return self.__thick_start

    @thick_start.setter
    def thick_start(self, value):
        try:
            value = int(value)
        except (ValueError, TypeError):
            raise ValueError("Thick start must be an integer, not {}! Value: {}".format(type(value), value))
        self.__thick_start = value
        del self.invalid

    @thick_start.deleter
    def thick_start(self):
        self.__thick_start = 0
        del self.invalid

    @property
    def thick_end(self):
        return self.__thick_end

    @thick_end.setter
    def thick_end(self, value):
        try:
            value = int(value)
        except (ValueError, TypeError):
            raise ValueError("Thick end must be an integer, not {}! Value: {}".format(type(value), value))
        self.__thick_end = value
        del self.invalid

    @thick_end.deleter
    def thick_end(self):
        self.__thick_end = 0
        del self.invalid

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
        del self.invalid
        self.__phase = val

    @phase.deleter
    def phase(self):
        self.__phase = None
        del self.invalid

    @property
    def block_count(self):
        return self.__block_count

    @block_count.setter
    def block_count(self, value):
        try:
            value = int(value)
        except (ValueError, TypeError):
            raise ValueError("Block count must be an integer, not {}! Value: {}".format(type(value), value))
        self.__block_count = value
        del self.invalid

    @property
    def block_sizes(self):
        return self.__block_sizes

    @block_sizes.setter
    def block_sizes(self, sizes):
        sizes = np.array(sizes)
        if not issubclass(sizes.dtype.type, np.int64):
            raise TypeError("Block sizes should be integers!")
        self.__block_sizes = sizes
        del self.invalid

    @block_sizes.deleter
    def block_sizes(self):
        self.__block_sizes = np.zeros(1, dtype=np.int64)
        del self.invalid

    @property
    def block_starts(self):
        return self.__block_starts

    @block_starts.setter
    def block_starts(self, starts):
        starts = np.array(starts)
        if not issubclass(starts.dtype.type, np.int64):
            raise TypeError("Block sizes should be integers! Dtype: {}; array: {}".format(
                starts.dtype, starts
            ))
        self.__block_starts = starts
        del self.invalid

    @block_starts.deleter
    def block_starts(self):
        self.__block_starts = np.zeros(1, dtype=np.int64)
        del self.invalid

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

    def expand(self, sequence, upstream, downstream, expand_orf=False, logger=create_null_logger()):

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
            raise NotImplementedError("I can only expand ORFs on the sense strand")

        old_sequence = sequence[upstream:len(self) + upstream]
        assert len(old_sequence) + upstream + downstream == len(sequence)
        self.fasta_length = len(sequence)

        # I presume that the sequence is already in the right orientation
        old_start_pos = self.thick_start + self.phase - 1
        old_end_pos = self.thick_end - (self.thick_end - old_start_pos) % 3
        old_orf = old_sequence[old_start_pos:old_end_pos].upper()
        logger.debug("Old sequence of %s (%s bps): %s[...]%s", self.id, len(old_sequence),
                     old_sequence[:10], old_sequence[-10:])
        logger.debug("Old ORF of %s (%s bps, phase %s): %s[...]%s", self.id, len(old_orf), self.phase,
                     old_orf[:10], old_orf[-10:])
        assert len(old_orf) > 0, (old_start_pos, old_end_pos)
        assert len(old_orf) % 3 == 0, (old_start_pos, old_end_pos)

        old_pep = _translate_str(old_orf, self.table, gap="N")
        if "*" in old_pep and old_pep.find("*") < len(old_pep) - 1:
            logger.error("Stop codon found within the ORF of %s (pos %s of %s; phase %s). This is invalid!",
                         self.id, old_pep.find("*"), len(old_pep), self.phase)

        self.start_codon = old_orf[:3]
        self.stop_codon = old_orf[-3:]
        logger.debug("%s: start codon %s, old start %s (%s); stop codon %s, old stop %s (%s)",
                     self.name, self.start_codon, self.thick_start + self.phase,
                     (self.thick_start + self.phase + upstream),
                     self.stop_codon, self.thick_end, (self.thick_end + upstream))
        # Now expand
        self.end = len(sequence)
        self.thick_start += upstream
        self.thick_end += upstream
        start_codon = str(self.start_codon).upper()
        stop_codon = str(self.stop_codon).upper()
        self.has_start_codon = (start_codon in self.table.start_codons)
        self.has_stop_codon = (stop_codon in self.table.stop_codons)
        self.logger.debug("%s has start codon (%s): %s", self.chrom, start_codon, self.has_start_codon)
        self.logger.debug("%s has stop codon (%s): %s", self.chrom, stop_codon, self.has_stop_codon)
        if expand_orf is True and not (self.has_start_codon and self.has_stop_codon):
            if not self.has_start_codon:
                for pos in range(old_start_pos + upstream,
                                 0,
                                 -3):
                    codon = sequence[pos:pos + 3].upper()

                    self.thick_start = pos + 1
                    if codon in self.table.start_codons:
                        # self.thick_start = pos
                        self.start_codon = codon
                        self.__has_start = True
                        logger.debug("Position %d, codon %s. Start codon found.", pos, codon)
                        break
                if self.start_codon not in self.table.start_codons:
                    self.phase = (self.thick_start - 1) % 3
                    logger.debug("No start codon found for %s. Thick start %s, new phase: %s",
                                 self.id, self.thick_start, self.phase)
                    self.thick_start = 1
                else:
                    self.phase = 0
                    self.__has_start = True

            coding_seq = sequence[self.thick_start + self.phase - 1:self.end]
            if len(coding_seq) % 3 != 0:
                # Only get a multiple of three
                coding_seq = coding_seq[:-((len(coding_seq)) % 3)]
            prot_seq = _translate_str(coding_seq, table=self.table, gap="N")
            # print(coding_seq, prot_seq, self.table, sep="\n")
            # raise ValueError()
            if "*" in prot_seq:
                self.thick_end = self.thick_start + self.phase - 1 + (1 + prot_seq.find("*")) * 3
                self.stop_codon = coding_seq[prot_seq.find("*") * 3:(1 + prot_seq.find("*")) * 3].upper()
                self.__has_stop = True
                logger.debug("New stop codon for %s: %s", self.name, self.thick_end)

            if self.stop_codon not in self.table.stop_codons:
                logger.debug("No valid stop codon found for %s", self.name)
                self.thick_end = self.end

        self.block_sizes = [self.thick_end - self.thick_start]
        self.block_starts = [self.thick_start]
        return

    @property
    def blocks(self):

        """This will return the coordinates of the blocks, with a 1-offset (as in GFF3)"""

        # First thing: calculate where each start point will be
        starts = self.block_starts + self.start - 1
        _bstarts = starts + 1
        _bends = starts + self.block_sizes

        return list(zip(_bstarts, _bends))

    def to_transcriptomic(self, sequence=None, fasta_index=None, start_adjustment=False,
                          lenient=False, alias=None, coding=True):

        """This method will return a transcriptomic version of the BED12. If the object is already transcriptomic,
        it will return itself."""

        if self.transcriptomic is True:
            return self
        self.coding = coding

        # First six fields of a BED object

        # Now we have to calculate the thickStart
        # block_count = self.block_count

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
            bsizes = np.flip(self.block_sizes)
            tStart, tEnd = self.block_sizes.sum() - tEnd, self.block_sizes.sum() - tStart

        bstarts = np.concatenate([np.zeros(1, dtype=np.int64), bsizes[:-1].cumsum()])
        # bstarts = [0]
        # for bs in bsizes[:-1]:
        #     bstarts.append(bs + bstarts[-1])
        assert len(bstarts) == len(bsizes) == self.block_count, (bstarts, bsizes, self.block_count)

        if self.coding:
            new_name = "ID={};coding={};phase={}".format(self.name.split(";")[0],
                                                         self.coding,
                                                         self.phase if self.phase is not None else 0)
        else:
            new_name = "ID={};coding={}".format(self.name.split(";")[0], self.coding)

        if alias is not None:
            new_name += ";alias={}".format(alias)

        new = list((self.name.split(";")[0],
                    0,
                    self.block_sizes.sum(),
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
                    lenient=lenient,
                    start_adjustment=start_adjustment)
        # assert new.invalid is False
        assert isinstance(new, type(self)), type(new)
        return new

    @property
    def logger(self):
        """
        Property. It returns the logger instance attached to the class.
        :rtype : logging.Logger | None
        """

        return self.__logger

    @logger.setter
    def logger(self, logger):
        """Set a logger for the instance.
        :param logger: a Logger instance
        :type logger: logging.Logger | None
        """
        if logger is None:
            if self.__logger is None:
                logger = create_null_logger()
                self.__logger = logger
            else:
                pass
        else:
            assert isinstance(logger, logging.Logger), type(logger)
            self.__logger = logger
        self.__logger.propagate = False
    #
    # @logger.deleter
    # def logger(self):
    #     self.__logger = create_null_logger()


class Bed12Parser(Parser):
    """Parser class for a Bed12Parser file.
    It accepts optionally a fasta index which is used to
    determine whether an ORF has start/stop codons."""

    __annot_type__ = "bed12"

    def __init__(self, handle,
                 fasta_index=None,
                 transcriptomic=False,
                 max_regression=0,
                 start_adjustment=True,
                 is_gff=False,
                 coding=False,
                 logger=create_null_logger(),
                 table=0):
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
        self.coding = coding
        self.start_adjustment = start_adjustment
        self.logger = logger
        self.fasta_index = self.__set_fasta_index(fasta_index)
        self.__closed = False
        self.header = False
        self.__table = table
        self._is_bed12 = (not is_gff)

    @staticmethod
    def __set_fasta_index(fasta_index):
        if isinstance(fasta_index, dict):
            # check that this is a bona fide dictionary ...
            assert isinstance(
                fasta_index[random.choice(fasta_index.keys())],
                # fasta_index[numpy.random.choice(fasta_index.keys(), 1)],
                Bio.SeqRecord.SeqRecord)
        elif fasta_index is not None:
            if isinstance(fasta_index, (str, bytes)):
                if isinstance(fasta_index, bytes):
                    fasta_index = fasta_index.decode()
                assert os.path.exists(fasta_index)
                fasta_index = pysam.FastaFile(fasta_index)
            else:
                assert isinstance(fasta_index, pysam.FastaFile), type(fasta_index)
        return fasta_index

    def __iter__(self):
        return self

    def __next__(self, seq=None):

        if self._is_bed12 is True:
            return self.bed_next()
        else:
            return self.gff_next()

    def __getstate__(self):
        state = super().__getstate__()
        # Now let's remove the fasta index
        if state["fasta_index"] is not None:
            if isinstance(state["fasta_index"], pysam.FastaFile):
                state["fasta_index"] = state["fasta_index"].filename
        return state

    def __setstate__(self, state):
        fasta_index = state.pop("fasta_index", None)
        super().__setstate__(state)
        self.logger = create_null_logger()
        self.__set_fasta_index(fasta_index)

    def bed_next(self):
        """

        :return:
        """

        bed12 = None
        while bed12 is None:
            line = next(self._handle)
            try:
                bed12 = BED12(line,
                              fasta_index=self.fasta_index,
                              transcriptomic=self.transcriptomic,
                              max_regression=self._max_regression,
                              coding=self.coding,
                              table=self.__table,
                              logger=self.logger,
                              start_adjustment=self.start_adjustment)
            except Exception:
                error = "Invalid line for file {}, position {}:\n{}".format(
                    self.name, self._handle.tell(), line)
                raise ValueError(error)
        return bed12

    def gff_next(self):
        """

        :return:
        """

        bed12 = None
        while bed12 is None:
            line = next(self._handle)
            try:
                gff_line = GffLine(line)
            except Exception:
                error = "Invalid line for file {}, position {}:\n{}".format(
                    self.name, self._handle.tell(), line)
                raise ValueError(error)

            if gff_line.feature != "CDS":
                continue
            # Compatibility with BED12
            try:
                bed12 = BED12(gff_line,
                              fasta_index=self.fasta_index,
                              transcriptomic=self.transcriptomic,
                              max_regression=self._max_regression,
                              table=self.__table,
                              start_adjustment=self.start_adjustment,
                              logger=self.logger)
            except Exception:
                error = "Invalid line for file {}, position {}:\n{}".format(
                    self.name, self._handle.tell(), line)
                raise ValueError(error)
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

    @property
    def coding(self):
        if not hasattr(self, "__coding"):
            self.__coding = False
        return self.__coding

    @coding.setter
    def coding(self, coding):
        if coding not in (False, True):
            raise ValueError(coding)
        self.__coding = coding


class Bed12ParseWrapper(mp.Process):

    def __init__(self,
                 identifier=None,
                 rec_queue=None,
                 return_queue=None,
                 log_queue=None,
                 level="DEBUG",
                 fasta_index=None,
                 transcriptomic=False,
                 max_regression=0,
                 is_gff=False,
                 coding=False,
                 start_adjustment=True,
                 table=0):

        """
        :param send_queue:
        :type send_queue: mp.Queue
        :param return_queue:
        :type send_queue: mp.Queue
        :param kwargs:
        """

        super().__init__()
        if identifier is None:
            raise ValueError
        self.__identifier = identifier
        self.rec_queue = rec_queue
        self.return_queue = return_queue
        self.logging_queue = log_queue
        self.transcriptomic = transcriptomic
        self.__max_regression = 0
        self._max_regression = max_regression
        self.coding = coding
        self.start_adjustment = start_adjustment

        if isinstance(fasta_index, dict):
            # check that this is a bona fide dictionary ...
            assert isinstance(
                # fasta_index[numpy.random.choice(fasta_index.keys(), 1)],
                fasta_index[random.choice(fasta_index.keys())],
                Bio.SeqRecord.SeqRecord)
        elif fasta_index is not None:
            if isinstance(fasta_index, (str, bytes)):
                if isinstance(fasta_index, bytes):
                    fasta_index = fasta_index.decode()
                assert os.path.exists(fasta_index)
                # fasta_index = pysam.FastaFile(fasta_index)
                fasta_index = pyfaidx.Fasta(fasta_index)
            else:
                assert isinstance(fasta_index, pysam.FastaFile), type(fasta_index)
        self._level = level
        self.fasta_index = fasta_index
        self.__closed = False
        self.header = False
        self.__table = table
        self._is_bed12 = (not is_gff)

    def bed_next(self, line, sequence=None):
        """

        :return:
        """

        try:
            bed12 = BED12(line,
                          logger=self.logger,
                          sequence=sequence,
                          transcriptomic=self.transcriptomic,
                          max_regression=self._max_regression,
                          start_adjustment=self.start_adjustment,
                          coding=self.coding,
                          table=self.__table)
        except Exception:
            raise ValueError("Invalid line: {}".format(line))
        return bed12

    def gff_next(self, line, sequence):
        """

        :return:
        """

        try:
            line = GffLine(line)
        except Exception:
            error = "Invalid line:\n{}".format(line)
            raise ValueError(error)

        if line.feature != "CDS":
            return None
            # Compatibility with BED12
        bed12 = BED12(line,
                      logger=self.logger,
                      sequence=sequence,
                      transcriptomic=self.transcriptomic,
                      max_regression=self._max_regression,
                      start_adjustment=self.start_adjustment,
                      table=self.__table)
        # raise NotImplementedError("Still working on this!")
        return bed12

    def run(self, *args, **kwargs):
        self.handler = logging_handlers.QueueHandler(self.logging_queue)
        self.logger = logging.getLogger(self.name)
        self.logger.addHandler(self.handler)
        self.logger.setLevel(self._level)
        self.logger.propagate = False

        self.logger.info("Started %s", self.__identifier)
        if self.rec_queue is None:
            self.return_queue.put(b"FINISHED")
            raise ValueError
        while True:
            if self.rec_queue.empty():
                sleep(0.1)
                continue
            line = self.rec_queue.get()
            if line in ("EXIT", b"EXIT"):
                self.rec_queue.put(b"EXIT")
                self.return_queue.put(b"FINISHED")
                break
            try:
                num, line, seq = line
                if seq is not None:
                    seq = zlib.decompress(seq).decode()
                if not self._is_bed12:
                    row = self.gff_next(line, seq)
                else:
                    row = self.bed_next(line, seq)

                if not row or row.header is True:
                    continue
                if row.invalid is True:
                    self.logger.warning("Invalid entry, reason: %s\n%s",
                                        row.invalid_reason,
                                        row)
                    continue
                # self.cache[num] = 
                self.return_queue.put((num, msgpack.dumps(row.as_simple_dict())))
            except AttributeError:
                pass
            except ValueError:
                raise ValueError(line)
