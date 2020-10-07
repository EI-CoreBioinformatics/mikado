# coding: utf-8

"""
Pretty basic class that defines a reference gene with its transcripts.
Minimal checks.
"""

import re
import copy
import logging
import operator
from sys import intern
from .transcript import Transcript
from .transcriptcomputer import TranscriptComputer
from ..exceptions import InvalidTranscript, InvalidCDS
from ..parsers.GFF import GffLine
from ..parsers.GTF import GtfLine
from ..utilities.log_utils import create_null_logger


class Gene:

    """
    :param transcr: a transcript used to initialize the container.
    :param gid:Id of the gene.
    :param logger: an optional Logger from the logging module.
    """

    __name__ = "gene"

    def __init__(self, transcr: [None, Transcript], gid=None, logger=create_null_logger(),
                 only_coding=False, use_computer=False):

        self.transcripts = dict()
        self.__logger = None
        self.logger = logger
        self.__introns = None
        self.exception_message = ''
        self.chrom, self.source, self.start, self.end, self.strand = [None] * 5
        self.only_coding = only_coding
        self.coding_transcripts = set()
        self.id = None
        self.attributes = dict()
        self.feature = "gene"
        self.__from_gene = False
        self.__use_computer = use_computer

        if transcr is not None:
            if isinstance(transcr, Transcript):
                if isinstance(transcr, TranscriptComputer):
                    self.__use_computer = True
                else:
                    self.__use_computer = False

                self.transcripts[transcr.id] = transcr
                if transcr.parent:
                    self.id = transcr.parent[0]
                else:
                    self.logger.debug("No gene ID found for %s, creating a mock one.", transcr.id)
                    transcr.parent = f"{transcr.id}.gene"
                    self.id = transcr.parent[0]
                self.transcripts[transcr.id] = transcr
            elif isinstance(transcr, GffLine):
                if transcr.is_gene is True:
                    self.__from_gene = True
                    self.id = transcr.id
                    self.attributes = transcr.attributes.copy()
                    self.feature = transcr.feature
                elif transcr.is_exon is True:
                    self.__from_gene = False
                    if transcr.parent:
                        self.id = transcr.parent[0]
                    else:
                        self.id = transcr.id
                    self.attributes = transcr.attributes.copy()
                    self.feature = "gene"
                    if self._use_computer is True:
                        ntranscr = TranscriptComputer(transcr,
                                                      trust_orf=True,
                                                      logger=self.logger,
                                                      accept_undefined_multi=True,
                                                      source=transcr.source,
                                                      is_reference=True)
                    else:
                        ntranscr = Transcript(transcr,
                                              trust_orf=True,
                                              logger=self.logger,
                                              accept_undefined_multi=True,
                                              source=transcr.source,
                                              is_reference=True)
                    ntranscr.add_exon(transcr)
                    self.add(ntranscr)
                    self.logger.debug("New transcript for %s: %s", self.id, ntranscr.id)
            elif isinstance(transcr, GtfLine):
                self.id = transcr.gene

            self.chrom, self.source, self.start, self.end, self.strand = (transcr.chrom,
                                                                          transcr.source,
                                                                          transcr.start,
                                                                          transcr.end,
                                                                          transcr.strand)
        if gid is not None:
            self.id = gid
        # Internalize in memory for less memory usage
        [intern(str(_)) for _ in [self.chrom, self.source, self.id]
         if _ is not None]

    def __contains__(self, item):

        if isinstance(item, (str, bytes)):
            return item in self.transcripts
        elif isinstance(item, Transcript):
            return item in self.transcripts.values()

    def keys(self):
        return self.transcripts.keys()

    @property
    def logger(self):

        """
        Logger instance for the class.
        :rtype : logging.Logger
        """
        return self.__logger

    @logger.setter
    def logger(self, logger):
        """Set a logger for the instance.
        :param logger
        :type logger: logging.Logger | None
        """
        if isinstance(logger, logging.Logger):
            self.__logger = logger
        elif logger is None:
            name = "gene_{0}".format(self.id if self.id else "generic")
            self.__logger = create_null_logger(name)
        else:
            raise TypeError("Invalid object for logger: {0}, (type {1})".format(
                logger, type(logger)))

        for tid in self.transcripts:
            self.transcripts[tid].logger = logger

    @logger.deleter
    def logger(self):
        """
        Destroyer for the logger. It sets the internal __logger attribute to None.
        """
        self.__logger = None

    def add(self, transcr: Transcript):
        """
        This method adds a transcript to the storage.
        :param transcr: the transcript to be added.
        """

        self.transcripts[transcr.id] = transcr
        if self.chrom is None:
            self.chrom = transcr.chrom
        elif self.chrom != transcr.chrom:
            raise AssertionError("Discrepant chromosome for gene {0} and transcript {1}".format(
                    self.id, transcr.id
                ))

        if transcr.strand != self.strand:
            if self.strand is None:
                self.strand = transcr.strand
            elif transcr.strand is None:
                transcr.strand = self.strand
            else:
                raise AssertionError("Discrepant strands for gene {0} ({2}) and transcript {1} ({3})".format(
                    self.id, transcr.id, self.strand, transcr.strand
                ))

        transcr.logger = self.logger

    _prot_pattern = re.compile("-Protein$")

    def add_exon(self, row):
        """

        :param row:
        :type row: (GtfLine | GffLine)
        :return:
        """

        found_tids = set()
        for parent in (_ for _ in row.parent if _ in self.transcripts):
            found_tids.add(parent)
            self.transcripts[parent].add_exon(row)

        for parent in (_ for _ in row.parent if _ not in self.transcripts):
            found = False
            if self._prot_pattern.search(parent) and self._prot_pattern.sub(r"", parent):
                continue
            
            for tid in self.transcripts:
                if parent in self.transcripts[tid].derived_children:
                    found = True
                    if tid not in found_tids:
                        self.transcripts[tid].add_exon(row)
                    break
            if not found:
                if self._use_computer is True:
                    self.transcripts[parent] = TranscriptComputer(row,
                                                                  trust_orf=True,
                                                                  logger=self.logger,
                                                                  accept_undefined_multi=True,
                                                                  source=row.source,
                                                                  is_reference=True,
                                                                  )
                else:
                    self.transcripts[parent] = Transcript(row,
                                                          trust_orf=True,
                                                          logger=self.logger,
                                                          accept_undefined_multi=True,
                                                          source=row.source,
                                                          is_reference=True,
                                                          )
                self.transcripts[parent].parent = self.id

        if row.id in self.transcripts:
            found_tids.add(row.id)
            self.transcripts[row.id].add_exon(row)

        self.logger.debug("Found transcripts: %s", found_tids)

    def __getitem__(self, tid: str) -> Transcript:
        return self.transcripts[tid]

    def finalize(self, exclude_utr=False):
        """
        This method will finalize the container by checking the consistency of all the
        transcripts and eventually removing incorrect ones.

        :param exclude_utr: boolean flag
        :return:
        """

        to_remove = set()
        for tid in self.transcripts:
            try:
                self.transcripts[tid].finalize()
                if self.only_coding is True and self.transcripts[tid].selected_cds_length == 0:
                    to_remove.add(tid)
                if self.transcripts[tid].selected_cds_length > 0:
                    self.coding_transcripts.add(tid)
                if exclude_utr is True:
                    self.transcripts[tid].remove_utrs()
            except InvalidCDS:
                self.transcripts[tid].strip_cds()
            except InvalidTranscript as err:
                self.exception_message += "{0}\n".format(err)
                to_remove.add(tid)
            except Exception as err:
                self.exception_message += "Error in gene {} for transcript {}".format(self.id, tid)
                self.exception_message += "{0}\n".format(err)
                self.logger.exception(self.exception_message)
                raise

        for k in to_remove:
            del self.transcripts[k]

        if len(self.transcripts) > 0:
            if self.source is None:
                self.source = set([_.source for _ in self.transcripts.values()]).pop()
            __new_start = min(_.start for _ in self)

            if __new_start != self.start:
                if self.__from_gene is True:
                    self.logger.warning("Resetting the start for %s from %s to %d",
                                        self.id, self.start, __new_start)
                else:
                    self.logger.debug("Resetting the start for %s from %s to %d",
                                      self.id, self.start, __new_start)

                self.start = __new_start

            __new_end = max(_.end for _ in self)
            if __new_end != self.end:
                if self.__from_gene is True:
                    self.logger.warning("Resetting the end for %s from %s to %d",
                                        self.id, self.end, __new_end)
                else:
                    self.logger.debug("Resetting the end for %s from %s to %d",
                                      self.id, self.end, __new_end)
                self.end = __new_end
        if self.exception_message:
            self.logger.exception(self.exception_message)
            self.exception_message = ""

    def as_dict(self, remove_attributes=True):

        """
        Method to transform the gene object into a JSON-friendly representation.
        :return:
        """

        state = dict()
        for key in ["chrom", "source", "start", "end", "strand", "id"]:
            state[key] = getattr(self, key)

        state["use_computer"] = self._use_computer
        state["transcripts"] = dict.fromkeys(self.transcripts.keys())

        for tid in state["transcripts"]:
            try:
                state["transcripts"][tid] = self.transcripts[tid].as_dict(
                    remove_attributes=remove_attributes)
            except TypeError as exc:
                self.logger.error(tid, repr(self.transcripts[tid]))
                self.logger.exception(exc)
                raise TypeError(exc)

        return state

    def load_dict(self, state, exclude_utr=False, protein_coding=False, trust_orf=False):

        for key in ["chrom", "source", "start", "end", "strand", "id"]:
            setattr(self, key, state[key])

        self.__use_computer = state["use_computer"]
        for tid, tvalues in state["transcripts"].items():
            if self._use_computer is True:
                transcript = TranscriptComputer(logger=self.logger)
            else:
                transcript = Transcript(logger=self.logger)
            transcript.load_dict(tvalues, trust_orf=trust_orf)
            transcript.finalize()
            if protein_coding is True and transcript.is_coding is False:
                self.logger.debug("{0} is non coding ({1}, {2})".format(
                    transcript.id,
                    transcript.combined_cds,
                    transcript.segments))
                continue
            if exclude_utr is True:
                has_utrs = (transcript.utr_length > 0)
                transcript.remove_utrs()
                if has_utrs is True and (transcript.utr_length > 0):
                    raise AssertionError("Failed to remove the UTRs!")
            self.transcripts[tid] = transcript

        self.chrom = intern(self.chrom)
        if self.source:
            self.source = intern(self.source)
        self.id = intern(self.id)

        return

    def remove(self, tid: str):
        """

        :param tid: name of the transcript to remove.

        This method will remove a transcript from the container, and recalculate the
         necessary instance attributes.

        """

        del self.transcripts[tid]
        if len(self.transcripts) == 0:
            self.end = None
            self.start = None
            self.chrom = None
        else:
            self.start = min(self.transcripts[tid].start for tid in self.transcripts)
            self.end = max(self.transcripts[tid].end for tid in self.transcripts)

    def __repr__(self):
        return " ".join(self.transcripts.keys())

    def __str__(self):
        return self.format("gff3")

    def __iter__(self) -> Transcript:
        """Iterate over the transcripts attached to the gene."""
        return iter(self.transcripts.values())

    def __len__(self) -> int:
        return len(self.transcripts)

    def __getstate__(self):

        logger = self.logger
        del self.logger
        state = self.__dict__.copy()
        self.logger = logger
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self.logger = None

    def __lt__(self, other):
        if self.chrom != other.chrom:
            return self.chrom < other.chrom
        else:
            if self.start != other.start:
                return self.start < other.start
            elif self.end != other.end:
                return self.end < other.end
            else:
                if self.strand is not None and other.strand is not None:
                    return self.strand < other.strand
                elif self.strand is None and other.strand is not None:
                    return False
                elif self.strand is not None and other.strand is None:
                    return True
                else:
                    return False

    def __eq__(self, other):
        if self.chrom == other.chrom and self.start == other.start and \
                self.end == other.end and self.strand == other.strand:
            return True
        return False

    def format(self, format_name, transcriptomic=False):

        if format_name not in ("gff", "gtf", "gff3", "bed12", "bed"):
            raise ValueError(
                "Invalid format: {0}. Accepted formats: bed/bed12 (equivalent), gff/gff3 (equivalent), gtf".format(
                    format_name))

        self.finalize()  # Necessary to sort the exons
        lines = []
        if format_name in ("gff", "gff3") and transcriptomic is False:
            line = GffLine(None)
            for attr in ["chrom",
                         "source",
                         "start",
                         "end",
                         "strand"]:
                setattr(line, attr, getattr(self, attr))
            line.feature = self.feature
            line.attributes = self.attributes.copy()
            line.id = self.id
            assert line.id is not None
            lines.append(str(line))

        for tid, transcript in sorted(self.transcripts.items(), key=operator.itemgetter(1)):
            lines.append(transcript.format(format_name, transcriptomic=transcriptomic))

        if format_name in ("gff", "gff3"):
            lines.append("###")
        return "\n".join(lines)

    def copy(self):

        """Method to return a deep copy of the gene."""

        logger = self.logger
        del self.logger
        state = copy.deepcopy(self)
        self.logger = logger
        return state

    @property
    def monoexonic(self):
        """
        Boolean property. False if at least one transcript is multiexonic,
        True otherwise.
        :return: bool
        :rtype: bool
        """

        return all(transcript.monoexonic is True for transcript in self.transcripts.values())

    @property
    def introns(self):
        """
        It returns the set of all introns in the container.
        :rtype : set
        """

        return set.union(*[self.transcripts[tid].introns for tid in self.transcripts])

    @property
    def exons(self):
        """
        It returns the set of all exons in the container.
        :rtype : set
        """
        return set.union(*[set(self.transcripts[tid].exons) for tid in self.transcripts])

    @property
    def has_monoexonic(self):
        """
        True if any of the transcripts is monoexonic.
        :rtype : bool
        """
        return any(len(self.transcripts[tid].introns) == 0 for tid in self.transcripts.keys())

    # @property
    # def monoexonic(self):
    #     return all(len(self.transcripts[tid].introns) == 0 for tid in self.transcripts.keys())

    @property
    def num_transcripts(self):
        """
        Number of transcripts.
        :rtype : int
        """
        return len(self.transcripts)

    @property
    def num_coding_transcripts(self):
        """
        Number of coding transcripts
        :return:
        """

        return len([_ for _ in self.transcripts if self.transcripts[_].is_coding is True])

    @property
    def is_coding(self):
        """
        Property. It evaluates to True if at least one transcript is coding, False otherwise.
        """
        return any(self.transcripts[tid].selected_cds_length > 0 for tid in self.transcripts.keys())

    @property
    def _use_computer(self):
        return self.__use_computer