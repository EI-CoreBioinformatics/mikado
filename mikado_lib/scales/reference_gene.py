# coding: utf-8

"""
Pretty basic class that defines a reference gene with its transcripts.
Minimal checks.
"""

import logging
from mikado_lib.loci_objects import transcript
from mikado_lib.exceptions import InvalidTranscript, InvalidCDS


class Gene:
    
    """
    :param tr: a transcript used to initialize the container.
    :param gid:Id of the gene.
    :param logger: an optional Logger from the logging module.
    """

    def __init__(self, tr: transcript, gid=None, logger=None):

        self.chrom, self.start, self.end, self.strand = tr.chrom, tr.start, tr.end, tr.strand
        self.id = gid
        self.transcripts = dict()
        self.transcripts[tr.id] = tr
        self.logger = logging.NullHandler()
        self.set_logger(logger)
        self.exception_message = ''
    
    def set_logger(self, logger):
        """
        :param logger: a Logger instance.
        :type logger: None | logging.Logger

        """
        if logger is None:
            self.logger = logging.getLogger("null_logger")
            handler = logging.NullHandler
            self.logger.addHandler(handler)
        else:
            self.logger = logger
        for tid in self.transcripts:
            self.transcripts[tid].logger = logger

    def add(self, tr: transcript):
        """
        :param tr:

        This method adds a transcript to the Locus.
        """
        self.start = min(self.start, tr.start)
        self.end = max(self.end, tr.end)
        self.transcripts[tr.id] = tr
        assert self.strand == tr.strand
        
    def __getitem__(self, tid: str) -> transcript:
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
                if exclude_utr is True:
                    self.transcripts[tid].remove_utrs()
            except InvalidCDS:
                self.transcripts[tid].strip_cds()
            except InvalidTranscript as err:
                self.exception_message += "{0}\n".format(err)
                to_remove.add(tid)
            except Exception as err:
                print(err)
                raise
        for k in to_remove:
            del self.transcripts[k]
    
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
        self.start = min(self.transcripts[tid].start for tid in self.transcripts)
        self.end = max(self.transcripts[tid].end for tid in self.transcripts)
    
    def __str__(self):
        return " ".join(self.transcripts.keys())
    
    def __iter__(self) -> transcript:
        """Iterate over the transcripts attached to the gene."""
        return iter(self.transcripts.values())
    
    def __len__(self) -> int:
        return len(self.transcripts)
    
    def __getstate__(self):
        self.logger = logging.NullHandler()
        
    def __setstate__(self, state):
        self.__dict__.update(state)
        self.set_logger(None)

    def __lt__(self, other):
        if self.chrom != other.chrom:
            return self.chrom < other.chrom
        else:
            if self.start != other.start:
                return self.start < other.start
            elif self.end != other.end:
                return self.end < other.end
            else:
                return self.strand < other.strand

    def __eq__(self, other):
        if self.chrom == other.chrom and self.start == other.start and \
                self.end == other.end and self.strand == other.strand:
            return True
        return False

