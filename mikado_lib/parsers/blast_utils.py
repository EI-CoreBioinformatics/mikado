#!/usr/bin/env python3

"""
This module contains generic-purpose utilities to deal with BLAST XML files.
"""

import os
import subprocess
import gzip
import operator
import multiprocessing
import io
import collections
import time
import threading
import queue
import logging
from mikado_lib.parsers import HeaderError

__author__ = 'Luca Venturini'


def create_opener(filename):

    """
    Function to create the appropriate opener for a BLAST file.
    If a handle is given instead of a filename, the function returns the input immediately.

    :param filename
    :return:
    """

    if isinstance(filename, (gzip.GzipFile, io.TextIOWrapper)):
        return filename
    elif not isinstance(filename, str) or not os.path.exists(filename):
        raise OSError("Non-existent file: {0}".format(filename))

    if filename.endswith(".gz"):
        if filename.endswith(".xml.gz"):
            return gzip.open(filename, "rt")
        elif filename.endswith(".asn.gz"):
            # I cannot seem to make it work with gzip.open
            zcat = subprocess.Popen(["zcat", filename], shell=False,
                                    stdout=subprocess.PIPE)
            blast_formatter = subprocess.Popen(
                ['blast_formatter', '-outfmt', '5', '-archive', '-'],
                shell=False, stdin=zcat.stdout, stdout=subprocess.PIPE)
            return io.TextIOWrapper(blast_formatter.stdout, encoding="UTF-8")
    elif filename.endswith(".xml"):
        return open(filename)
    else:
        raise ValueError("Unrecognized file format: %s", filename)


def check_beginning(handle, filename, previous_header):

    """
    Static method to check that the beginning of the XML file is actually correct.

    :return
    """

    exc = None
    header = []
    try:
        first_line = next(handle)
        first_line = first_line.strip()
        if first_line != '<?xml version="1.0"?>':
            exc = ValueError("Invalid header for {0}!\n\t{1}".format(filename, first_line))
            raise exc
        second_line = next(handle)
        second_line = second_line.strip()
        valid_schemas = [
            " ".join(['<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN"',
                      '"http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">']),
            '<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "NCBI_BlastOutput.dtd">'
        ]

        if second_line not in valid_schemas:
            exc = ValueError("Invalid XML type for {0}!\n\t{1}".format(filename, second_line))
            raise exc
        header = [first_line, second_line]
    except StopIteration:
        exc = OSError("Empty file: {0}".format(filename))
    except ValueError:
        pass

    if exc is not None:
        return handle, header, exc

    while True:
        line = next(handle)

        if "<Iteration>" in line:
            break
        line = line.rstrip()
        if not line:
            exc = HeaderError("Invalid header for {0}:\n\n{1}".format(
                filename,
                "\n".join(header)
            ))
            break
        if len(header) > 10**3:
            exc = HeaderError("Abnormally long header ({0}) for {1}:\n\n{2}".format(
                len(header),
                filename,
                "\n".join(header)
            ))
            break
        header.append(line)

    if not any(iter(True if "BlastOutput" in x else False for x in header)):
        exc = HeaderError("Invalid header for {0}:\n\n{1}".format(filename, "\n".join(header)))

    if previous_header is not None:
        checker = [header_line for header_line in header if
                   "BlastOutput_query" not in header_line]
        previous_header = [header_line for header_line in previous_header if
                           "BlastOutput_query" not in header_line]
        if checker != previous_header:
            exc = HeaderError("BLAST XML header does not match for {0}".format(
                filename))

    return handle, header, exc


def merge(intervals: [(int, int)]):
    """
    :param intervals: a list of integer duplexes
    :type intervals: list
    This function is used to merge together intervals, which have to be supplied as a list
    of duplexes - (start,stop). The function will then merge together overlapping tuples and
    return a list of non-overlapping tuples.
    If the list is composed by only one element, the function returns immediately.
    """

    # Assume tuple of the form (start,end)
    # And return 0- and 1-length intervals
    new_intervals = []
    for interval in intervals:
        new_intervals.append(tuple(sorted(interval)))

    intervals = new_intervals[:]
    if len(intervals) < 2:
        return intervals

    # Sort according to start, end
    intervals = sorted(intervals, key=operator.itemgetter(0, 1))
    final_list = [intervals[0]]

    for start, end in intervals[1:]:
        if start > final_list[-1][1]:
            final_list.append(tuple([start, end]))
        elif end > final_list[-1][1]:
            final_list[-1] = tuple([final_list[-1][0], end])
    return final_list


# pylint: disable=no-member
class _Merger(multiprocessing.Process):

    """
    This private class acts as a background process behind the XMLMerger class.
    This allows the XMLMerger class to appear like a normal read-only
    file-like objects, allowing compatibility with parsers such as
    e.g. Bio.Blast.NCBIXML.parse
    """

    def __init__(self, filenames, header, other_queue, logger=None):
        # pylint: disable=no-member
        multiprocessing.Process.__init__(self)
        # pylint: enable=no-member
        self.queue = other_queue
        self.filenames = filenames
        self.header = header
        self.logger = logger

    def run(self):
        """
        Implementation of the "run" method of the Process mother class.
        During the running, _Merger will perform the following:

        - check that the XML header is compatible
        - check that the file ends correctly
        - append the new lines to the stream

        WARNING: as it is impossible to look for the end of a
        gzipped file or of a stream, we have to keep all lines in memory
        to ascertain that the file we are trying to merge is not corrupt.
        This makes unfeasible to use the current implementation for
        merging large XML files.
        """

        # self.logger.info("Merger running")
        print_header = True

        for filename in self.filenames:
            self.logger.debug("Begun %s", filename)
            if print_header is True:
                # self.logger.info("Printing header")
                self.queue.put_nowait(self.header)
                print_header = False
            try:
                handle = create_opener(filename)
            except OSError as exc:
                self.logger.exception(exc)
                continue
            except ValueError as exc:
                self.logger.exception(exc)
                continue

            handle, _, exc = check_beginning(handle,
                                             filename,
                                             self.header)
            if exc is not None:
                self.logger.exception(exc)
                self.logger.error("Skipped %s", filename)

                continue
            self.logger.debug("Finished parsing header for %s", filename)
            lines = collections.deque()
            lines.append("<Iteration>")
            bo_found = False
            for line in handle:
                if line.strip() == "":
                    continue
                if "BlastOutput" in line:
                    bo_found = True
                    break
                lines.append(line.rstrip())

            if bo_found is False:
                exc = ValueError("{0} is an invalid XML file".format(filename))
                self.logger.exception(exc)
                continue

            self.logger.debug("Finished parsing lines for %s", filename)
            self.queue.put_nowait(lines)
            self.logger.debug("Sent %d lines for %s", len(lines), filename)

        # We HAVE  to wait some seconds, otherwise the XML parser
        # might miss the end of the file.
        time.sleep(5)
        self.queue.put(["</BlastOutput_iterations>\n</BlastOutput>"])

        self.queue.put("Finished")
        return
# pylint: enable=no-member


class XMLMerger(threading.Thread):

    """
    This class has the purpose of merging on the fly multiple BLAST alignment
    files, be they in XML, XML.gz, or ASN.gz format. It uses the _Merger
    private class as a background process which bothers itself with the real
    work, while this class provides a file-like interface for external
    applications such Bio.Blast.NCBIXML.parse.
    """

    logger = logging.getLogger("blast_merger")
    logger.propagate = False
    logger.setLevel(logging.INFO)
    __stream_handler = logging.StreamHandler()
    __formatter = logging.Formatter("{asctime} - {name} - {levelname} - {message}", style='{')
    __stream_handler.setFormatter(__formatter)
    logger.addHandler(__stream_handler)

    def __init__(self, filenames, log_level=logging.WARNING, log=None):

        threading.Thread.__init__(self)
        self.logger.setLevel(log_level)
        if log is not None:
            self.file_handler = logging.FileHandler(log, "w")
            self.file_handler.setFormatter(self.__formatter)
            self.logger.addHandler(self.file_handler)
            self.logger.removeHandler(self.__stream_handler)

        self.__filenames = collections.deque(filenames)
        self.__event = threading.Event()
        # pylint: disable=no-member
        manager = multiprocessing.Manager()
        # pylint: enable=no-member
        self.__queue = manager.Queue()
        self.lines = collections.deque()
        self.started = False
        self.finished = False
        if len(self.__filenames) > 0:
            while True:
                _, header, exc = check_beginning(
                    create_opener(self.__filenames[0]),
                    self.__filenames[0],
                    None)
                if exc is not None:
                    self.logger.exception(exc)
                    _ = self.__filenames.popleft()
                    continue
                self.header = header
                break
        self.logger.debug("Header has %d lines", len(self.header))
        if len(self.__filenames) == 0:
            raise IndexError("No files left!")

        self.merger = _Merger(self.__filenames, self.header, self.__queue, logger=self.logger)
        self.start()

    def run(self):
        # pylint: disable=no-member
        self.merger.start()
        # pylint: enable=no-member

        self.logger.info("Reader started")
        while True:
            try:
                lines = self.__queue.get_nowait()
            except queue.Empty:
                self.logger.debug("Queue was empty, waiting")
                time.sleep(0.01)
                continue
            if lines == "Finished":
                break
            self.logger.debug("Received %d lines", len(lines))
            self.lines.extend(lines)

        self.finished = True
        # pylint: disable=no-member
        self.merger.join()
        # pylint: enable=no-member
        return

    def read(self, size=None):
        """
        This method allows the XMLMerger class to act, for all purposes,
        like a file-like interface. The bytes are read from the queue,
        delivered by the background _Merger instance.
        :param size: optional parameter to indicate how much we want to read
        from the file-like interface.
        :return:
        """

        total = ""
        while len(self.lines) == 0:
            if self.finished is True:
                break
        if len(self.lines) == 0:
            raise StopIteration
        if size is None:
            while self.finished is False:
                time.sleep(0.1)
            total = self.lines
        else:
            previous_val = None
            import sys
            while sys.getsizeof(total) < size and len(self.lines) > 0:
                val = self.lines.popleft()
                if val == "</BlastOutput_iterations>\n</BlastOutput>":
                    self.logger.info("Received termination lines")
                    if previous_val is not None:
                        self.logger.info("Previous val: %s",
                                         previous_val[-2000:])
                total += val
                if previous_val is not None:
                    previous_val += val
                else:
                    previous_val = val

        total = total.rstrip()
        return total

    def __next__(self):

        while len(self.lines) == 0:
            if self.finished is True:
                break
        if len(self.lines) == 0:
            # Thread is finished and deque is exhausted. Finished
            raise StopIteration
        return self.lines.popleft()

    def __iter__(self):
        return self

    def __exit__(self, *args):
        _ = args
        self.join()

    def __enter__(self):
        pass
