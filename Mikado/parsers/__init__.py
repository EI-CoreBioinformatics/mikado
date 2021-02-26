#!/usr/bin/env python3
# coding: utf_8

"""
    This module defines the iterators that will parse BED12, GTF, GFF files.
"""

import io
import abc
import gzip
import bz2
from functools import partial
import tempfile
import os
from ..exceptions import InvalidParsingFormat
from ..utilities.file_type import filetype
from itertools import chain
import sys


class HeaderError(Exception):
    """
    Mock exception which is raised when a header/comment line (e.g. starting with "#") is found.
    """
    pass


class Parser(metaclass=abc.ABCMeta):
    """Generic parser iterator. Base parser class."""

    def __init__(self, handle):
        self.__closed = False
        self._handle = self.__get_handle(handle)
        self.closed = False

    def __iter__(self):
        return self

    def __get_handle(self, handle, position=None):
        if not isinstance(handle, io.IOBase):
            if handle.endswith(".gz") or filetype(handle) == b"application/gzip":
                opener = gzip.open
            elif handle.endswith(".bz2") or filetype(handle) == b"application/x-bzip2":
                opener = bz2.open
            else:
                opener = partial(open, **{"buffering": 1})
            try:
                handle = opener(handle, "rt")
            except FileNotFoundError:
                raise FileNotFoundError("File not found: {0}".format(handle))
        if position is not None:
            handle.seek(position)
        return handle

    def __next__(self):
        return next(self._handle)

    def __enter__(self):
        if self.closed is True:
            raise ValueError('I/O operation on closed file.')
        return self

    def __exit__(self, *args):
        _ = args
        self._handle.close()
        self.closed = True

    def close(self):
        """
        Alias for __exit__
        """
        self.__exit__()

    @property
    def name(self):
        """
        Return the filename.
        """
        if hasattr(self._handle, "name"):
            return self._handle.name
        return None

    @property
    def closed(self):
        """
        Boolean flag. If True, the file has been closed already.
        """
        return self.__closed

    @closed.setter
    def closed(self, *args):
        """
        :param args: boolean flag

        This sets the closed flag of the file.

        """
        if not isinstance(args[0], bool):
            raise TypeError("Invalid value: {0}".format(args[0]))

        self.__closed = args[0]

    def __getstate__(self):
        try:
            position = self._handle.tell()
        except:
            position = None
        state = dict()
        state.update(self.__dict__)
        state["position"] = position
        if hasattr(self._handle, "filename"):
            _handle = self._handle.filename
            if isinstance(_handle, bytes):
                _handle = _handle.decode()
            state["_handle"] = _handle
        elif hasattr(self._handle, "name"):
            _handle = self._handle.name
            if isinstance(_handle, bytes):
                _handle = _handle.decode()
            state["_handle"] = _handle
        else:
            raise TypeError("Unknown handle: {}".format(self._handle))
        state.pop("logger", None)
        return state

    def __setstate__(self, state):
        position = state.get("position")
        del state["position"]
        self.__dict__.update(state)
        self._handle = self.__get_handle(state["_handle"], position=position)


# noinspection PyPep8
from . import GFF
# noinspection PyPep8
from . import GTF
# noinspection PyPep8
from . import bed12
# noinspection PyPep8
from . import blast_utils
from . import bam_parser


def to_gff(string, input_format=None):
    """
    Function to recognize the input file type (GFF or GTF).
    :param string:
    :type string: (str|io.TextIOWrapper|io.BytesIO|io.BufferedReader|IO)
    :rtype: (Mikado.parsers.GTF.GTF | Mikado.parsers.GFF.GFF3)
    """

    streaming = False
    if isinstance(string, (io.BytesIO, io.BufferedReader, io.TextIOWrapper)):
        fname = "-"
        streaming = True
    elif isinstance(string, (bytes, str)):
        if isinstance(string, bytes):
            string = string.decode()
        fname = string
    else:
        raise ValueError("Invalid input type: {}".format(type(string)))

    order = [GTF.GTF, GFF.GFF3, bed12.Bed12Parser, bam_parser.BamParser]

    if input_format == "bam" or fname.endswith(".bam"):
        first = bam_parser.BamParser
        input_format = "bam"
    elif input_format == "gtf" or".gtf" in fname:
        first = GTF.GTF
        input_format = "gtf"
    elif input_format == "gff3" or ".gff" in fname or ".gff3" in fname:
        first = GFF.GFF3
        input_format = "gff3"
    elif input_format == "bed12" or ".bed12" in fname or ".bed" in fname:
        first = bed12.Bed12Parser
        input_format = "bed12"
    elif streaming is False:
        # We will have to impute.
        first = GTF.GTF
        input_format = "gtf"
    else:
        raise InvalidParsingFormat(
            "I cannot infer the correct file type from streaming. Please provide the correct file type.")

    if streaming is False and os.stat(fname).st_size == 0 or streaming is True:
        return first(fname)

    found = False
    raised = dict()
    for test in chain([first], [_ for _ in order if _ != first]):
        try:
            parser = test(fname)
            for row in parser:
                if test.__annot_type__ == "bam" or row.header is False:
                    found = True
                    break
                else:
                    continue
            if found:
                break
        except InvalidParsingFormat as exc:
            raised[test.__annot_type__] = exc
            continue

    if found and test.__annot_type__ == "bam":
        return test(fname)
    elif found:
        return test(string)
    else:
        raise InvalidParsingFormat(
            "Invalid file specified: {} should have been of format {}, but it could not be verified. Error:\n{}".format(
             fname if fname != "-" else "stream", input_format, raised[input_format]
            ))
