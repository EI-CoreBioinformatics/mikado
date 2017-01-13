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
import magic


class HeaderError(Exception):
    """
    Mock exception which is raised when a header/comment line (e.g. starting with "#") is found.
    """
    pass


class Parser(metaclass=abc.ABCMeta):
    """Generic parser iterator. Base parser class."""

    wizard = magic.Magic(mime=True)

    def __init__(self, handle):
        self.__closed = False
        if not isinstance(handle, io.IOBase):
            if handle.endswith(".gz") or self.wizard.from_file(handle) == b"application/gzip":
                opener = gzip.open
            elif handle.endswith(".bz2") or self.wizard.from_file(handle) == b"application/x-bzip2":
                opener = bz2.open
            else:
                opener = partial(open, **{"buffering": 1})
            try:
                handle = opener(handle, "rt")
            except FileNotFoundError:
                raise FileNotFoundError("File not found: {0}".format(handle))

        self._handle = handle
        self.closed = False

    def __iter__(self):
        return self

    def __next__(self):
        line = self._handle.readline()
        return line

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
        return self._handle.name

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


# noinspection PyPep8
from . import GFF
# noinspection PyPep8
from . import GTF
# noinspection PyPep8
from . import bed12
# noinspection PyPep8
from . import blast_utils


def to_gff(string):
    """
    Function to recognize the input file type (GFF or GTF).
    :param string:
    :type string: str
    :rtype: (Mikado.parsers.GTF.GTF | Mikado.parsers.GFF.GFF3)
    """

    # handle = open(string)
    if ".gtf" in string:
        return GTF.GTF(string)
    elif ".gff" in string or ".gff3" in string:
    # elif string.endswith('gff') or string.endswith('gff3'):
        return GFF.GFF3(string)
    else:
        raise ValueError('Unrecognized format for {}'.format(string))
