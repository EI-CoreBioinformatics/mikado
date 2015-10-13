#!/usr/bin/env python3
# coding: utf_8

"""
    This module defines the iterators that will parse BED12, GTF, GFF files.
"""

import io
import abc


class HeaderError(Exception):
    """
    Mock exception which is raised when a header/comment line (e.g. starting with "#") is found.
    """
    pass


class SizeError(Exception):
    """
    Custom exception
    """
    def __init__(self, value=None):
        self.value = value
        Exception.__init__(self)

    def __str__(self):
        return str(self.value)


class Parser(metaclass=abc.ABCMeta):
    """Generic parser iterator. Base parser class."""

    def __init__(self, handle):
        self.__closed = False
        if not isinstance(handle, io.IOBase):
            try:
                handle = open(handle, "rt", buffering=1)
            except:
                raise TypeError

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
import mikado_lib.parsers.GFF
# # noinspection PyPep8
import mikado_lib.parsers.GTF
# # noinspection PyPep8
import mikado_lib.parsers.bed12
# # noinspection PyPep8
# import mikado_lib.parsers.blast_utils