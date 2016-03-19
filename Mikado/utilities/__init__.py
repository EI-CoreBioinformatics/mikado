"""
This module contains basic utilities for the suite, like e.g. database connection
and log creation.
"""

import os
import functools
from . import dbutils
from . import log_utils
import collections
from ..parsers import to_gff
from itertools import zip_longest
from .overlap import c_overlap, c_overlap_positive

__author__ = 'Luca Venturini'


def path_join(output_dir, output_file):

    """Small utility to join together a directory path and
    an output file, checking first that the output file is not
    an absolute path.

    :param output_dir: the output directory
    :param output_file: the output file
    """

    if os.path.isabs(output_file):
        return output_file
    else:
        return os.path.join(output_dir,
                            output_file)


def memoize(obj):

    """
    Function to memorize the results of functions/properties in memory for fast access.
    Source: https://wiki.python.org/moin/PythonDecoratorLibrary#Memoize

    :param obj: the object/function to memoize
    :return:
    """

    cache = obj.cache = {}

    @functools.wraps(obj)
    def memoizer(*args, **kwargs):
        key = str(args) + str(kwargs)
        if key not in cache:
            cache[key] = obj(*args, **kwargs)
        return cache[key]
    return memoizer


def merge_partial(filenames, handle):

    """This function merges the partial files created by the multiprocessing into a single
    sorted file.

    :param filenames: the filenames to merge into the handle
    :type filenames: list[str]

    :param handle: the handle to use for printing
    :type handle: io.TextIOWrapper
    """

    current_lines = collections.defaultdict(list)

    fnames = [open(_) for _ in filenames]

    for lines in zip_longest(*fnames):
        for line in iter(_ for _ in lines if _ is not None):
            _ = line.split("/")
            index = int(_[0])
            current_lines[index].append("/".join(_[1:]))

    total = max(current_lines.keys())

    [_.close() for _ in fnames]
    [os.remove(_) for _ in filenames]

    for index in sorted(current_lines.keys()):
        for line in current_lines[index]:
            print(line, file=handle, end="")
        del current_lines[index]

    return total


def overlap(first_interval: tuple([int, int]),
            second_interval: tuple([int, int]), flank=0,
            positive=False) -> int:

    """
    :param first_interval: a tuple of integers
    :type first_interval: (int,int)

    :param second_interval: a tuple of integers
    :type second_interval: (int,int | intervaltree.Interval)

    :param flank: an optional extending parameter to check for neighbours
    :type flank: int

    :param positive: boolean flag. If set to true, the max between overlap and 0 will be returned.
    :type positive: bool

    This static method returns the overlap between two intervals.

    Values<=0 indicate no overlap.

    The optional "flank" argument (default 0) allows to expand a locus
    upstream and downstream.
    As a static method, it can be used also outside of any instance -
    "abstractlocus.overlap()" will function.
    Input: two 2-tuples of integers.
    """

    if positive is False:
        return c_overlap(first_interval[0], first_interval[1],
                         second_interval[0], second_interval[1],
                         flank)
    else:
        return c_overlap_positive(first_interval[0], first_interval[1],
                                  second_interval[0], second_interval[1],
                                  flank)


def grouper(iterable, n):
    """
    Function to chunk an iterable into slices of at most n elements.
    :param iterable:
    :param n:
    :return:
    """

    temp = []
    for val in iterable:
        temp.append(val)
        if len(temp) >= n:
            yield temp
            temp = []

    yield temp