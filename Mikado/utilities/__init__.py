"""
This module contains basic utilities for the suite, like e.g. database connection
and log creation.
"""

import os
import functools
from . import dbutils
from . import log_utils
import collections
import gzip
import numpy
try:
    import simplejson as json
except ImportError:
    import json
from itertools import zip_longest
from .overlap import overlap
# from ..parsers import to_gff

__author__ = 'Luca Venturini'


def comma_split(string):
    """Small utility to split a string based on comma. Useful for parsers."""

    return string.split(",")


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


def merge_partial(filenames, handle, logger=None, gzipped=False):

    """This function merges the partial files created by the multiprocessing into a single
    sorted file.

    :param filenames: the filenames to merge into the handle
    :type filenames: list[str]

    :param handle: the handle to use for printing
    :type handle: io.TextIOWrapper

    :param logger: logger to be used for the merging

    """

    if logger is None:
        logger = log_utils.create_null_logger("merger")
    if len(filenames) == 0:
        logger.warning("Nothing to merge. Exiting")

    logger.debug("Starting to merge %d files (root: %s)",
                 len(filenames), "-".join(filenames[0].split("-")[:-1]))

    current_lines = collections.defaultdict(list)

    try:
        if gzipped is False:
            fnames = [open(_) for _ in filenames if os.stat(_).st_size > 0]
        else:
            fnames = [gzip.open(_, "rt") for _ in filenames if os.stat(_).st_size > 0]
    except FileNotFoundError:
        raise FileNotFoundError(os.listdir(os.path.dirname(filenames[0])))

    if len(fnames) == 0:
        logger.warning("All the files to merge (root %s) are empty. Exiting.",
                       "-".join(filenames[0].split("-")[:-1]))
        [_.close() for _ in fnames]
        # [os.remove(_) for _ in filenames]

        return 0

    for lines in zip_longest(*fnames):
        for line in iter(_ for _ in lines if _ is not None):
            _ = line.split("/")
            index = int(_[0])
            current_lines[index].append("/".join(_[1:]))

    if len(current_lines) == 0:
        logger.exception("Nothing found to merge  for root %s. ERROR!.",
                         "-".join(filenames[0].split("-")[:-1]))
        [_.close() for _ in fnames]
        [os.remove(_) for _ in filenames]

        raise IndexError

    total = max(current_lines.keys())

    [_.close() for _ in fnames]
    [os.remove(_) for _ in filenames]

    for index in sorted(current_lines.keys()):
        for line in current_lines[index]:
            print(line, file=handle, end="")
        del current_lines[index]

    return total


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
    if temp:
        yield temp


def merge_dictionaries(dict_a, dict_b, path=None):
    """Recursive function to merge two dictionaries.

    :param dict_a: first dictionary
    :type dict_a: dict

    :param dict_b: second dictionary
    :type dict_b: dict

    :param path: list to be updated during recursion to indicate
                 that we are in a secondary node
    :type path: list(str)

    Source: http://stackoverflow.com/questions/7204805/dictionaries-of-dictionaries-merge
    """

    if path is None:
        path = []
    for key in dict_b:
        if key in dict_a and isinstance(dict_a[key], dict) and isinstance(dict_b[key], dict):
            merge_dictionaries(
                dict_a[key],
                dict_b[key], path + [str(key)])
        else:
            dict_a[key] = dict_b[key]
    return dict_a


def merge_ranges(ranges):
    """
    Merge overlapping and adjacent ranges and yield the merged ranges
    in order. The argument must be an iterable of pairs (start, stop).

    >>> list(merge_ranges([(5,7), (3,5), (-1,3)]))
    [(-1, 7)]
    >>> list(merge_ranges([(5,6), (3,4), (1,2)]))
    [(1, 2), (3, 4), (5, 6)]
    >>> list(merge_ranges([]))
    []
    """
    ranges = iter(sorted(ranges))
    try:
        _ = next(ranges)
        current_start, current_stop = _
    except ValueError:
        raise ValueError((_, list(ranges)))
    except StopIteration:
        return
    for start, stop in ranges:
        if start > current_stop:
            # Gap between segments: output current segment and start a new one.
            yield current_start, current_stop
            current_start, current_stop = start, stop
        else:
            # Segments adjacent or overlapping: merge.
            current_stop = max(current_stop, stop)
    yield current_start, current_stop


class NumpyEncoder(json.JSONEncoder):
    """Necessary to avoid crashes with numpy integers"""
    def default(self, obj):
        if isinstance(obj, numpy.integer):
            return int(obj)
        elif isinstance(obj, numpy.floating):
            return float(obj)
        elif isinstance(obj, numpy.ndarray):
            return obj.tolist()
        else:
            return super(NumpyEncoder, self).default(obj)
