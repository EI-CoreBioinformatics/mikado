"""
This module contains basic utilities for the suite, like e.g. database connection
and log creation.
"""

import os
import functools
from . import dbutils
from . import log_utils
from ..parsers import GTF, GFF

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


def to_gff(string):
    """
    Function to recognize the input file type (GFF or GTF).
    :param string:
    :type string: str
    """
    handle = open(string)
    if string.endswith("gtf"):
        return GTF.GTF(handle)
    elif string.endswith('gff') or string.endswith('gff3'):
        return GFF.GFF3(handle)
    else:
        raise ValueError('Unrecognized format')