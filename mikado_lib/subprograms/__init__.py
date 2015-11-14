# coding: utf-8

"""This module contains the ORM modules necessary to create the starting DB from the input data."""

from ..parsers import GTF, GFF


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


# noinspection PyPep8
from . import configure
# noinspection PyPep8
from . import compare
# noinspection PyPep8
from . import pick
# noinspection PyPep8
from . import prepare
# noinspection PyPep8
from . import serialise
# noinspection PyPep8
from . import util
