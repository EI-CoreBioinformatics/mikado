# coding: utf-8

"""This module contains the ORM modules necessary to create the starting DB from the input data."""

from mikado_lib.parsers import GTF,GFF


def to_gff(string):
    """
    Function to recognize the input file type (GFF or GTF).
    :param string:
    :type string: str
    """
    f = open(string)
    if string.endswith("gtf"):
        return GTF.GTF(f)
    elif string.endswith('gff') or string.endswith('gff3'):
        return GFF.GFF3(f)
    else:
        raise ValueError('Unrecognized format')


import mikado_lib.subprograms.compare
import mikado_lib.subprograms.pick
import mikado_lib.subprograms.prepare
import mikado_lib.subprograms.serialise
import mikado_lib.subprograms.util
