#!/usr/bin/env python3
# coding: utf_8

"""
    This module defines the iterators that will parse BED12, GTF, GFF files.
"""

import io
import os
from ..exceptions import InvalidParsingFormat
from itertools import chain
from .parser import Parser
from . import GFF
from . import GTF
from . import bed12
from . import blast_utils
from . import bam_parser


def parser_factory(string, input_format=None):
    """
    Function to recognize the input file type (GFF, GTF, BED12, BAM).
    :param string:
    :type string: (str|io.TextIOWrapper|io.BytesIO|io.BufferedReader|IO)
    :rtype: (GTF.GTF | GFF.GFF3 | bam_parser.BamParser | bed12.Bed12Parser)
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
            else:
                raised[test.__annot_type__] = "No valid line found."
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
