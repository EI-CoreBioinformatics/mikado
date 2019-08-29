#!/usr/bin/env python3

"""
Collection of utilities that are useful for managing Mikado-related files,
including e.g. a statistical calculator for GFF/GTFs, an awk-like utility for
GTFs, and a grep-like utility for annotation files.
"""


from . import awk_gtf
from . import metrics
from . import stats
from . import trim
from . import grep
from . import convert
from . import class_codes
import argparse

__author__ = 'Luca Venturini'


def util_parser():
    """
    Function to return all available utility parsers.
    :rtype: argparse.ArgumentParser
    """

    desc = """Collection of utilities for managing GTF/GFF files."""
    parser = argparse.ArgumentParser(prog="util",
                                     description=desc)
    utils = parser.add_subparsers()

    utils.add_parser("awk_gtf",
                     description="Script to retrieve specific feature slices from a GTF file.")
    utils.choices["awk_gtf"] = awk_gtf.awk_parser()
    utils.choices["awk_gtf"].prog = "mikado util awk_gtf"
    utils.choices["awk_gtf"].description = "Script to retrieve specific feature slices from a GTF file."

    utils.add_parser("class_codes",
                     description="Script to print out the class codes.")
    utils.choices["class_codes"] = class_codes.code_parser()
    utils.choices["class_codes"].prog = "mikado util class_codes"
    utils.choices["class_codes"].description = "Script to print out the class codes."

    utils.add_parser("convert",
                     description="Script to do GTF <-> GFF3 > BED12 conversions.")
    utils.choices["convert"] = convert.convert_parser()
    utils.choices["convert"].prog = "mikado util convert"
    utils.choices["convert"].description = "Script to do GTF <-> GFF3 > BED12 conversions."

    utils.add_parser("grep",
                     description="Script to extract specific models from GFF/GTF files.")
    utils.choices["grep"] = grep.grep_parser()
    utils.choices["grep"].prog = "mikado util grep"
    utils.choices["grep"].description = "Script to extract specific models from GFF/GTF files."

    utils.add_parser("metrics",
                     description="Simple script to obtain the documentation on the transcript metrics.")
    utils.choices["metrics"] = metrics.metric_parser()
    utils.choices["metrics"].prog = "mikado util metrics"
    utils.choices["metrics"].description = "Simple script to obtain the documentation on \
the transcript metrics."

    utils.add_parser("stats", description="""GFF/GTF statistics script.
    It will compute median/average length of RNAs, exons, CDS features, etc.""")
    utils.choices["stats"] = stats.stats_parser()
    utils.choices["stats"].prog = "mikado util stats"
    utils.choices["stats"].description = "GFF/GTF statistics script. \
It will compute median/average length of RNAs, exons, CDS features, etc."

    utils.add_parser("trim",
                     description="Script to remove up to N bps from terminal exons in an annotation file.")
    utils.choices["trim"] = trim.trim_parser()
    utils.choices["trim"].prog = "mikado util trim"
    utils.choices["trim"].description = "Script to remove up to N bps from terminal\
exons in an annotation file."

    parser.add_help = True

    return parser
