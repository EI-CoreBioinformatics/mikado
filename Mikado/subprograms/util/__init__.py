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
from . import merge_blast
from . import convert
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
                     help="Script to retrieve specific feature slices from a GTF file.")
    utils.choices["awk_gtf"] = awk_gtf.awk_parser()
    utils.choices["awk_gtf"].prog = "mikado.py util awk_gtf"

    utils.add_parser("convert",
                     help="Script to do GTF <-> GFF3 > BED12 conversions.")
    utils.choices["convert"] = convert.convert_parser()
    utils.choices["convert"].prog = "mikado.py util convert"

    utils.add_parser("grep",
                     help="Script to extract specific models from GFF/GTF files.")
    utils.choices["grep"] = grep.grep_parser()
    utils.choices["grep"].prog = "mikado.py util grep"

    utils.add_parser("metrics",
                     help="Simple script to obtain the documentation on the transcript metrics.")
    utils.choices["metrics"] = metrics.metric_parser()
    utils.choices["metrics"].prog = "mikado.py util metrics"

    utils.add_parser("stats", help="""GFF/GTF statistics script.
    It will compute median/average length of RNAs, exons, CDS features, etc.""")
    utils.choices["stats"] = stats.stats_parser()
    utils.choices["stats"].prog = "mikado.py util stats"

    utils.add_parser("trim",
                     help="Script to remove up to N bps from terminal exons in an annotation file.")
    utils.choices["trim"] = trim.trim_parser()
    utils.choices["trim"].prog = "mikado.py util trim"

    utils.add_parser("merge_blast",
                     help="""Script to merge together multiple BLAST XML files.
                     It also converts them on the fly if they are in ASN/compressed XML.""")
    utils.choices["merge_blast"] = merge_blast.merger_parser()
    utils.choices["merge_blast"].prog = "mikado.py util merge_blast"

    parser.add_help = True

    return parser
