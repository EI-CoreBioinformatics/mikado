#!/usr/bin/env python3

"""
Collection of utilities that are useful for managing mikado_lib-related files,
including e.g. a statistical calculator for GFF/GTFs, an awk-like utility for
GTFs, and a grep-like utility for annotation files.
"""


import mikado_lib.subprograms.util.awk_gtf
import mikado_lib.subprograms.util.metrics
import mikado_lib.subprograms.util.stats
import mikado_lib.subprograms.util.trim
import mikado_lib.subprograms.util.grep
import mikado_lib.subprograms.util.merge_blast
import argparse

__author__ = 'Luca Venturini'


def util_parser():
    """
    Function to return all available utility parsers.
    :rtype: argparse.Namespace
    """

    desc = """Collection of utilities for managing GTF/GFF files."""
    parser = argparse.ArgumentParser(prog="util",
                                     description=desc)
    utils = parser.add_subparsers()

    utils.add_parser("awk_gtf",
                     help="Script to retrieve specific feature slices from a GTF file.")
    utils.choices["awk_gtf"] = mikado_lib.subprograms.util.awk_gtf.awk_parser()
    utils.choices["awk_gtf"].prog = "mikado_lib util awk_gtf"

    utils.add_parser("grep",
                     help="Script to extract specific models from GFF/GTF files.")
    utils.choices["grep"] = mikado_lib.subprograms.util.grep.grep_parser()
    utils.choices["grep"].prog = "mikado_lib util grep"

    utils.add_parser("metrics",
                     help="Simple script to obtain the documentation on the transcript metrics.")
    utils.choices["metrics"] = mikado_lib.subprograms.util.metrics.metric_parser()
    utils.choices["metrics"].prog = "mikado_lib util metrics"

    utils.add_parser("stats", help="""GFF/GTF statistics script.
    It will compute median/average length of RNAs, exons, CDS features, etc.""")
    utils.choices["stats"] = mikado_lib.subprograms.util.stats.stats_parser()
    utils.choices["stats"].prog = "mikado_lib util stats"

    utils.add_parser("trim",
                     help="Script to remove up to N bps from terminal exons in an annotation file.")
    utils.choices["trim"] = mikado_lib.subprograms.util.trim.trim_parser()
    utils.choices["trim"].prog = "mikado_lib util trim"

    utils.add_parser("merge_blast",
                     help="""Script to merge together multiple BLAST XML files.
                     It also converts them on the fly if they are in ASN/compressed XML.""")
    utils.choices["merge_blast"] = mikado_lib.subprograms.util.merge_blast.merger_parser()
    utils.choices["merge_blast"].prog = "mikado_lib util merge_blast"

    parser.add_help = True

    return parser
