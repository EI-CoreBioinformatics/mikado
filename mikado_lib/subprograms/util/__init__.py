__author__ = 'luca'

import mikado_lib.subprograms.util.awk_gtf
import mikado_lib.subprograms.util.metrics
import mikado_lib.subprograms.util.stats
import mikado_lib.subprograms.util.trim
import mikado_lib.subprograms.util.grep
import argparse


def default(args):

    util_parser().print_help()


def util_parser():
    """
    Function to return all available utility parsers.
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(prog="util",
                                     description="""Collection of utilities for managing GTF/GFF files.""")
    parser.set_defaults(func=default)
    utils = parser.add_subparsers()

    utils.add_parser("awk_gtf", help="Script to retrieve specific feature slices from a GTF file.")
    utils.choices["awk_gtf"] = mikado_lib.subprograms.util.awk_gtf.awk_parser()
    utils.choices["awk_gtf"].prog = "mikado util awk_gtf"

    utils.add_parser("grep", help="Script to extract specific features from GFF/GTF files, using a list of IDs.")
    utils.choices["grep"] = mikado_lib.subprograms.util.grep.grep_parser()
    utils.choices["grep"].prog = "mikado util grep"

    utils.add_parser("metrics", help="Simple script to obtain the documentation on the metrics.")
    utils.choices["metrics"] = mikado_lib.subprograms.util.metrics.metric_parser()
    utils.choices["metrics"].prog = "mikado util metrics"

    utils.add_parser("stats", help="""GFF statistics script. It will compute median/average length of RNAs, exons,"
                                   CDS features, etc.""")
    utils.choices["stats"] = mikado_lib.subprograms.util.stats.stats_parser()
    utils.choices["stats"].prog = "mikado util stats"

    utils.add_parser("trim", help="Script to remove up to N bps from terminal exons in an annotation file.")
    utils.choices["trim"] = mikado_lib.subprograms.util.trim.trim_parser()
    utils.choices["trim"].prog = "mikado util trim"

    return parser
