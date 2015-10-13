#!/usr/bin/env python3

"""
Main launcher of the pipeline.
"""

import argparse
import sys
import Mikado.subprograms


def main(call_args=None):

    """
    Main launcher function for the pipeline.
    """
    if call_args is None:
        call_args = sys.argv[1:]

    parser = argparse.ArgumentParser(prog="Mikado",
                                     description="""Mikado is a program to analyse RNA-Seq data
                                     and determine the best transcript for each locus in accordance
                                     to user-specified criteria.""")
    subparsers = parser.add_subparsers(
        title="Components",
        help="These are the various components of Mikado.")

    subparsers.add_parser("prepare", description="GTF preparation script",
                          help="""Mikado prepare analyses an input GTF file and
                          prepares it for the picking analysis by sorting its transcripts
                          and performing some simple consistency checks.""")
    subparsers.choices["prepare"] = Mikado.subprograms.prepare.prepare_parser()
    subparsers.choices["prepare"].prog = "Mikado prepare"

    subparsers.add_parser("serialise", description="Data serialisation script",
                          help="""Mikado serialise creates the database used
                          by the pick program. It handles Junction and ORF BED12
                          files as well as BLAST XML results.""")
    subparsers.choices["serialise"] = Mikado.subprograms.serialise.serialise_parser()
    subparsers.choices["serialise"].prog = "Mikado serialise"

    subparsers.add_parser("pick", description="Comparison script",
                          help="""Mikado pick analyses a sorted GTF/GFF files in order
                          to identify its loci and choose the best transcripts according
                          to user-specified criteria. It is dependent on files produced
                          by the "prepare" and "serialise" components.""")
    subparsers.choices["pick"] = Mikado.subprograms.pick.pick_parser()
    subparsers.choices["pick"].prog = "Mikado pick"

    subparsers.add_parser("compare", description="Comparison between reference and prediction",
                          help="""Mikado compare produces a detailed comparison of
                          reference and prediction files. It has been directly inspired
                          by Cufflinks's cuffcompare and ParsEval.""")
    subparsers.choices["compare"] = Mikado.subprograms.compare.compare_parser()
    subparsers.choices["compare"].prog = "Mikado compare"

    subparsers.add_parser("util", description="Miscellaneous utilities",
                          help="Subparser holding various utilities of the suite.")
    subparsers.choices["util"] = Mikado.subprograms.util.util_parser()
    subparsers.choices["util"].prog = "Mikado util"

    args = parser.parse_args(call_args)
    if hasattr(args, "func"):
        args.func(args)
    elif len(call_args) > 0 and call_args[0] == "util":
        Mikado.subprograms.util.util_parser().print_help()
    else:
        parser.print_help()

if __name__ == '__main__':
    main()
