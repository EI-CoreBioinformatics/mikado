#!/usr/bin/env python3

"""
Main launcher of the pipeline.
"""

import argparse
import sys
import mikado_lib.subprograms


def main(call_args=None):

    """
    Main launcher function for the pipeline.
    """
    if call_args is None:
        call_args = sys.argv[1:]

    parser = argparse.ArgumentParser(prog="mikado_lib",
                                     description="""mikado_lib is a program to analyse RNA-Seq data
                                     and determine the best transcript for each locus in accordance
                                     to user-specified criteria.""")
    subparsers = parser.add_subparsers(
        title="Components",
        help="These are the various components of mikado_lib.")
    subparsers.add_parser("configure", description="Command to create the necessary configuration files.",
                          help="""This utility guides the user throught the process of creating
                           a configuration file for mikado_lib.""")
    subparsers.choices["configure"] = mikado_lib.subprograms.configure.configure_parser()
    subparsers.choices["configure"].prog = "mikado_lib configure"


    subparsers.add_parser("prepare", description="GTF preparation script",
                          help="""mikado_lib prepare analyses an input GTF file and
                          prepares it for the picking analysis by sorting its transcripts
                          and performing some simple consistency checks.""")
    subparsers.choices["prepare"] = mikado_lib.subprograms.prepare.prepare_parser()
    subparsers.choices["prepare"].prog = "mikado_lib prepare"

    subparsers.add_parser("serialise", description="Data serialisation script",
                          help="""mikado_lib serialise creates the database used
                          by the pick program. It handles Junction and ORF BED12
                          files as well as BLAST XML results.""")
    subparsers.choices["serialise"] = mikado_lib.subprograms.serialise.serialise_parser()
    subparsers.choices["serialise"].prog = "mikado_lib serialise"

    subparsers.add_parser("pick", description="Comparison script",
                          help="""mikado_lib pick analyses a sorted GTF/GFF files in order
                          to identify its loci and choose the best transcripts according
                          to user-specified criteria. It is dependent on files produced
                          by the "prepare" and "serialise" components.""")
    subparsers.choices["pick"] = mikado_lib.subprograms.pick.pick_parser()
    subparsers.choices["pick"].prog = "mikado_lib pick"

    subparsers.add_parser("compare", description="Comparison between reference and prediction",
                          help="""mikado_lib compare produces a detailed comparison of
                          reference and prediction files. It has been directly inspired
                          by Cufflinks's cuffcompare and ParsEval.""")
    subparsers.choices["compare"] = mikado_lib.subprograms.compare.compare_parser()
    subparsers.choices["compare"].prog = "mikado_lib compare"

    subparsers.add_parser("util", description="Miscellaneous utilities",
                          help="Subparser holding various utilities of the suite.")
    subparsers.choices["util"] = mikado_lib.subprograms.util.util_parser()
    subparsers.choices["util"].prog = "mikado_lib util"

    args = parser.parse_args(call_args)
    if hasattr(args, "func"):
        args.func(args)
    elif len(call_args) > 0 and call_args[0] == "util":
        mikado_lib.subprograms.util.util_parser().print_help()
    else:
        parser.print_help()

if __name__ == '__main__':
    main()
