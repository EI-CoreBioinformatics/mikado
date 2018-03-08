#!/usr/bin/env python3
# coding: utf_8

"""
Mikado is a Python suite whose purpose is to find and resolve genic loci in a
genomic annotation. This is the library it relies onto.
"""

__title__ = "Mikado"
__author__ = 'Luca Venturini'
__license__ = 'GPL3'
__copyright__ = 'Copyright 2015-2016 Luca Venturini'
__version__ = "1.2"

__all__ = ["configuration",
           "exceptions",
           "loci",
           "parsers",
           "picking",
           "preparation",
           "scales",
           "serializers",
           "subprograms",
           "utilities"]

from . import configuration
from . import exceptions
from . import loci
from . import parsers
from . import picking
from . import preparation
from . import scales
from . import serializers
from . import subprograms
from . import utilities
import argparse
import sys
from . import subprograms
from .utilities.log_utils import create_default_logger
from multiprocessing import freeze_support
from numpy.testing import Tester
test = Tester().test
# import pkg_resources

# __spec__ = "Mikado"


# __main__ = pkg_resources.resource_filename(subprograms.__name__, "__main__.py")


def main(call_args=None):

    """
    Main launcher function for the pipeline.
    :param call_args: optional argument string to be passed to execute the commands.
    Otherwise, the string will be derived from sys.argv[1:]
    """

    if call_args is None:
        call_args = sys.argv[1:]

    parser = argparse.ArgumentParser(prog="Mikado",
                                     description="""Mikado is a program to analyse RNA-Seq data
and determine the best transcript for each locus in accordance to user-specified criteria.""")

    parser.add_argument("--version", default=False, action="store_true",
                        help="Print Mikado current version and exit.")

    subparsers = parser.add_subparsers(
        title="Components",
        help="""These are the various components of Mikado:

""")
    subparsers.add_parser("configure",
                          help="This utility guides the user through \
the process of creating a configuration file for Mikado.")
    subparsers.choices["configure"] = subprograms.configure.configure_parser()
    subparsers.choices["configure"].prog = "Mikado configure"

    subparsers.add_parser("prepare", help ="Mikado prepare analyses an input GTF file and \
prepares it for the picking analysis by sorting its transcripts \
and performing some simple consistency checks.")
    subparsers.choices["prepare"] = subprograms.prepare.prepare_parser()
    subparsers.choices["prepare"].prog = "Mikado prepare"

    subparsers.add_parser("serialise", help="Mikado serialise creates the database used \
by the pick program. It handles Junction and ORF BED12 files as well as BLAST XML results.")
    subparsers.choices["serialise"] = subprograms.serialise.serialise_parser()
    subparsers.choices["serialise"].prog = "Mikado serialise"

    subparsers.add_parser("pick", help="Mikado pick analyses a sorted GTF/GFF files in order \
to identify its loci and choose the best transcripts according to user-specified criteria. \
It is dependent on files produced by the \"prepare\" and \"serialise\" components.")
    subparsers.choices["pick"] = subprograms.pick.pick_parser()
    subparsers.choices["pick"].prog = "Mikado pick"

    subparsers.add_parser("compare", help="Mikado compare produces a detailed comparison of \
reference and prediction files. It has been directly inspired \
by Cufflinks's cuffcompare and ParsEval.")
    subparsers.choices["compare"] = subprograms.compare.compare_parser()
    subparsers.choices["compare"].prog = "Mikado compare"

    subparsers.add_parser("util", help="Miscellaneous utilities")
    subparsers.choices["util"] = subprograms.util.util_parser()
    subparsers.choices["util"].prog = "Mikado util"

    try:
        args = parser.parse_args(call_args)
        if hasattr(args, "func"):
            args.func(args)
        elif len(call_args) > 0 and call_args[0] == "util":
            subprograms.util.util_parser().print_help()
        elif args.version is True:
            print("Mikado v{}".format( __version__))
            sys.exit(0)
        else:
            parser.print_help()
    except KeyboardInterrupt:
        raise KeyboardInterrupt
    except BrokenPipeError:
        pass
    except Exception as exc:
        logger = create_default_logger("main")
        logger.error("Mikado crashed, cause:")
        logger.exception(exc)
        sys.exit(1)

if __name__ == '__main__':
    # __spec__ = "Mikado"
    freeze_support()
    main()
