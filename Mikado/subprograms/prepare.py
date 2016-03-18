#!/usr/bin/env python3
# coding: utf-8

"""
Subprogram that constitutes the first step of the Mikado.py pipeline.
"""

import sys
import os
import argparse
import logging
import logging.handlers
from ..utilities import path_join
from ..utilities.log_utils import formatter
from ..preparation.prepare import prepare
from ..configuration.configurator import to_json


__author__ = 'Luca Venturini'


def setup(args):
    """Method to set up the analysis using the JSON configuration
    and the command line options.

    :param args: the ArgumentParser-derived namespace.
    """

    logger = logging.getLogger("prepare")
    logger.setLevel(logging.INFO)

    if args.start_method is not None:
        args.json_conf["multiprocessing_method"] = args.start_method

    if args.output_dir is not None:
        args.json_conf["prepare"]["output_dir"] = getattr(args,
                                                          "output_dir")
    if not os.path.exists(args.json_conf["prepare"]["output_dir"]):
        try:
            os.makedirs(args.json_conf["prepare"]["output_dir"])
        except (OSError, PermissionError) as exc:
            logger.error("Failed to create the output directory!")
            logger.exception(exc)
            raise
    elif not os.path.isdir(args.json_conf["prepare"]["output_dir"]):
        logger.error(
            "The specified output directory %s exists and is not a file; aborting",
            args.json_conf["prepare"]["output_dir"])
        raise OSError("The specified output directory %s exists and is not a file; aborting" %
                      args.json_conf["prepare"]["output_dir"])

    if args.log is not None:
        args.log.close()
        args.json_conf["prepare"]["log"] = args.log.name

    if args.json_conf["prepare"]["log"]:
        handler = logging.FileHandler(
            path_join(
                args.json_conf["prepare"]["output_dir"],
                args.json_conf["prepare"]["log"]),
            "w")
    else:
        handler = logging.StreamHandler()

    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.info("Command line: %s",  " ".join(sys.argv))
    logger.propagate = False

    if args.verbose is True:
        args.json_conf["log_settings"]["log_level"] = "DEBUG"
    elif args.quiet is True:
        args.json_conf["log_settings"]["log_level"] = "WARN"

    args.level = args.json_conf["log_settings"]["log_level"]
    logger.setLevel(args.level)

    if args.gff:
        args.json_conf["prepare"]["gff"] = args.gff
    else:
        if not args.json_conf["prepare"]["gff"]:
            parser = prepare_parser()
            print(parser.format_help())
            sys.exit(0)

    if args.labels != '':
        args.labels = args.labels.split(",")
        # Checks labels are unique
        assert len(set(args.labels)) == len(args.labels)
        assert not any([True for _ in args.labels if _.strip() == ''])
        if len(args.labels) != len(args.json_conf["prepare"]["gff"]):
            raise ValueError("Incorrect number of labels specified")
        args.json_conf["prepare"]["labels"] = args.labels
    else:
        if not args.json_conf["prepare"]["labels"]:
            args.labels = [""] * len(args.json_conf["prepare"]["gff"])
            args.json_conf["prepare"]["labels"] = args.labels

    for option in ["out", "out_fasta", "fasta",
                   "minimum_length", "procs", "single"]:
        if ((getattr(args, option) or getattr(args, option) == 0) and
                getattr(args, option) is not False):
            args.json_conf["prepare"][option] = getattr(args, option)

    if args.lenient is not None:
        args.json_conf["prepare"]["lenient"] = True

    if args.strand_specific is True:
        args.json_conf["prepare"]["strand_specific"] = True

    if args.strip_cds is True:
        args.json_conf["prepare"]["strip_cds"] = True

    return args, logger


def prepare_launcher(args):

    args, logger = setup(args)
    prepare(args, logger)


def prepare_parser():
    """
    This function defines the parser for the command line interface
    of the program.
    :return: an argparse.Namespace object
    :rtype: argparse.ArgumentParser
    """

    def to_cpu_count(string):
        """
        :param string: cpu requested
        :rtype: int
        """
        try:
            string = int(string)
        except:
            raise
        return max(1, string)

    def positive(string):
        """
        Simple function to return the absolute value of the integer of the input string.
        :param string:
        :return:
        """

        return abs(int(string))

    parser = argparse.ArgumentParser("""Script to prepare a GTF for the pipeline;
    it will perform the following operations:
    1- add the "transcript" feature
    2- sort by coordinates
    3- check the strand""")
    parser.add_argument("--fasta", type=argparse.FileType(),
                        help="Genome FASTA file. Required.")
    verbosity = parser.add_mutually_exclusive_group()
    verbosity.add_argument("-v", "--verbose", action="store_true", default=False)
    verbosity.add_argument("-q", "--quiet", action="store_true", default=False)
    parser.add_argument("--start-method", dest="start_method",
                        choices=["fork", "spawn", "forkserver"],
                        default=None, help="Multiprocessing start method.")
    parser.add_argument("-s", "--strand-specific", dest="strand_specific",
                        action="store_true", default=False,
                        help="""Flag. If set, monoexonic transcripts
                        will be left on their strand rather than being
                        moved to the unknown strand.""")
    parser.add_argument("-l", "--log", type=argparse.FileType("w"), default=None,
                        help="Log file. Optional.")
    parser.add_argument("--lenient", action="store_true", default=None,
                        help="""Flag. If set, transcripts with only non-canonical
                        splices will be output as well.""")
    parser.add_argument("-m", "--minimum_length", default=200, type=positive,
                        help="Minimum length for transcripts. Default: 200 bps.")
    parser.add_argument("-p", "--procs",
                        help="Number of processors to use (default %(default)s)",
                        type=to_cpu_count, default=1)
    parser.add_argument("-scds", "--strip_cds", action="store_true", default=False,
                        help="Boolean flag. If set, ignores any CDS/UTR segment.")
    parser.add_argument("--labels", type=str, default="",
                        help="""Labels to attach to the IDs of the transcripts of the input files,
                        separated by comma.""")
    parser.add_argument("--single", action="store_true", default=False,
                        help="Disable multi-threading. Useful for debugging.")
    parser.add_argument("-od", "--output-dir", dest="output_dir",
                        type=str, default=".",
                        help="Output directory. Default: current working directory")
    parser.add_argument("-o", "--out", default=None,
                        help="Output file. Default: mikado_prepared.fasta.")
    parser.add_argument("-of", "--out_fasta", default=None,
                        help="Output file. Default: mikado_prepared.fasta.")
    parser.add_argument("--json-conf", dest="json_conf",
                        type=to_json, default="",
                        help="Configuration file.")
    parser.add_argument("gff", help="Input GFF/GTF file(s).", nargs="*")
    parser.set_defaults(func=prepare_launcher)
    return parser
