#!/usr/bin/env python3

"""
This subprogram is the second step in the pipeline. Its purpose is to
create the database that is needed by pick for fast access to the data
necessary during the analysis.
"""

import argparse
import functools
import glob
import os
import sys
import logging
import logging.handlers
import sqlalchemy
from ..configuration import configurator
from ..utilities import path_join, comma_split
from ..utilities.log_utils import create_default_logger, formatter
from ..utilities import dbutils
from ..serializers import orf, blast_serializer, junction
from ..exceptions import InvalidJson
import pyfaidx

__author__ = 'Luca Venturini'


def xml_launcher(xml_candidate=None, json_conf=None, logger=None):

    """
    Thin rapper around blast_utils.XmlSerializer. Its purpose is
    to create a standard serializer object (use it with partial
     from functools)

    :param xml_candidate: An XML or ASN BLAST file name

    :param json_conf: the configuration dictionary, if available
    :type json_conf: (None | dict)

    :param logger: the logger instance.
    :type logger: (None | logging.Logger)

    :return:
    """

    xml_serializer = blast_serializer.XmlSerializer(
        xml_candidate,
        json_conf=json_conf,
        logger=logger)
    xml_serializer()


def load_junctions(args, logger):
    """
    Function that performs the loading of the junctions.

    :param args: the Namespace with all the details from the command line.

    :param logger: the logging instance.
    :type logger: (None | logging.Logger)

    :return:
    """

    if args.json_conf["serialise"]["files"]["junctions"] is None:
        logger.info("Skipping junction loading as no junctions have been provided.")
        return

    if not args.json_conf["reference"]["genome"]:
        exc = InvalidJson(
            "Missing the genome FAI file for serialising the junctions. \
I cannot proceed with this step!")
        logger.exception(exc)
        raise exc

    logger.info("Starting to load junctions: %s",
                args.json_conf["serialise"]["files"]["junctions"])
    for junction_file in iter(
            j_file for j_file in args.json_conf["serialise"]["files"]["junctions"]
            if j_file != ''):
        logger.info("Loading junctions: %s", junction_file)
        serializer = junction.JunctionSerializer(
            junction_file,
            json_conf=args.json_conf,
            logger=logger)
        serializer()
    logger.info("Loaded junctions")


def load_blast(args, logger):

    """
    Function to load the BLAST data into the chosen database.

    :param args: the Namespace with all the details from the command line.

    :param logger: the logging instance.
    :type logger: (None | logging.Logger)

    """
    if args.json_conf["serialise"]["files"]["xml"]:
        logger.info("Starting to load BLAST data")
        filenames = []

        part_launcher = functools.partial(
            xml_launcher,
            **{"json_conf": args.json_conf, "logger": logger})

        for xml in args.json_conf["serialise"]["files"]["xml"]:
            if os.path.isdir(xml):
                filenames.extend(
                    [os.path.join(xml, _xml) for _xml in
                     os.listdir(xml) if (_xml.endswith(".xml") or
                                         _xml.endswith(".xml.gz") or
                                         _xml.endswith(".asn.gz")) is True])
            else:
                filenames.extend(glob.glob(xml))

        if len(filenames) > 0:
            part_launcher(filenames)
            logger.info("Finished to load BLAST data")
        else:
            logger.warning("No valid BLAST file specified, skipping this phase")


def load_orfs(args, logger):

    """
    Function to load the ORFs into the DB.
    :param args:
    :param logger:
    :return:
    """

    if args.json_conf["serialise"]["files"]["orfs"] is not None:
        logger.info("Starting to load ORF data")
        for orf_file in args.json_conf["serialise"]["files"]["orfs"]:
            serializer = orf.OrfSerializer(orf_file,
                                           json_conf=args.json_conf,
                                           logger=logger)
            serializer.serialize()
        logger.info("Finished loading ORF data")


def setup(args):

    """
    Function to set up everything for the serialisation.
    :param args:
    :return:
    """

    logger = create_default_logger("serialiser")
    # Get the log level from general settings
    if args.start_method is not None:
        args.json_conf["multiprocessing_method"] = args.start_method

    if args.log_level is None:
        args.log_level = args.json_conf["log_settings"]["log_level"]
    else:
        args.json_conf["log_settings"]["log_level"] = args.log_level

    args.json_conf["serialise"]["log_level"] = args.json_conf["log_settings"]["log_level"]

    # Retrieve data from the argparse and put it into the configuration
    for key in args.json_conf["serialise"]:
        if key == "files":
            for file_key in args.json_conf["serialise"]["files"]:
                if getattr(args, file_key):
                    if file_key in ("xml", "junctions", "orfs"):
                        setattr(args, file_key, getattr(args, file_key).split(","))
                    args.json_conf["serialise"]["files"][file_key] = getattr(args, file_key)
        elif key in ("SimpleComment", "Comment"):
            # Necesarry for JSON configurations
            continue
        else:
            if getattr(args, key) or getattr(args, key) == 0:
                if getattr(args, key) is False or getattr(args, key) is None:
                    continue
                else:
                    args.json_conf["serialise"][key] = getattr(args, key)

    if args.db is not None:
        args.json_conf["db_settings"]["db"] = args.db
        args.json_conf["dbtype"] = "sqlite"

    if args.output_dir is not None:
        args.json_conf["serialise"]["files"]["output_dir"] = args.output_dir
        if args.json_conf["db_settings"]["dbtype"] == "sqlite":
            args.json_conf["db_settings"]["db"] = os.path.basename(
                args.json_conf["db_settings"]["db"])

    if not os.path.exists(args.json_conf["serialise"]["files"]["output_dir"]):
        try:
            os.makedirs(args.json_conf["serialise"]["files"]["output_dir"])
        except (OSError, PermissionError) as exc:
            logger.error("Failed to create the output directory!")
            logger.exception(exc)
            raise
    elif not os.path.isdir(args.json_conf["serialise"]["files"]["output_dir"]):
        logger.error(
            "The specified output directory %s exists and is not a file; aborting",
            args.json_conf["serialise"]["files"]["output_dir"])
        raise OSError("The specified output directory %s exists and is not a file; aborting" %
                      args.json_conf["prepare"]["files"]["output_dir"])

    if args.json_conf["db_settings"]["dbtype"] == "sqlite":
        args.json_conf["db_settings"]["db"] = path_join(
            args.json_conf["serialise"]["files"]["output_dir"],
            args.json_conf["db_settings"]["db"])

    if args.json_conf["serialise"]["files"]["log"] is not None:
        if not isinstance(args.json_conf["serialise"]["files"]["log"], str):
            args.json_conf["serialise"]["files"]["log"].close()
            args.json_conf["serialise"]["files"]["log"] = args.json_conf[
                "serialise"]["files"]["log"].name
        args.json_conf["serialise"]["files"]["log"] = \
            path_join(
                args.json_conf["serialise"]["files"]["output_dir"],
                args.json_conf["serialise"]["files"]["log"])
        logger.removeHandler(logger.handlers[0])
        handler = logging.FileHandler(
            args.json_conf["serialise"]["files"]["log"], "w")
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    if args.json_conf["serialise"]["files"]["junctions"] is not None:
        if args.genome_fai is not None:
            args.json_conf["reference"]["genome_fai"] = args.genome_fai
        elif args.json_conf["reference"]["genome_fai"] is None:
            if args.json_conf["reference"]["genome"] is not None:
                _ = pyfaidx.Fasta(args.json_conf["reference"]["genome"])
                args.json_conf["reference"]["genome_fai"] = _.faidx.indexname
            else:
                logger.critical("Missing FAI file for junction loading!")
                sys.exit(1)

    if args.max_regression is not None:
        args.json_conf["serialise"]["max_regression"] = args.max_regression

    logger.setLevel(args.log_level)
    logger.info("Command line: %s",
                " ".join(sys.argv))

    # Add sqlalchemy logging
    sql_logger = logging.getLogger("sqlalchemy.engine")
    if args.log_level == "DEBUG":
        level = args.json_conf["serialise"]["log_level"]
    else:
        level = args.json_conf["log_settings"]["sql_level"]

    sql_logger.setLevel(level)
    sql_logger.addHandler(logger.handlers[0])

    logger.info("Requested %d threads, forcing single thread: %s",
                args.json_conf["serialise"]["procs"],
                args.json_conf["serialise"]["single_thread"])

    return args, logger


def serialise(args):

    """
    Wrapper around the serializers objects. It uses the configuration
    supplied through the command line to launch the necessary tools.

    :param args: namespace with the necessary information for the serialisation
    :return:
    """

    args, logger = setup(args)

    if args.json_conf["serialise"]["force"] is True:
        logger.warn("Removing old data because force option in place")
        engine = dbutils.connect(args.json_conf)
        meta = sqlalchemy.MetaData(bind=engine)
        meta.reflect(engine)
        for tab in reversed(meta.sorted_tables):
            logger.warn("Dropping %s", tab)
            tab.drop()
        dbutils.DBBASE.metadata.create_all(engine)

    load_junctions(args, logger)
    load_orfs(args, logger)
    load_blast(args, logger)
    logger.info("Finished")
    try:
        return 0
    except KeyboardInterrupt:
        raise
    except Exception as exc:
        logger.exception(exc)
    finally:
        return 0


def serialise_parser():
    """
    Parser function for the serialisation step.
    :return: argparse.Namespace
    """

    parser = argparse.ArgumentParser("Serialisation utility of the Mikado suite.")
    parser.add_argument("--start-method", dest="start_method",
                        choices=["fork", "spawn", "forkserver"],
                        default=None, help="Multiprocessing start method.")
    orfs = parser.add_argument_group()
    orfs.add_argument("--orfs", type=str, default=None,
                      help="ORF BED file(s), separated by commas")
    orfs.add_argument("--transcripts", default=None,
                      help="""Transcript FASTA file(s) used for ORF calling and BLAST queries,
                      separated by commas.
                      If multiple files are given, they must be in the same order of the
                      ORF files.
                      E.g. valid command lines are:

                      --transcript_fasta all_seqs1.fasta --orfs all_orfs.bed
                      --transcript_fasta seq1.fasta,seq2.fasta --orfs orfs1.bed,orf2.bed
                      --transcript_fasta all_seqs.fasta --orfs orfs1.bed,orf2.bed

                      These are invalid instead:

                      # Inverted order
                      --transcript_fasta seq1.fasta,seq2.fasta --orfs orfs2.bed,orf1.bed
                      #Two transcript files, one ORF file
                      --transcript_fasta seq1.fasta,seq2.fasta --orfs all_orfs.bed
                      """)
    orfs.add_argument("-mr", "--max-regression", dest="max_regression",
                      type=float, default=None,
                      help=""""Amount of sequence in the ORF (in %%) to backtrack
                      in order to find a valid START codon, if one is absent. Default: %(default)s""")

    blast = parser.add_argument_group()
    blast.add_argument("--max_target_seqs", type=int, default=None,
                       help="Maximum number of target sequences.")
    blast.add_argument("--blast_targets", default=[], type=comma_split,
                       help="Target sequences")
    blast.add_argument("--discard-definition", action="store_true", default=False,
                       help="""Flag. If set, the sequences IDs instead of their definition
                       will be used for serialisation.""")
    blast.add_argument("--xml", type=str, help="""XML file(s) to parse.
    They can be provided in three ways:
    - a comma-separated list
    - as a base folder
    - using bash-like name expansion (*,?, etc.). In this case, you have to
    enclose the filename pattern in double quotes.

    Multiple folders/file patterns can be given, separated by a comma.
    """, default=[])
    blast.add_argument("-p", "--procs", type=int,
                       help="""Number of threads to use for
    analysing the BLAST files. This number should not be higher than the total number of XML files.
    """, default=None)
    blast.add_argument("--single-thread", action="store_true",
                       default=None, dest="single_thread",
                       help="""Force serialise to run with a single thread, irrespective of
                       other configuration options.""")

    junctions = parser.add_argument_group()
    junctions.add_argument("--genome_fai", default=None)
    junctions.add_argument("--junctions", type=str, default=None)
    generic = parser.add_argument_group()
    generic.add_argument("-mo", "--max-objects", dest="max_objects",
                         type=int, default=None,  # So it can actually be set through the JSON
                         help="""Maximum number of objects to cache in memory before
                         committing to the database. Default: 100,000 i.e.
                         approximately 450MB RAM usage for Drosophila.""")
    generic.add_argument("-f", "--force", action="store_true", default=False,
                         help="""Flag. If set, an existing databse will be deleted (sqlite)
                         or dropped (MySQL/PostGreSQL) before beginning the serialisation.""")
    # If None, the default configuration will be used (from the blueprint)
    generic.add_argument("--json-conf", default=None,
                         dest="json_conf", type=configurator.to_json,
                         required=True)
    generic.add_argument("-l", "--log", type=str, default=None,
                         nargs='?',
                         help="Optional log file. Default: stderr")
    parser.add_argument("-od", "--output-dir", dest="output_dir",
                        type=str, default=None,
                        help="Output directory. Default: current working directory")
    generic.add_argument("-lv", "--log_level", default=None,
                         choices=["DEBUG", "INFO", "WARN", "ERROR"],
                         help="Log level. Default: derived from the configuration; if absent, INFO")

    generic.add_argument("db", type=str, default=None,
                         nargs='?',
                         help="Optional output database. Default: derived from json_conf")

    parser.set_defaults(func=serialise)
    return parser
