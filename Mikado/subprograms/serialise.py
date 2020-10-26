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
from ..utilities import path_join, comma_split
from ..utilities.log_utils import create_default_logger, formatter
from ..utilities import dbutils
from ..exceptions import InvalidJson
from ..exceptions import InvalidSerialization
import random
from ..utilities import blast_keys


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

    from ..serializers import blast_serializer
    xml_serializer = blast_serializer.BlastSerializer(
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

    if not args.json_conf["serialise"]["files"]["junctions"]:
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
        logger.debug("Loading junctions: %s", junction_file)
        from ..serializers import junction
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
                     os.listdir(xml) if _xml.endswith(
                        (".xml", ".asn", ".tsv", ".txt", ".xml.gz", ".asn.gz", ".tsv.gz", ".txt.gz")) is True])
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

    if len(args.json_conf["serialise"]["files"]["orfs"]) > 0:
        logger.info("Starting to load ORF data")
        for orf_file in args.json_conf["serialise"]["files"]["orfs"]:
            logger.debug("Starting to load ORFs from %s", orf_file)
            from ..serializers import orf
            try:
                serializer = orf.OrfSerializer(orf_file,
                                               json_conf=args.json_conf,
                                               logger=logger)
                serializer()
            except InvalidSerialization:
                logger.critical("Mikado serialise failed due to problems with the input data. Please check the logs.")
                os.remove(args.json_conf["db_settings"]["db"])
                sys.exit(1)
        logger.info("Finished loading ORF data")
    else:
        logger.info("No ORF data provided, skipping")


def load_external(args, logger):

    """Function to load external data from."""

    if args.json_conf["serialise"]["files"]["external_scores"] in (None, ""):
        logger.debug("No external scores to load, returning")
        return
    else:
        logger.info("Starting to load external data")
        from ..serializers import external
        with external.ExternalSerializer(
                args.json_conf["serialise"]["files"]["external_scores"],
                json_conf=args.json_conf,
                logger=logger) as serializer:
            serializer()
        logger.info("Finished loading external data")


def setup(args):

    """
    Function to set up everything for the serialisation.
    :param args:
    :return:
    """

    import pyfaidx
    import sqlalchemy
    from ..configuration import configurator
    args.json_conf = configurator.to_json(args.json_conf)
    logger = create_default_logger("serialiser")
    # Get the log level from general settings
    if args.start_method is not None:
        args.json_conf["multiprocessing_method"] = args.start_method

    if args.log_level is None:
        args.log_level = args.json_conf["log_settings"]["log_level"]
    else:
        args.json_conf["log_settings"]["log_level"] = args.log_level

    args.json_conf["serialise"]["log_level"] = args.json_conf["log_settings"]["log_level"]

    if args.procs is not None and args.procs > 0:
        args.json_conf["threads"] = args.procs

    # Retrieve data from the argparse and put it into the configuration
    for key in args.json_conf["serialise"]:
        if key == "files":
            for file_key in args.json_conf["serialise"]["files"]:
                if getattr(args, file_key, None):
                    if file_key in ("xml", "junctions", "orfs"):
                        setattr(args, file_key, getattr(args, file_key).split(","))
                    args.json_conf["serialise"]["files"][file_key] = getattr(args, file_key)
        elif key in ("SimpleComment", "Comment"):
            # Necessary for JSON configurations
            continue
        else:
            if getattr(args, key, None) or getattr(args, key, None) == 0:
                if getattr(args, key) is False or getattr(args, key) is None:
                    continue
                else:
                    args.json_conf["serialise"][key] = getattr(args, key)

    if args.db is not None:
        args.json_conf["db_settings"]["db"] = args.db
        args.json_conf["dbtype"] = "sqlite"

    if args.seed is not None:
        args.json_conf["seed"] = args.seed
        # numpy.random.seed((args.seed) % (2 ** 32 - 1))
        random.seed((args.seed) % (2 ** 32 - 1))
    else:
        # numpy.random.seed(None)
        random.seed(None)

    if args.output_dir is not None:
        args.json_conf["serialise"]["files"]["output_dir"] = args.output_dir
        if args.json_conf["db_settings"]["dbtype"] == "sqlite":
            args.json_conf["db_settings"]["db"] = os.path.basename(
                args.json_conf["db_settings"]["db"])

    args.json_conf["serialise"]["start_adjustment"] = args.start_adjustment

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

    if args.log is not None:
        args.json_conf["serialise"]["files"]["log"] = args.log

    args.json_conf["log_settings"]["log"] = args.json_conf["serialise"]["files"]["log"][:]

    if args.json_conf["serialise"]["files"]["log"] is not None and args.json_conf["serialise"]["files"]["log"] != "":
        if args.log != args.json_conf["serialise"]["files"]["log"] and args.log is not None:
            args.json_conf["serialise"]["files"]["log"] = args.log
        if not isinstance(args.json_conf["serialise"]["files"]["log"], str):
            args.json_conf["serialise"]["files"]["log"].close()
            args.json_conf["serialise"]["files"]["log"] = args.json_conf["serialise"]["files"]["log"].name

        log = args.json_conf["serialise"]["files"]["log"]
        if os.path.dirname(log) == "":
            args.json_conf["serialise"]["files"]["log"] = \
                os.path.join(args.json_conf["serialise"]["files"]["output_dir"],
                             os.path.basename(log))
        else:
            logdir = os.path.dirname(log).rstrip(os.path.sep)
            logdir = os.path.relpath(logdir, args.json_conf["serialise"]["files"]["output_dir"])
            args.json_conf["serialise"]["files"]["log"] = \
                os.path.join(args.json_conf["serialise"]["files"]["output_dir"],
                             logdir,
                             os.path.basename(log))
        # path_join(args.json_conf["serialise"]["files"]["output_dir"], args.json_conf["serialise"]["files"]["log"])
        handlers = logger.handlers[:]
        for handler in handlers:
            # if hasattr(handler, "baseFilename"):
            logger.removeHandler(handler)

        os.makedirs(os.path.dirname(args.json_conf["serialise"]["files"]["log"]), exist_ok=True)
        open(args.json_conf["serialise"]["files"]["log"], "wt").close()
        handler = logging.FileHandler(args.json_conf["serialise"]["files"]["log"], mode="wt", delay=False)
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    logger.setLevel("INFO")
    try:
        logger.info("Command line: %s", " ".join(sys.argv))
    except FileNotFoundError:
        for handler in logger.handlers:
            print("Handler:", handler.baseFilename)
        raise

    logger.info("Random seed: %s", args.json_conf["seed"])
    logger.setLevel(args.log_level)

    if args.json_conf["serialise"]["files"]["junctions"]:
        if args.genome_fai is not None:
            args.json_conf["reference"]["genome_fai"] = args.genome_fai
        elif args.json_conf["reference"]["genome_fai"] in (None, ""):
            if args.json_conf["reference"]["genome"] not in (None, ""):
                _ = pyfaidx.Fasta(args.json_conf["reference"]["genome"])
                args.json_conf["reference"]["genome_fai"] = _.faidx.indexname
            else:
                logger.critical("Missing FAI file for junction loading!")
                sys.exit(1)

    # File with the external scores
    if args.external_scores is not None:
        args.json_conf["serialise"]["files"]["external_scores"] = args.external_scores

    if args.max_regression is not None:
        args.json_conf["serialise"]["max_regression"] = args.max_regression

    if args.codon_table is not None:
        try:
            args.codon_table = int(args.codon_table)
        except ValueError:
            pass
        args.json_conf["serialise"]["codon_table"] = args.codon_table
    else:
        assert "codon_table" in args.json_conf["serialise"]

    # Add sqlalchemy logging
    sql_logger = logging.getLogger("sqlalchemy.engine")
    if args.log_level == "DEBUG":
        level = args.json_conf["serialise"]["log_level"]
    else:
        level = args.json_conf["log_settings"]["sql_level"]

    sql_logger.setLevel(level)
    sql_logger.addHandler(logger.handlers[0])

    logger.info("Using a %s database (location: %s)",
                args.json_conf["db_settings"]["dbtype"],
                args.json_conf["db_settings"]["db"])

    logger.info("Requested %d threads, forcing single thread: %s",
                args.json_conf["threads"],
                args.json_conf["serialise"]["single_thread"])

    if args.json_conf["serialise"]["force"] is True:
        if (args.json_conf["db_settings"]["dbtype"] == "sqlite" and
                os.path.exists(args.json_conf["db_settings"]["db"])):
            logger.warn("Removing old data from %s because force option in place",
                        args.json_conf["db_settings"]["db"])
            os.remove(args.json_conf["db_settings"]["db"])

        engine = dbutils.connect(args.json_conf)
        meta = sqlalchemy.MetaData(bind=engine)
        meta.reflect(engine)
        for tab in reversed(meta.sorted_tables):
            logger.debug("Dropping %s", tab)
            tab.drop()
            if args.json_conf["db_settings"]["dbtype"] == "mysql":
                engine.execute("OPTIMIZE TABLE {}".format(tab.name))
        if args.json_conf["db_settings"]["dbtype"] == "mysql":
            engine.execute("")
        # This would fail in MySQL as it uses the OPTIMIZE TABLE syntax above
        elif args.json_conf["db_settings"]["dbtype"] != "sqlite":
            engine.execute("VACUUM")
        dbutils.DBBASE.metadata.create_all(engine)

    return args, logger, sql_logger


def serialise(args):

    """
    Wrapper around the serializers objects. It uses the configuration
    supplied through the command line to launch the necessary tools.

    :param args: namespace with the necessary information for the serialisation
    :return:
    """

    args, logger, sql_logger = setup(args)

    # logger.info("Command line: %s",  " ".join(sys.argv))
    load_orfs(args, logger)
    load_blast(args, logger)
    load_external(args, logger)
    load_junctions(args, logger)
    logger.info("Finished")
    try:
        return 0
    except KeyboardInterrupt:
        raise
    except Exception as exc:
        logger.exception(exc)
    finally:
        logging.shutdown()
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
    orfs.add_argument("--codon-table", dest="codon_table", default=None,
                      help="""Codon table to use. Default: 0 (ie Standard, NCBI #1, but only ATG is considered \
a valid start codon.""")
    orfs.add_argument("-nsa", "--no-start-adjustment", default=True,
                      action="store_false",
                      dest="start_adjustment",
                      help="Disable the start adjustment algorithm. Useful when using e.g. TransDecoder vs 5+.")

    blast = parser.add_argument_group()
    blast.add_argument("--max-target-seqs", type=int, default=None,
                       help="Maximum number of target sequences.")
    blast.add_argument("-bt", "--blast-targets", "--blast_targets", default=[], type=comma_split,
                       help="Target sequences")
    blast.add_argument("--xml", "--tsv", type=str, dest="xml", help="""BLAST file(s) to parse.
    They can be provided in three ways:
    - a comma-separated list
    - as a base folder
    - using bash-like name expansion (*,?, etc.). In this case, you have to
    enclose the filename pattern in double quotes.

    Multiple folders/file patterns can be given, separated by a comma.
    BLAST files must be either of two formats:
- BLAST XML
- BLAST tabular format, with the following **custom** fields:
    {fields}
    """.format(fields=" ".join(blast_keys)), default=[])
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

    external_args = parser.add_argument_group()
    external_args.add_argument(
        "--external-scores",
        dest="external_scores",
        help="""Tabular file containing external scores for the transcripts.
        Each column should have a distinct name, and transcripts have to be listed on the first column.""")

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
                         dest="json_conf", type=str,
                         required=False)
    generic.add_argument("-l", "--log", type=str, default=None, nargs='?', help="Optional log file. Default: stderr")
    parser.add_argument("-od", "--output-dir", dest="output_dir",
                        type=str, default=None,
                        help="Output directory. Default: current working directory")
    generic.add_argument("-lv", "--log-level", default=None,
                         choices=["DEBUG", "INFO", "WARN", "ERROR"],
                         help="Log level. Default: derived from the configuration; if absent, INFO")
    generic.add_argument("db", type=str, default=None,
                         nargs='?',
                         help="Optional output database. Default: derived from json_conf")
    generic.add_argument("--seed", type=int, default=None,
                         help="Random seed number.")
    parser.set_defaults(func=serialise)
    return parser
