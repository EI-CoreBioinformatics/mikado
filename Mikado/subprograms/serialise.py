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
import shutil
import sys
import logging
import logging.handlers
from copy import copy
from typing import Union

import tempfile

import io

from Mikado import version

from ..utilities import path_join, comma_split, create_shm_handle
from ..utilities.dbutils import DBBASE
from ..utilities.log_utils import create_default_logger, create_null_logger
from ._utils import check_log_settings_and_create_logger
from ..utilities import dbutils
from ..exceptions import InvalidConfiguration
from ..exceptions import InvalidSerialization
import random
from ..utilities import blast_keys
from ..configuration import MikadoConfiguration, DaijinConfiguration
import pysam
import sqlalchemy
from ..configuration import configurator


__author__ = 'Luca Venturini'


def setup_shm_db(configuration: MikadoConfiguration, logger=create_null_logger()) -> (MikadoConfiguration,
                                                                                      Union[None, io.FileIO]):
    """
    This method will copy the SQLite input DB into memory.
    """

    if configuration.serialise.shm is True and configuration.db_settings.dbtype == "sqlite":
        handle = create_shm_handle()
        if handle is None:
            logger.info("Mikado was asked to copy the database into /dev/shm, but it is either " 
                        "not available or not writable for this user.")
            configuration.serialise.shm = False
            return configuration, None
        else:
            logger.info(
                f"Using {handle.name} as temporary database before copying into {configuration.db_settings.db}.")
            configuration.db_settings.db = handle.name
            return configuration, handle


def xml_launcher(xml_candidate=None, configuration=None, logger=None):

    """
    Thin rapper around blast_utils.XmlSerializer. Its purpose is
    to create a standard serializer object (use it with partial
     from functools)

    :param xml_candidate: An XML or ASN BLAST file name

    :param configuration: the configuration dictionary, if available
    :type configuration: (None | MikadoConfiguration|DaijinConfiguration)

    :param logger: the logger instance.
    :type logger: (None | logging.Logger)

    :return:
    """

    from ..serializers import blast_serializer
    xml_serializer = blast_serializer.BlastSerializer(
        xml_candidate,
        configuration=configuration,
        logger=logger)
    xml_serializer()


def load_junctions(mikado_configuration: Union[MikadoConfiguration, DaijinConfiguration], logger):
    """
    Function that performs the loading of the junctions.

    :param mikado_configuration: Mikado configuration object

    :param logger: the logging instance.
    :type logger: (None | logging.Logger)

    :return:
    """

    assert isinstance(mikado_configuration, (MikadoConfiguration, DaijinConfiguration))

    if not mikado_configuration.serialise.files.junctions:
        logger.info("Skipping junction loading as no junctions have been provided.")
        return

    if not mikado_configuration.reference.genome:
        exc = InvalidConfiguration(
            "Missing the genome FAI file for serialising the junctions. \
I cannot proceed with this step!")
        logger.exception(exc)
        raise exc

    logger.info("Starting to load junctions: %s",
                mikado_configuration.serialise.files.junctions)
    for junction_file in iter(
            j_file for j_file in mikado_configuration.serialise.files.junctions
            if j_file != ''):
        logger.debug("Loading junctions: %s", junction_file)
        from ..serializers import junction
        serializer = junction.JunctionSerializer(
            junction_file,
            configuration=mikado_configuration,
            logger=logger)
        serializer()
    logger.info("Loaded junctions")


def load_blast(mikado_configuration: Union[MikadoConfiguration, DaijinConfiguration], logger):

    """
    Function to load the BLAST data into the chosen database.

    :param mikado_configuration: Mikado configuration object

    :param logger: the logging instance.
    :type logger: (None | logging.Logger)

    """
    if mikado_configuration.serialise.files.xml:
        logger.info("Starting to load BLAST data")
        filenames = []

        part_launcher = functools.partial(
            xml_launcher,
            **{"configuration": mikado_configuration, "logger": logger})

        for xml in mikado_configuration.serialise.files.xml:
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


def load_orfs(mikado_configuration: Union[DaijinConfiguration, MikadoConfiguration], logger):

    """
    Function to load the ORFs into the DB.
    :param mikado_configuration: configuration object
    :param logger:
    :return:
    """

    if len(mikado_configuration.serialise.files.orfs) > 0:
        from ..serializers import orf
        logger.info("Starting to load ORF data")
        for orf_file in mikado_configuration.serialise.files.orfs:
            logger.debug("Starting to load ORFs from %s", orf_file)
            try:
                serializer = orf.OrfSerializer(orf_file,
                                               configuration=mikado_configuration,
                                               logger=logger)
                serializer()
            except InvalidSerialization as exc:
                logger.critical("Mikado serialise failed due to problems with the input data. Please check the logs.")
                raise exc
        logger.info("Finished loading ORF data")
    else:
        logger.info("No ORF data provided, skipping")


def load_external(mikado_configuration: Union[MikadoConfiguration, DaijinConfiguration], logger):

    """Function to load external data from."""

    if mikado_configuration.serialise.files.external_scores in (None, ""):
        logger.debug("No external scores to load, returning")
        return
    else:
        logger.info("Starting to load external data")
        from ..serializers import external
        with external.ExternalSerializer(
                mikado_configuration.serialise.files.external_scores,
                configuration=mikado_configuration,
                logger=logger) as serializer:
            serializer()
        logger.info("Finished loading external data")


def _set_serialise_files_options(conf: Union[MikadoConfiguration, DaijinConfiguration],
                                 args, logger=create_null_logger()) -> Union[MikadoConfiguration, DaijinConfiguration]:
    conf.serialise.files.orfs = args.orfs.split(",") if args.orfs else conf.serialise.files.orfs
    conf.serialise.files.xml = args.xml.split(",") if args.xml else conf.serialise.files.xml
    conf.serialise.files.junctions = args.junctions.split(",") if args.junctions else conf.serialise.files.junctions
    conf.serialise.files.transcripts = args.transcripts if args.transcripts else conf.serialise.files.transcripts
    conf.serialise.files.blast_targets = args.blast_targets if args.blast_targets else \
        conf.serialise.files.blast_targets
    conf.reference.genome_fai = args.genome_fai if args.genome_fai else \
        conf.reference.genome_fai
    conf.reference.genome = args.genome if args.genome else conf.reference.genome

    conf.serialise.files.output_dir = args.output_dir if args.output_dir is not None else \
        conf.serialise.files.output_dir

    try:
        os.makedirs(conf.serialise.files.output_dir, exist_ok=True)
    except (OSError, PermissionError, FileExistsError) as exc:
        logger.error("Failed to create the output directory!")
        logger.exception(exc)
        raise exc
        
    if args.db is not None:
        conf.db_settings.db = args.db
        conf.db_settings.dbtype = "sqlite"    

    if conf.serialise.files.output_dir != "." and conf.db_settings.dbtype == "sqlite":
        conf.db_settings.db = path_join(conf.serialise.files.output_dir,
                                        os.path.basename(conf.db_settings.db))

    if conf.serialise.files.junctions:
        if conf.reference.genome_fai in (None, ""):
            if conf.reference.genome not in (None, ""):
                _ = pysam.Fastafile(conf.reference.genome)
                conf.reference.genome_fai = conf.reference.genome + ".fai"
            else:
                exc = InvalidConfiguration("Missing FAI file for junction loading!")
                logger.critical(exc)
                raise exc

    return conf


def _set_serialise_run_options(conf: Union[MikadoConfiguration, DaijinConfiguration],
                               args, logger=create_null_logger()) -> Union[MikadoConfiguration, DaijinConfiguration]:
    conf.serialise.max_regression = args.max_regression or conf.serialise.max_regression
    conf.serialise.start_adjustment = args.start_adjustment
    conf.serialise.max_target_seqs = args.max_target_seqs or conf.serialise.max_target_seqs
    conf.serialise.single_thread = args.single_thread or conf.serialise.single_thread
    conf.serialise.shm = args.use_shm if args.use_shm is not None else conf.serialise.shm
    conf.serialise.files.blast_loading_debug = True if args.blast_loading_debug else \
        conf.serialise.files.blast_loading_debug
    conf.serialise.force = args.force if args.force is not None else conf.serialise.force
    # File with the external scores
    conf.serialise.files.external_scores = args.external_scores if args.external_scores else \
        conf.serialise.files.external_scores
    conf.serialise.codon_table = args.codon_table if isinstance(args.codon_table, (int, str)) else \
        conf.serialise.codon_table
    return conf


def _execute_force_removal(configuration: Union[MikadoConfiguration, DaijinConfiguration],
                           logger=create_null_logger()):

    if configuration.db_settings.dbtype == "sqlite":
        if os.path.exists(configuration.db_settings.db):
            logger.warn("Removing old data from %s because force option in place",
                        configuration.db_settings.db)
            os.remove(configuration.db_settings.db)
            assert not os.path.exists(configuration.db_settings.db)
        else:
            pass
        return

    engine = dbutils.connect(configuration)
    meta = sqlalchemy.MetaData(bind=engine)
    meta.reflect(engine)
    for tab in reversed(meta.sorted_tables):
        logger.debug("Dropping %s", tab)
        tab.drop()
        if configuration.db_settings.dbtype == "mysql":
            engine.execute("OPTIMIZE TABLE {}".format(tab.name))
    if configuration.db_settings.dbtype != "mysql":
        engine.execute("VACUUM")
    dbutils.DBBASE.metadata.create_all(engine)
    return


def setup(args):

    """
    Function to set up everything for the serialisation.
    :param args:
    :return:
    """

    mikado_configuration = configurator.load_and_validate_config(args.configuration)
    logger = create_default_logger("serialiser")
    mikado_configuration.multiprocessing_method = args.start_method if args.start_method is not None else \
        mikado_configuration.multiprocessing_method
    mikado_configuration.threads = args.procs if (args.procs is not None and args.procs > 0) else \
        mikado_configuration.threads

    # Retrieve data from the argparse and put it into the configuration
    mikado_configuration = _set_serialise_files_options(args=args, conf=mikado_configuration, logger=logger)
    mikado_configuration = _set_serialise_run_options(mikado_configuration, args, logger=logger)
    mikado_configuration = configurator.load_and_validate_config(mikado_configuration)
    mikado_configuration, logger = check_log_settings_and_create_logger(mikado_configuration,
                                                                        args.log, args.log_level,
                                                                        section="serialise")
    logger.setLevel("INFO")
    logger.info(f"Mikado version: {version.__version__}")
    logger.info("Command line: %s", " ".join(sys.argv))
    if args.random_seed is True:
        mikado_configuration.seed = None
    elif args.seed is not None:
        mikado_configuration.seed = args.seed
    else:
        pass

    mikado_configuration.check()
    random.seed(mikado_configuration.seed)
    logger.info("Random seed: %s", mikado_configuration.seed)
    logger.setLevel(mikado_configuration.log_settings.log_level)

    # Add sqlalchemy logging
    sql_logger = logging.getLogger("sqlalchemy.engine")
    if args.log_level == "DEBUG":
        level = mikado_configuration.log_settings.log_level
    else:
        level = mikado_configuration.log_settings.sql_level

    sql_logger.setLevel(level)
    sql_logger.addHandler(logger.handlers[0])

    logger.info("Using a %s database (location: %s)", mikado_configuration.db_settings.dbtype,
                mikado_configuration.db_settings.db)

    logger.info("Requested %d threads, forcing single thread: %s", mikado_configuration.threads,
                mikado_configuration.serialise.single_thread)

    if mikado_configuration.serialise.force is True:
        _execute_force_removal(mikado_configuration, logger=logger)
        if (args.force is True or mikado_configuration.serialise.force is True) and mikado_configuration.db_settings.dbtype == "sqlite":
            assert not os.path.exists(mikado_configuration.db_settings.db)

    return mikado_configuration, logger, sql_logger


def serialise(args):

    """
    Wrapper around the serializers objects. It uses the configuration
    supplied through the command line to launch the necessary tools.

    :param args: namespace with the necessary information for the serialisation
    :return:
    """

    mikado_configuration, logger, sql_logger = setup(args)
    intended_db = copy(mikado_configuration.db_settings.db)
    pre_existing = os.path.exists(intended_db)
    if mikado_configuration.serialise.shm is True and mikado_configuration.db_settings.dbtype == "sqlite":
        mikado_configuration, temp_db = setup_shm_db(mikado_configuration, logger=logger)
        if pre_existing and mikado_configuration.serialise.force is False:
            shutil.copy(intended_db, temp_db)
    elif mikado_configuration.db_settings.dbtype == "sqlite":
        temp_db = tempfile.NamedTemporaryFile(suffix=".db")
        mikado_configuration.db_settings.db = temp_db.name
    else:
        temp_db = None

    try:
        load_orfs(mikado_configuration, logger)
        load_blast(mikado_configuration, logger)
        load_external(mikado_configuration, logger)
        load_junctions(mikado_configuration, logger)
    except (ValueError,OSError,PermissionError,IndexError,InvalidSerialization) as exc:
        logger.error("Mikado crashed due to an error. Please check the logs for hints on the cause of the error; if "
                     "it is a bug, please report it to https://github.com/EI-CoreBioinformatics/mikado/issues.")
        logger.error(exc)
        assert (temp_db is None and mikado_configuration.db_settings.dbtype != "sqlite") or (
                temp_db is not None and mikado_configuration.db_settings.dbtype == "sqlite")
        if mikado_configuration.db_settings.dbtype == "sqlite":
            temp_db.close()
            if pre_existing and os.path.exists(intended_db):
                os.remove(intended_db)
        else:
            engine = dbutils.connect(mikado_configuration)
            connection = engine.connect()
            for tbl in reversed(DBBASE.metadata.sorted_tables):
                connection.execute(tbl.delete())
        logging.shutdown()
        sys.exit(1)

    logger.info("Finished")

    if temp_db is not None:
        logger.info(f"Copying from the temporary database {mikado_configuration.db_settings.db} "
                    f"back to disk {intended_db}.")
        shutil.copy(mikado_configuration.db_settings.db, intended_db)
        logger.info(f"Copied from the temporary database {mikado_configuration.db_settings.db} "
                    f"back to disk {intended_db}.")
        temp_db.close()

    logging.shutdown()


def serialise_parser():
    """
    Parser function for the serialisation step.
    :return: argparse.Namespace
    """

    parser = argparse.ArgumentParser("Serialisation utility of the Mikado suite.")
    parser.add_argument("--start-method", dest="start_method",
                        choices=["fork", "spawn", "forkserver"],
                        default=None, help="Multiprocessing start method.")
    shm = parser.add_mutually_exclusive_group()
    shm.add_argument("--shm", action="store_true", default=None, dest="use_shm",
                     help="Use /dev/shm (if available) for faster database building.")
    shm.add_argument("--no-shm", action="store_false", default=None, dest="use_shm",
                     help="Force building the database on its final location, even if /dev/shm is available.")
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
    junctions.add_argument("--genome", default=None)
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
committing to the database. Default: 100,000 i.e. approximately 450MB RAM usage for Drosophila.""")
    forcing = parser.add_mutually_exclusive_group()
    forcing.add_argument("--no-force", action="store_false", default=None, dest="force",
                         help="""Flag. If set, do not drop the contents of an existing Mikado DB
before beginning the serialisation.""")
    forcing.add_argument("--force", action="store_true", default=None, dest="force",
                         help="""Flag. If set, an existing databse will be deleted (sqlite)
or dropped (MySQL/PostGreSQL) before beginning the serialisation.""")
    # If None, the default configuration will be used (from the blueprint)
    generic.add_argument("--configuration", "--json-conf", default=None,
                         dest="configuration", type=str,
                         required=False)
    generic.add_argument("-l", "--log", type=str, default=None, nargs='?', help="Optional log file. Default: stderr")
    parser.add_argument("-od", "--output-dir", dest="output_dir",
                        type=str, default=None,
                        help="Output directory. Default: current working directory")
    log_arguments = parser.add_mutually_exclusive_group()
    log_arguments.add_argument("-lv", "--log-level", default=None,
                         choices=["DEBUG", "INFO", "WARN", "ERROR"],
                         help="Log level. Default: derived from the configuration; if absent, INFO")
    log_arguments.add_argument("--verbose", default=None, dest="log_level", action="store_const", const="DEBUG")
    log_arguments.add_argument("--quiet", default=None, dest="log_level", action="store_const", const="WARNING")
    log_arguments.add_argument("--blast-loading-debug", default=None,
                               dest="blast_loading_debug", action="store_true",
                               help="Flag. If set, Mikado will switch on the debug mode for the XML/TSV loading.")
    generic.add_argument("db", type=str, default=None,
                         nargs='?',
                         help="Optional output database. Default: derived from configuration")
    seed_group = parser.add_mutually_exclusive_group()
    seed_group.add_argument("--seed", type=int, default=None, help="Random seed number. Default: 0.")
    seed_group.add_argument("--random-seed", action="store_true", default=False,
                            help="Generate a new random seed number (instead of the default of 0)")
    parser.set_defaults(func=serialise)
    return parser
