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

from Bio import SeqIO

from mikado_lib.configuration import json_utils
from mikado_lib.configuration.log_utils import create_default_logger
from mikado_lib.serializers import orf, blast_serializer, junction, dbutils

__author__ = 'Luca Venturini'


def to_seqio(string):
    """
    Convert a string to a SeqIO index.

    :param string
    :type string: str
    """
    if not (os.path.exists(string) and os.path.isfile(string) and os.stat(string).st_size > 0):
        raise OSError("Invalid sequence file: {0}".format(string))
    return SeqIO.index(string, "fasta")


def xml_launcher(xml_candidate=None, args=None, logger=None):

    """
    Thin rapper around blast_utils.XmlSerializer. Its purpose is
    to launch the

    :param xml_candidate: An XML or ASN BLAST file name
    :param args: namespace with the parameters
    :return:
    """

    # If we have a value in the JSON, use it
    max_target_conf = float("Inf")
    if "blast" in args.json_conf["blast"]:
        if "max_target_seqs" in args.json_conf["blast"]:
            max_target_conf = args.json_conf["blast"]["max_target_seqs"]

    args.max_target_seqs = min(args.max_target_seqs, max_target_conf)

    xml_serializer = blast_serializer.XmlSerializer(
        xml_candidate,
        discard_definition=args.discard_definition,
        max_target_seqs=args.max_target_seqs,
        maxobjects=args.max_objects,
        target_seqs=args.target_seqs,
        query_seqs=args.transcript_fasta,
        json_conf=args.json_conf,
        logger=logger)
    xml_serializer()


def serialise(args):

    """
    Wrapper around the serializers objects. It uses the configuration
    supplied through the command line to launch the necessary tools.

    :param args: namespace with the necessary information for the serialisation
    :return:
    """

    # Provisional logging. To be improved
    logger = create_default_logger("serialiser")
    if args.log is not None:
        args.log.close()
        formatter = logger.handlers[0].formatter
        logger.removeHandler(logger.handlers[0])
        handler = logging.FileHandler(args.log.name, "w")
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    # Add sqlalchemy logging
    sql_logger = logging.getLogger("sqlalchemy.engine")
    sql_logger.setLevel(args.log_level)
    sql_logger.addHandler(logger.handlers[0])

    logger.setLevel(args.log_level)

    if args.force is True:
        logger.warn("Removing old data because force option in place")
        engine = dbutils.connect(args.json_conf)
        meta = sqlalchemy.MetaData(bind=engine)
        meta.reflect(engine)
        for tab in reversed(meta.sorted_tables):
            tab.drop()
        dbutils.DBBASE.metadata.create_all(engine)

    if args.junctions is not None:
        logger.info("Starting to load junctions")
        for junction_file in args.junctions.split(","):
            serializer = junction.JunctionSerializer(junction_file,
                                                     fai=args.genome_fai,
                                                     json_conf=args.json_conf,
                                                     maxobjects=args.max_objects,
                                                     logger=logger)
            serializer()
        logger.info("Loaded junctions")

    if args.xml is not None:
        logger.info("Starting to load BLAST data")
        filenames = []

        part_launcher = functools.partial(xml_launcher, **{"args": args,
                                                           "logger": logger})

        for xml in args.xml.split(","):
            if os.path.isdir(xml):
                filenames.extend([os.path.join(xml, _xml) for _xml in
                                  os.listdir(xml) if (_xml.endswith(".xml") or
                                                      _xml.endswith(".xml.gz") or
                                                      _xml.endswith(".asn.gz")) is True])
            else:
                filenames.extend(glob.glob(xml))

        if len(filenames) == 0:
            raise ValueError("No valid BLAST file specified!")

        part_launcher(filenames)
        logger.info("Finished to load BLAST data")

    if args.orfs is not None:
        logger.info("Starting to load ORF data")
        for orf_file in args.orfs.split(","):
            serializer = orf.OrfSerializer(orf_file,
                                           fasta_index=args.transcript_fasta,
                                           maxobjects=args.max_objects,
                                           json_conf=args.json_conf,
                                           logger=logger)
            serializer.serialize()
        logger.info("Finished loading ORF data")


def serialise_parser():
    """
    Parser function for the serialisation step.
    :return: argparse.Namespace
    """

    parser = argparse.ArgumentParser("Serialisation utility of the Mikado suite.")
    orfs = parser.add_argument_group()
    orfs.add_argument("--orfs", type=str, default=None,
                      help="ORF BED file(s), separated by commas")
    orfs.add_argument("--transcript_fasta", default=None,
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

    blast = parser.add_argument_group()
    blast.add_argument("--max_target_seqs", type=int, default=sys.maxsize,
                       help="Maximum number of target sequences.")
    blast.add_argument("--target_seqs", default=None, type=to_seqio, help="Target sequences")
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
    """)

    junctions = parser.add_argument_group()
    junctions.add_argument("--genome_fai", default=None)
    junctions.add_argument("--junctions", type=str)
    generic = parser.add_argument_group()
    generic.add_argument("-mo", "--max-objects", dest="max_objects",
                         type=int, default=10 ** 5,
                         help="Maximum number of objects to cache in memory.")
    generic.add_argument("-f", "--force", action="store_true", default=False,
                         help="""Flag. If set, an existing databse will be deleted (sqlite)
                         or dropped (MySQL/PostGreSQL) before beginning the serialisation.""")
    generic.add_argument("--json-conf", default=None,
                         dest="json_conf", type=json_utils.to_json,
                         required=True)
    generic.add_argument("-l", "--log", type=argparse.FileType("w"), default=None,
                         nargs='?',
                         help="Optional log file. Default: stderr")
    generic.add_argument("-lv", "--log_level", default="INFO",
                         choices=["DEBUG", "INFO", "WARN", "ERROR"],
                         nargs='?',
                         help="Log level. Default: %(default)s")

    generic.add_argument("db", type=str, default=None,
                         nargs='?',
                         help="Optional output database. Default: derived from json_conf")

    parser.set_defaults(func=serialise)
    return parser
