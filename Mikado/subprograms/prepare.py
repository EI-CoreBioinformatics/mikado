#!/usr/bin/env python3
# coding: utf-8

"""
Subprogram that constitutes the first step of the Mikado pipeline.
"""

import sys
import os
import argparse
import logging
import logging.handlers
from ..configuration import MikadoConfiguration, DaijinConfiguration
from ..utilities import path_join
from ..utilities.log_utils import formatter
from ..exceptions import InvalidJson
from collections import Counter
from typing import Union


__author__ = 'Luca Venturini'


def parse_list_file(cfg, list_file):
    configuration = {
        "pick": {
            "chimera_split": {
                "skip": []
            }
        },
        "prepare": {
            "files": {
                "gff": [],
                "labels": [],
                "strand_specific_assemblies": [],
                "source_score": {},
                "reference": [],
                "exclude_redundant": [],
                "strip_cds": [],
                "skip": []
            }
        }
    }

    files_counter = Counter()

    if isinstance(list_file, str):
        list_file = open(list_file)

    for line in list_file:
        fields = line.rstrip().split("\t")
        gff_name, label, stranded = fields[:3]
        if not os.path.exists(gff_name):
            raise ValueError("Invalid file name: {}".format(gff_name))
        if label in configuration["prepare"]["files"]["labels"]:
            raise ValueError("Non-unique label specified: {}".format(label))
        if stranded.lower() not in ("true", "false"):
            raise ValueError("Malformed line for the list: {}".format(line))
        if gff_name in configuration["prepare"]["files"]["gff"]:
            raise ValueError("Repeated prediction file: {}".format(line))
        elif label != '' and label in configuration["prepare"]["files"]["labels"]:
            raise ValueError("Repeated label: {}".format(line))
        configuration["prepare"]["files"]["gff"].append(gff_name)
        configuration["prepare"]["files"]["labels"].append(label)
        if stranded.capitalize() == "True":
            configuration["prepare"]["files"]["strand_specific_assemblies"].append(gff_name)
        if len(fields) >= 4:
            try:
                score = float(fields[3])
            except ValueError:
                score = 0
            configuration["prepare"]["files"]["source_score"][label] = score
        for arr, pos, default in [("reference", 4, False), ("exclude_redundant", 5, False),
                                  ("strip_cds", 6, False), ("skip_split", 7, False)]:
            try:
                val = fields[pos]
                if val.lower() in ("false", "true"):
                    val = eval(val.capitalize())
                else:
                    raise ValueError("Malformed line. The last two fields should be either True or False.")
            except IndexError:
                val = default
            if arr == "skip_split":
                configuration["pick"]["chimera_split"]["skip"].append(val)
            else:
                configuration["prepare"]["files"][arr].append(val)

    files_counter.update(configuration["prepare"]["files"]["gff"])
    if files_counter.most_common()[0][1] > 1:
        raise InvalidJson(
            "Repeated elements among the input GFFs! Duplicated files: {}".format(
                ", ".join(_[0] for _ in files_counter.most_common() if _[1] > 1)))

    assert "exclude_redundant" in configuration["prepare"]["files"]

    cfg.prepare.files.gff = configuration["prepare"]["files"]["gff"]
    cfg.prepare.files.labels = configuration["prepare"]["files"]["labels"]
    cfg.prepare.files.strand_specific_assemblies = configuration["prepare"]["files"]["strand_specific_assemblies"]
    cfg.prepare.files.source_score = configuration["prepare"]["files"]["source_score"]
    cfg.prepare.files.reference = configuration["prepare"]["files"]["reference"]
    cfg.prepare.files.exclude_redundant = configuration["prepare"]["files"]["exclude_redundant"]
    cfg.prepare.files.strip_cds = configuration["prepare"]["files"]["strip_cds"]
    cfg.pick.chimera_split.skip = configuration["pick"]["chimera_split"]["skip"]

    return cfg


def parse_prepare_options(args, mikado_config) -> Union[DaijinConfiguration, MikadoConfiguration]:
    if args.codon_table not in (None, False, True):
        mikado_config.serialise.codon_table = str(args.codon_table)

    assert isinstance(mikado_config.reference.genome, str)
    if getattr(args, "minimum_cdna_length", None) not in (None, False):
        mikado_config.prepare.minimum_cdna_length = args.minimum_cdna_length
    if getattr(args, "max_intron_length", None) not in (None, False):
        mikado_config.prepare.max_intron_length = args.max_intron_length
    if getattr(args, "single", None) not in (None, False):
        mikado_config.prepare.single = args.single

    if args.reference is not None:
        if hasattr(args.reference, "close") and hasattr(args.reference, "name"):
            args.reference.close()
            mikado_config.reference.genome = args.reference.name
        elif hasattr(args.reference, "close") and hasattr(args.reference, "filename"):
            args.reference.close()
            mikado_config.reference.genome = args.reference.filename.decode()
        elif isinstance(args.reference, bytes):
            mikado_config.reference.genome = args.reference.decode()
        elif isinstance(args.reference, str):
            mikado_config.reference.genome = args.reference
        
    assert isinstance(mikado_config.reference.genome, str)
    if getattr(args, "exclude_redundant", None) is not None:
        mikado_config.prepare.exclude_redundant = args.exclude_redundant
    if getattr(args, "lenient", None) is not None:
        mikado_config.prepare.lenient = True
    if getattr(args, "strip_faulty_cds", None) is not None:
        mikado_config.prepare.strip_faulty_cds = True
    if getattr(args, "strip_cds", False) is True:
        mikado_config.prepare.strip_cds = True
    if args.list:
        mikado_config = parse_list_file(mikado_config, args.list)
    elif args.gff and args.gff != [""] and args.gff != []:
        __gff_counter = Counter()
        __gff_counter.update(args.gff)
        if __gff_counter.most_common()[0][1] > 1:
            raise InvalidJson(
                "Repeated elements among the input GFFs! Duplicated files: {}".format(
                    ", ".join(_[0] for _ in __gff_counter.most_common() if _[1] > 1)
                ))
        mikado_config.prepare.files.gff = args.gff
        num_files = len(mikado_config.prepare.files.gff)
        if args.strand_specific:
            mikado_config.prepare.strand_specific = True
        elif args.strand_specific_assemblies:
            args.strand_specific_assemblies = args.strand_specific_assemblies.split(",")
            if len(args.strand_specific_assemblies) > num_files:
                raise ValueError("Incorrect number of strand-specific assemblies specified!")
            for member in args.strand_specific_assemblies:
                if member not in mikado_config.prepare.files.gff:
                    raise ValueError("Incorrect assembly file specified as strand-specific")
            mikado_config.prepare.files.strand_specific_assemblies = args.strand_specific_assemblies
        if args.labels:
            args.labels = args.labels.split(",")
            # Checks labels are unique
            assert len(set(args.labels)) == len(args.labels)
            assert not any([True for _ in args.labels if _.strip() == ''])
            if len(args.labels) != num_files:
                raise ValueError("Incorrect number of labels specified")
            mikado_config.prepare.files.labels = args.labels
        else:
            if not mikado_config.prepare.files.labels:
                args.labels = list(range(1, 1 + num_files))
                mikado_config.prepare.files.labels = args.labels
        mikado_config.prepare.files.exclude_redundant = [False] * len(mikado_config.prepare.files.gff)
        mikado_config.prepare.files.reference = [False] * len(mikado_config.prepare.files.gff)

    if not mikado_config.prepare.files.exclude_redundant:
        mikado_config.prepare.files.exclude_redundant = [False] * len(mikado_config.prepare.files.gff)
    elif len(mikado_config.prepare.files.exclude_redundant) != len(mikado_config.prepare.files.gff):
        raise ValueError("Mismatch between exclude_redundant and gff files")
    if not mikado_config.prepare.files.reference:
        mikado_config.prepare.files.reference = [False] * len(mikado_config.prepare.files.gff)
    elif len(mikado_config.prepare.files.reference) != len(mikado_config.prepare.files.gff):
        raise ValueError("Mismatch between is_reference and gff files")
    if args.minimum_cdna_length:
        mikado_config.prepare.minimum_cdna_length = args.minimum_cdna_length
    if args.max_intron_length:
        mikado_config.prepare.max_intron_length = args.max_intron_length
    if getattr(args, "single", None) not in (None, False):
        mikado_config.prepare.single = args.single

    assert isinstance(mikado_config.reference.genome, str)
    return mikado_config


def setup(args, logger=None):
    """Method to set up the analysis using the JSON configuration
    and the command line options.

    :param args: the ArgumentParser-derived namespace.
    """

    if logger is None or not isinstance(logger, logging.Logger):
        logger = logging.getLogger("prepare")
        logger.setLevel(logging.INFO)

    logger.debug("Starting to get prepare arguments")
    from ..configuration.configurator import load_and_validate_config
    mikado_config = load_and_validate_config(args.configuration)
    assert hasattr(mikado_config.reference, "genome"), mikado_config.reference

    if args.start_method:
        mikado_config.multiprocessing_method = args.start_method

    if args.output_dir is not None:
        mikado_config.prepare.files.output_dir = getattr(args, "output_dir")

    if not os.path.exists(mikado_config.prepare.files.output_dir):
        try:
            os.makedirs(mikado_config.prepare.files.output_dir)
        except (OSError, PermissionError) as exc:
            logger.error("Failed to create the output directory!")
            logger.exception(exc)
            raise
    elif not os.path.isdir(mikado_config.prepare.files.output_dir):
        logger.error(
            "The specified output directory %s exists and is not a folder; aborting",
            mikado_config.prepare.files.output_dir)
        raise OSError("The specified output directory %s exists and is not a folder; aborting" %
                      mikado_config.prepare.files.output_dir)

    parse_prepare_options(args, mikado_config)
    assert hasattr(mikado_config.reference, "genome"), mikado_config.reference

    if len(mikado_config.prepare.files.gff) == 0:
        parser = prepare_parser()
        logger.error("No input files found!")
        print(parser.format_help())
        sys.exit(0)

    if args.procs is not None and args.procs > 0:
        mikado_config.threads = args.procs

    if mikado_config.prepare.files.log:
        try:
            _ = open(path_join(
                    mikado_config.prepare.files.output_dir,
                    os.path.basename(mikado_config.prepare.files.log)), "wt")
        except TypeError:
            raise TypeError((mikado_config.prepare.files.output_dir, mikado_config.prepare.files.log))

        handler = logging.FileHandler(
            path_join(
                mikado_config.prepare.files.output_dir,
                os.path.basename(mikado_config.prepare.files.log)),
            mode="wt")
    else:
        handler = logging.StreamHandler()

    handler.setFormatter(formatter)
    while logger.handlers:
        logger.removeHandler(logger.handlers.pop())
    logger.addHandler(handler)
    assert logger.handlers == [handler]
    logger.propagate = False
    logger.info("Command line: %s",  " ".join(sys.argv))
    logger.info("Random seed: %s", mikado_config.seed)

    if args.seed is not None:
        mikado_config.seed = args.seed

    if args.verbose is True:
        mikado_config.log_settings.log_level = "DEBUG"
    elif args.quiet is True:
        mikado_config.log_settings.log_level = "WARN"

    logger.setLevel(mikado_config.log_settings.log_level)

    mikado_config.prepare.files.out = os.path.basename(mikado_config.prepare.files.out)
    if getattr(args, "out") not in (None, False):
        mikado_config.prepare.files.out = os.path.basename(args.out)

    mikado_config.prepare.files.out_fasta = os.path.basename(mikado_config.prepare.files.out_fasta)
    if getattr(args, "out_fasta") not in (None, False):
        mikado_config.prepare.files.out_fasta = os.path.basename(args.out_fasta)

    if isinstance(mikado_config.reference.genome, bytes):
        mikado_config.reference.genome = mikado_config.reference.genome.decode()

    return args, mikado_config, logger


def prepare_launcher(args):

    from ..preparation.prepare import prepare
    args, mikado_config, logger = setup(args)
    assert isinstance(mikado_config, (MikadoConfiguration, DaijinConfiguration))
    if not hasattr(mikado_config.reference, "genome"):
        raise InvalidJson("Invalid configuration; reference: {}".format(mikado_config.reference))
    try:
        prepare(mikado_config, logger)
        sys.exit(0)
    except Exception:
        raise


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
    parser.add_argument("--fasta", "--reference", dest="reference",
                        type=argparse.FileType(), help="Genome FASTA file. Required.")
    verbosity = parser.add_mutually_exclusive_group()
    verbosity.add_argument("-v", "--verbose", action="store_true", default=False)
    verbosity.add_argument("-q", "--quiet", action="store_true", default=False)
    parser.add_argument("--start-method", dest="start_method",
                        choices=["fork", "spawn", "forkserver"],
                        default=None, help="Multiprocessing start method.")
    strand = parser.add_mutually_exclusive_group()
    strand.add_argument("-s", "--strand-specific", dest="strand_specific",
                        action="store_true", default=False,
                        help="""Flag. If set, monoexonic transcripts
                        will be left on their strand rather than being
                        moved to the unknown strand.""")
    strand.add_argument("-sa", "--strand-specific-assemblies",
                        default=None,
                        type=str,
                        dest="strand_specific_assemblies",
                        help="Comma-delimited list of strand specific assemblies.")
    parser.add_argument("--list", type=argparse.FileType("r"),
                        help="""Tab-delimited file containing rows with the following format:
    <file>  <label> <strandedness(def. False)> <score(optional, def. 0)> <is_reference(optional, def. False)> <exclude_redundant(optional, def. True)> <strip_cds(optional, def. False)> <skip_split(optional, def. False)>
    "strandedness", "is_reference", "exclude_redundant", "strip_cds" and "skip_split" must be boolean values (True, False)
    "score" must be a valid floating number."""
                        )
    parser.add_argument("-l", "--log", type=argparse.FileType("w"), default=None,
                        help="Log file. Optional.")
    parser.add_argument("--lenient", action="store_true", default=None,
                        help="""Flag. If set, transcripts with only non-canonical
                        splices will be output as well.""")
    parser.add_argument("-m", "--minimum-cdna-length", default=None, dest="minimum_cdna_length", type=positive,
                        help="Minimum length for transcripts. Default: 200 bps.")
    parser.add_argument("-MI", "--max-intron-size", default=None, type=positive, dest="max_intron_length",
                        help="Maximum intron length for transcripts. Default: 1,000,000 bps.")
    parser.add_argument("-p", "--procs",
                        help="Number of processors to use (default %(default)s)",
                        type=to_cpu_count, default=None)
    parser.add_argument("-scds", "--strip_cds", action="store_true", default=False,
                        help="Boolean flag. If set, ignores any CDS/UTR segment.")
    parser.add_argument("--labels", type=str, default="",
                        help="""Labels to attach to the IDs of the transcripts of the input files,
                        separated by comma.""")
    parser.add_argument("--codon-table", dest="codon_table", default=None,
                        help="""Codon table to use. Default: 0 (ie Standard, NCBI #1, but only ATG is considered \
    a valid start codon.""")
    parser.add_argument("--single", "--single-thread", action="store_true", default=False,
                        help="Disable multi-threading. Useful for debugging.")
    parser.add_argument("-od", "--output-dir", dest="output_dir",
                        type=str, default=None,
                        help="Output directory. Default: current working directory")
    parser.add_argument("-o", "--out", default=None,
                        help="Output file. Default: mikado_prepared.gtf.")
    parser.add_argument("-of", "--out_fasta", default=None,
                        help="Output file. Default: mikado_prepared.fasta.")
    parser.add_argument("--configuration", "--json-conf", dest="configuration",
                        type=str, default="",
                        help="Configuration file.")
    parser.add_argument("-er", "--exclude-redundant", default=None,
                        dest="exclude_redundant", action="store_true",
                        help="Boolean flag. If invoked, Mikado prepare will exclude redundant models,\
ignoring the per-sample instructions.")
    cds_stripping = parser.add_mutually_exclusive_group()
    cds_stripping.add_argument("--strip-faulty-cds", default=None, action="store_true",
                        help="Flag. If set, transcripts with an incorrect CDS will be retained but \
with their CDS stripped. Default behaviour: the whole transcript will be considered invalid and discarded.")
    parser.add_argument("--seed", type=int, default=None,
                        help="Random seed number.")
    parser.add_argument("gff", help="Input GFF/GTF file(s).", nargs="*")
    parser.set_defaults(func=prepare_launcher)
    return parser
