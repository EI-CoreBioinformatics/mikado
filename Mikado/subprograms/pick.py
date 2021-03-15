#!/usr/bin/env python3
# coding: utf-8

"""Launcher of the Mikado pick step."""

import argparse
import re
import sys
import os
from typing import Union, Dict
from ._utils import check_log_settings_and_create_logger, _set_pick_mode
import marshmallow
from ..configuration import DaijinConfiguration, MikadoConfiguration
from ..exceptions import InvalidConfiguration
from ..utilities.log_utils import create_default_logger, create_null_logger
from ..utilities import to_region, percentage
from ..utilities import IntervalTree, Interval
from ..configuration.configurator import load_and_validate_config
from ..picking import Picker


def _parse_regions(regions_string: Union[None,str]) -> Union[None, Dict[str, IntervalTree]]:
    if regions_string is None:
        return None

    regions = dict()
    if os.path.exists(regions_string):
        with open(regions_string) as f_regions:
            for counter, line in enumerate(f_regions, start=1):
                try:
                    chrom, start, end = to_region(line)
                except ValueError:
                    raise ValueError(f"Invalid region line, no. {counter}: {line}")
                if chrom not in regions:
                    regions[chrom] = IntervalTree()
                regions[chrom].add(Interval(start, end))
    else:
        chrom, start, end = to_region(regions_string)
        regions[chrom] = IntervalTree.from_intervals([Interval(start, end)])

    return regions


def _set_pick_output_options(conf: Union[DaijinConfiguration, MikadoConfiguration], args,
                             logger=create_null_logger()) -> Union[DaijinConfiguration, MikadoConfiguration]:
    conf.pick.output_format.source = args.source if args.source is not None else conf.pick.output_format.source
    conf.pick.output_format.id_prefix = args.prefix if args.prefix is not None else conf.pick.output_format.id_prefix
    conf.pick.output_format.report_all_external_metrics = True if args.report_all_external_metrics is True else \
        conf.pick.output_format.report_all_external_metrics
    conf.pick.output_format.report_all_orfs = True if args.report_all_orfs is True else \
        conf.pick.output_format.report_all_orfs
    conf.pick.files.log = args.log if args.log else conf.pick.files.log
    pat = re.compile(r"\.(gff3|gff)")
    if args.loci_out:
        conf.pick.files.loci_out = args.loci_out if pat.search(args.loci_out) else "{0}.gff3".format(args.loci_out)

    if args.monoloci_out:
        conf.pick.files.monoloci_out = args.monoloci_out if pat.search(args.monoloci_out) else "{0}.gff3".format(
            args.monoloci_out)

    if args.subloci_out:
        conf.pick.files.subloci_out = args.subloci_out if pat.search(args.subloci_out) else "{0}.gff3".format(
            args.subloci_out)

    return conf


def _set_pick_run_options(conf: Union[DaijinConfiguration, MikadoConfiguration], args,
                          logger=create_null_logger()) -> Union[DaijinConfiguration, MikadoConfiguration]:
    conf.pick.run_options.single_thread = args.single
    conf.pick.run_options.exclude_cds = True if args.no_cds is True else conf.pick.run_options.exclude_cds
    conf.pick.run_options.intron_range = tuple(sorted(args.intron_range)) if args.intron_range is not None \
        else conf.pick.run_options.intron_range
    conf.pick.run_options.shm = True if args.shm is not None else conf.pick.run_options.shm
    if args.only_reference_update is True:
        conf.pick.run_options.only_reference_update = True
        conf.pick.run_options.reference_update = True
    conf.pick.run_options.reference_update = True if args.reference_update is True else \
        conf.pick.run_options.reference_update
    conf.pick.run_options.check_references = True if args.check_references is True else \
        conf.pick.run_options.check_references

    return conf


def _set_pick_clustering_options(conf: Union[DaijinConfiguration, MikadoConfiguration],
                                 args) -> Union[DaijinConfiguration, MikadoConfiguration]:

    conf.pick.clustering.purge = False if args.no_purge is True else conf.pick.clustering.purge
    conf.pick.clustering.flank = args.flank if args.flank is not None else conf.pick.clustering.flank
    conf.pick.clustering.min_cds_overlap = args.min_clustering_cds_overlap if \
        args.min_clustering_cds_overlap else conf.pick.clustering.min_cds_overlap
    conf.pick.clustering.cds_only = True if args.cds_only else conf.pick.clustering.cds_only
    if args.min_clustering_cdna_overlap is not None:
        conf.pick.clustering.min_cdna_overlap = args.min_clustering_cdna_overlap
        if args.min_clustering_cds_overlap is None:
            conf.pick.clustering.min_cds_overlap = args.min_clustering_cdna_overlap

    return conf


def _set_pick_as_options(conf: Union[DaijinConfiguration, MikadoConfiguration],
                         args) -> Union[DaijinConfiguration, MikadoConfiguration]:

    conf.pick.alternative_splicing.pad = args.pad if args.pad is True else \
        conf.pick.alternative_splicing.pad
    conf.pick.alternative_splicing.ts_max_splices = True if args.pad_max_splices \
        else conf.pick.alternative_splicing.ts_max_splices
    conf.pick.alternative_splicing.ts_distance = True if args.pad_max_distance is not None else \
        conf.pick.alternative_splicing.ts_distance
    conf.pick.alternative_splicing.cds_only = True if args.as_cds_only is True else \
        conf.pick.alternative_splicing.cds_only
    conf.pick.alternative_splicing.keep_cds_disrupted_by_ri = True if args.keep_disrupted_cds is True \
        else conf.pick.alternative_splicing.keep_cds_disrupted_by_ri
    conf.pick.alternative_splicing.keep_retained_introns = False if args.exclude_retained_introns is True else \
        conf.pick.alternative_splicing.keep_retained_introns

    return conf


def _set_conf_values_from_args(conf: Union[DaijinConfiguration, MikadoConfiguration], args,
                               logger=create_null_logger()) -> Union[DaijinConfiguration, MikadoConfiguration]:

    conf.multiprocessing_method = args.start_method if args.start_method else conf.multiprocessing_method
    conf.threads = args.procs if args.procs is not None else conf.threads
    if args.random_seed is True:
        conf.seed = None
    elif args.seed is not None:
        conf.seed = args.seed
    else:
        pass

    conf.pick.scoring_file = args.scoring_file if args.scoring_file is not None else conf.pick.scoring_file

    conf.prepare.max_intron_length = args.max_intron_length if args.max_intron_length is not None else \
        conf.prepare.max_intron_length
    conf.serialise.codon_table = str(args.codon_table) if args.codon_table not in (False, None, True) \
        else conf.serialise.codon_table     

    conf = _set_pick_output_options(conf, args)
    conf = _set_pick_mode(conf, args.mode)
    conf = _set_pick_run_options(conf, args)
    conf = _set_pick_clustering_options(conf, args)
    conf = _set_pick_as_options(conf, args)

    try:
        conf = load_and_validate_config(conf, logger=logger)
    except marshmallow.exceptions.MarshmallowError as exc:
        logger.critical("Invalid options specified for the configuration: {}".format(exc))
        raise exc
    return conf


def _check_db(conf: Union[MikadoConfiguration, DaijinConfiguration], args,
              logger=create_null_logger()) -> Union[MikadoConfiguration, DaijinConfiguration]:

    logger.debug("Checking the database")
    if args.sqlite_db is not None:
        if not os.path.exists(args.sqlite_db):
            exc = InvalidConfiguration(f"Mikado database {args.sqlite_db} not found. Exiting.")
            logger.critical(exc)
            raise exc
        logger.debug(f"Setting the database from the CLI to {args.sqlite_db}")
        conf.db_settings.db = args.sqlite_db
        conf.db_settings.dbtype = "sqlite"

    if conf.db_settings.dbtype == "sqlite":
        raw = conf.db_settings.db
        db_basename = os.path.basename(conf.db_settings.db)
        __compound = os.path.join(conf.pick.files.output_dir, db_basename)
        __base = os.path.join(conf.pick.files.output_dir, db_basename)
        found = False
        for option in raw, __compound, __base:
            if os.path.exists(option):
                conf.db_settings.db = option
                found = True
                break
        if found is False:
            exc = InvalidConfiguration(f"Mikado database {conf.db_settings.db} not found. Exiting.")
            logger.critical(exc)
            raise exc

    logger.debug(f"Found database: {conf.db_settings.dbtype}:///{conf.db_settings.db}")
    return conf


def _check_pick_input(conf: Union[MikadoConfiguration, DaijinConfiguration], args,
                      logger=create_null_logger()) -> Union[MikadoConfiguration, DaijinConfiguration]:

    if args.gff:
        conf.pick.files.input = args.gff
        if not os.path.exists(args.gff):
            raise InvalidConfiguration("The input file {} does not exist. Please double check!".format(args.gff))

    prep_gtf = os.path.join(conf.prepare.files.output_dir, conf.prepare.files.out)
    if not os.path.exists(conf.pick.files.input):
        if os.path.exists(prep_gtf):
            conf.pick.files.input = prep_gtf
        elif os.path.exists(conf.prepare.files.out):
            conf.pick.files.input = conf.prepare.files.out
        else:
            exc = InvalidConfiguration("I tried to infer the input file from the prepare option, but failed. Please "
                                       "point me to the correct file through the command line or by correcting the "
                                       "configuration file.")
            logger.critical(exc)
            raise exc

    if args.genome:
        if not os.path.exists(args.genome):
            raise InvalidConfiguration(f"The requested genome FASTA file does not seem to exist: {args.genome}")
        conf.reference.genome = args.genome

    if conf.pick.alternative_splicing.pad and not os.path.exists(conf.reference.genome):
        exc = InvalidConfiguration("Transcript padding cannot function unless the genome file is specified. \
Please either provide a valid genome file or disable the padding.")
        logger.critical(exc)
        raise exc

    conf = _check_db(conf, args, logger)
    return conf


def check_run_options(mikado_configuration: Union[MikadoConfiguration, DaijinConfiguration],
                      args: argparse.Namespace, logger=create_null_logger()):
    """
    Quick method to check the consistency of run option settings
    from the namespace.
    :param args: a Namespace
    :param logger: a logger instance.
    :return: args
    """

    mikado_configuration = _set_conf_values_from_args(mikado_configuration, args, logger=logger)
    mikado_configuration = _check_pick_input(mikado_configuration, args, logger)
    mikado_configuration = load_and_validate_config(mikado_configuration, logger=logger)
    return mikado_configuration


def pick(args):

    """
    This function launches the pick step, using the data derived from the Namespace.
    :param args: argparse Namespace with the configuration for the run.

    """

    logger = create_default_logger("pick", level="WARNING")
    mikado_configuration = load_and_validate_config(args.configuration, logger=logger)
    # Create the output directory. Necessary to do it here to avoid the logger being put in the wrong place.
    if args.output_dir is not None:
        mikado_configuration.pick.files.output_dir = os.path.abspath(args.output_dir)
    else:
        mikado_configuration.pick.files.output_dir = os.path.abspath(mikado_configuration.pick.files.output_dir)
    try:
        os.makedirs(mikado_configuration.pick.files.output_dir, exist_ok=True)
    except OSError:
        exc = OSError("I cannot create the output directory {}. Aborting.".format(
            mikado_configuration.pick.files.output_dir))
        logger.critical(exc)
        raise exc
    mikado_configuration, logger = check_log_settings_and_create_logger(mikado_configuration, args.log, args.log_level,
                                                                        section="pick")
    mikado_configuration = check_run_options(mikado_configuration, args, logger=logger)
    regions = _parse_regions(args.regions)
    creator = Picker(mikado_configuration, commandline=" ".join(sys.argv), regions=regions)
    creator()


def pick_parser():
    """
    Parser for the picking step.
    """
    parser = argparse.ArgumentParser(description="Launcher of the Mikado pipeline.")
    parser.add_argument("--fasta", "--genome", default=None, dest="genome",
                        help="Genome FASTA file. Required for transcript padding.")
    parser.add_argument("--start-method", dest="start_method",
                        choices=["fork", "spawn", "forkserver"],
                        default=None, help="Multiprocessing start method.")
    parser.add_argument("--shm", default=False, action="store_true",
                        help="Flag. If switched, Mikado pick will copy the database to RAM (ie SHM) for faster access \
during the run.")
    parser.add_argument("-p", "--procs", type=int, default=None,
                        help="""Number of processors to use. \
Default: look in the configuration file (1 if undefined)""")
    parser.add_argument("--configuration", "--json-conf", dest="configuration", required=True,
                        help="Configuration file for Mikado.")
    parser.add_argument("--scoring-file", dest="scoring_file",
                        type=str, default=None,
                        required=False,
                        help="Optional scoring file for the run. It will override the value set in the configuration.")
    parser.add_argument("-i", "--intron-range",
                        dest="intron_range", type=int, nargs=2,
                        default=None,
                        help="""Range into which intron lengths should fall, as a couple of integers. \
Transcripts with intron lengths outside of this range will be penalised. Default: (60, 900)""")
    padding = parser.add_mutually_exclusive_group()
    padding.add_argument("--no-pad", dest="pad", default=None, action="store_false", help="Disable transcript padding.")
    padding.add_argument("--pad", default=None,
                         action="store_true",
                         help="Whether to pad transcripts in loci.")
    padding.add_argument("--codon-table", dest="codon_table", default=None,
                         help="""Codon table to use. Default: 0 (ie Standard, NCBI #1, but only ATG is considered \
        a valid start codon.""")
    parser.add_argument("--pad-max-splices", default=None, dest="pad_max_splices",
                        type=int, help="Maximum splice sites that can be crossed during transcript padding.")
    parser.add_argument("--pad-max-distance", default=None, dest="pad_max_distance",
                        type=int, help="Maximum amount of bps that transcripts can be padded with (per side).")
    parser.add_argument("-r", "--regions",
                        help="""Either a single region on the CLI or a file listing a series of target regions.
Mikado pick will only consider regions included in this string/file.
Regions should be provided in a WebApollo-like format: <chrom>:<start>..<end>""")
    output = parser.add_argument_group("Options related to the output files.")
    output.add_argument("-od", "--output-dir", dest="output_dir",
                        type=str, default=None,
                        help="Output directory. Default: current working directory")
    output.add_argument("--subloci-out", type=str, default=None, dest="subloci_out")
    output.add_argument("--monoloci-out", type=str, default=None, dest="monoloci_out")
    output.add_argument("--loci-out", type=str, default=None, dest="loci_out",
                        help="""This output file is mandatory.
                        If it is not specified in the configuration file,
                        it must be provided here.""")
    output.add_argument("--prefix", type=str, default=None,
                        help="Prefix for the genes. Default: Mikado")
    output.add_argument('--source', type=str, default=None,
                        help='Source field to use for the output files.')
    output.add_argument("--report-all-external-metrics", default=None,
                        action="store_true",
                        help="Boolean switch. If activated, Mikado will report all available external metrics, not just \
those requested for in the scoring configuration. This might affect speed in Minos analyses.")
    parser.add_argument("--no_cds", action="store_true", default=False,
                        help="""Flag. If set, not CDS information will be printed out in the GFF output files.""")

    parser.add_argument("--flank", default=None, type=int,
                        help="""Flanking distance (in bps) to group non-overlapping transcripts into a single \
superlocus. Default: determined by the configuration file.""")
    parser.add_argument("--max-intron-length", default=None, type=int,
                        help="""Maximum intron length for a transcript. Default: inferred from the configuration \
file (default value there is 1,000,000 bps).""")
    parser.add_argument('--no-purge', action='store_true', default=False,
                        dest="no_purge",
                        help='''Flag. If set, the pipeline will NOT suppress any loci \
whose transcripts do not pass the requirements set in the JSON file.''')
    parser.add_argument("--cds-only", dest="cds_only",
                        default=None, action="store_true",
                        help=""""Flag. If set, Mikado will only look for overlap in the coding features \
when clustering transcripts (unless one transcript is non-coding, in which case  the whole transcript will \
be considered). Please note that Mikado will only consider the **best** ORF for this. \
Default: False, Mikado will consider transcripts in their entirety.""")
    parser.add_argument("--as-cds-only", dest="as_cds_only", default=None, action="store_true",
                        help="""Flag. If set, Mikado will only consider the CDS to determine whether a transcript
                        is a valid alternative splicing event in a locus.""")
    parser.add_argument("--reference-update", dest="reference_update", default=None,
                         action="store_true",
                         help="""Flag. If switched on, Mikado will prioritise transcripts marked as reference and will \
    consider any other transcipt within loci only in reference to these reference transcripts. Novel loci will still be reported.""")
    parser.add_argument("--report-all-orfs", default=False, action="store_true",
                        help="Boolean switch. If set to true, all ORFs will be reported, not just the primary.")
    parser.add_argument("--only-reference-update", dest="only_reference_update", default=None,
                        action="store_true",
                        help="""Flag. If switched on, Mikado will only keep loci where at least one of the transcripts \
is marked as "reference". CAUTION: if no transcript has been marked as reference, the output will be completely empty!""")
    parser.add_argument("-eri", "--exclude-retained-introns", default=None, action="store_true",
                        help="""Exclude all retained intron alternative splicing events from the final output. \
Default: False. Retained intron events that do not dirsupt the CDS are kept by Mikado in the final output.""")
    parser.add_argument("-kdc", "--keep-disrupted-cds", default=None, action="store_true",
                        help="""Keep in the final output transcripts whose CDS is most probably disrupted by a \
retained intron event. Default: False. Mikado will try to detect these instances and exclude them from the \
final output.""")
    parser.add_argument("-mco", "--min-clustering-cdna-overlap", default=None, type=percentage,
                         help="Minimum cDNA overlap between two transcripts for them to be considered part of the same \
locus during the late picking stages. \
NOTE: if --min-cds-overlap is not specified, it will be set to this value! \
Default: 20%%.")
    parser.add_argument("-mcso", "--min-clustering-cds-overlap", default=None, type=percentage,
                         help="Minimum CDS overlap between two transcripts for them to be considered part of the same \
locus during the late picking stages. \
NOTE: if not specified, and --min-cdna-overlap is specified on the command line, min-cds-overlap will be set to this value! \
Default: 20%%.")
    parser.add_argument("--check-references", dest="check_references", default=None,
                        action="store_true",
                        help="""Flag. If switched on, Mikado will also check reference models against the general
transcript requirements, and will also consider them as potential fragments. This is useful in the context of e.g.
updating an *ab-initio* results with data from RNASeq, protein alignments, etc. 
""")
    parser.add_argument("-db", "--sqlite-db", dest="sqlite_db",
                        default=None, type=str,
                        help="Location of an SQLite database to overwrite what is specified \
in the configuration file.")
    parser.add_argument("--single", action="store_true", default=False,
                        help="""Flag. If set, Creator will be launched with a single process, without involving the
multithreading apparatus. Useful for debugging purposes only.""")
    log_options = parser.add_argument_group("Log options")
    log_options.add_argument("-l", "--log", default=None,
                             help="""File to write the log to.
                             Default: decided by the configuration file.""")
    verbosity = log_options.add_mutually_exclusive_group()
    verbosity.add_argument("--verbose", default=None, dest="log_level", action="store_const", const="DEBUG")
    verbosity.add_argument("--quiet", default=None, dest="log_level", action="store_const", const="WARNING")
    verbosity.add_argument("-lv", "--log-level", dest="log_level",
                           choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], default=None,
                           help="Logging level. Default: retrieved by the configuration file.")
    # parser.formatter_class = argparse.RawTextHelpFormatter
    parser.add_argument("--mode", default=None,
                        choices=["nosplit", "stringent", "lenient", "permissive", "split"],
                        help="""Mode in which Mikado will treat transcripts with multiple ORFs.
                        - nosplit: keep the transcripts whole.
                        - stringent: split multi-orf transcripts if two consecutive ORFs have both BLAST hits
                        and none of those hits is against the same target.
                        - lenient: split multi-orf transcripts as in stringent, and additionally, also when
                         either of the ORFs lacks a BLAST hit (but not both).
                        - permissive: like lenient, but also split when both ORFs lack BLAST hits
                        - split: split multi-orf transcripts regardless of what BLAST data is available.""")
    seed_group = parser.add_mutually_exclusive_group()
    seed_group.add_argument("--seed", type=int, default=None, help="Random seed number. Default: 0.")
    seed_group.add_argument("--random-seed", action="store_true", default=False,
                            help="Generate a new random seed number (instead of the default of 0)")
    parser.add_argument("gff", nargs="?", default=None)
    parser.set_defaults(func=pick)
    return parser
