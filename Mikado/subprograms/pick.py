#!/usr/bin/env python3
# coding: utf-8

"""Launcher of the Mikado pick step."""

import argparse
import sys
import os

# from ..configuration.configurator import check_and_load_scoring
from ..utilities.log_utils import create_default_logger, create_null_logger
import random
from ..utilities import to_region, percentage
from ..utilities import IntervalTree, Interval
from ..configuration.configurator import load_and_validate_config


def check_log_settings(args):

    """
    Quick method to check the consistency of log settings
    from the namespace.
    :param args: a Namespace
    :return: args
    """

    if args.log == "stderr":
        args.configuration.log_settings.log = None
    elif args.log is not None:
        args.configuration.log_settings.log = args.log

    if args.log_level is not None:
        args.configuration.log_settings.log_level = args.log_level
    elif args.verbose is True:
        args.configuration.log_settings.log_level = "DEBUG"
    elif args.noverbose is True:
        args.configuration.log_settings.log_level = "ERROR"

    return args


def check_run_options(args, logger=create_null_logger()):
    """
    Quick method to check the consistency of run option settings
    from the namespace.
    :param args: a Namespace
    :param logger: a logger instance.
    :return: args
    """

    if args.start_method is not None:
        args.configuration.multiprocessing_method = args.start_method

    if args.procs is not None:
        args.configuration.threads = args.procs

    args.configuration.pick.run_options.single_thread = args.single

    if args.seed is not None:
        args.configuration.seed = args.seed

    if args.no_cds is not False:
        args.configuration.pick.run_options.exclude_cds = True
    if args.no_purge is True:
        args.configuration.pick.run_options.purge = False

    if args.flank is not None:
        args.configuration.pick.run_options.flank = args.flank

    if args.output_dir is not None:
        args.configuration.pick.files.output_dir = os.path.abspath(args.output_dir)
    else:
        args.configuration.pick.files.output_dir = os.path.abspath(args.configuration.pick.files.output_dir)

    if args.source is not None:
        args.configuration.pick.output_format.source = args.source
    if args.prefix is not None:
        args.configuration.pick.output_format.id_prefix = args.prefix

    if args.sqlite_db is not None:
        if not os.path.exists(args.sqlite_db):
            logger.critical("Mikado database %s not found. Exiting.", args.sqlite_db)
            sys.exit(1)
        args.configuration.db_settings.db = args.sqlite_db
        args.configuration.db_settings.dbtype = "sqlite"

    elif not (args.configuration.db_settings.dbtype == "sqlite" and
              not os.path.exists(args.configuration.db_settings.db) and
              os.path.abspath(args.configuration.pick.files.output_dir) != os.path.dirname(
                args.configuration.db_settings.db)
    ):
        __compound = os.path.join(args.configuration.pick.files.output_dir,
                                  args.configuration.db_settings.db)
        __base = os.path.join(args.configuration.pick.files.output_dir,
                              args.configuration.db_settings.db)
        if os.path.exists(__compound):
            args.configuration.db_settings.db = __compound
        elif os.path.exists(__base):
            args.configuration.db_settings.db = __base
        else:
            logger.critical("Mikado database {} not found. Exiting.", args.sqlite_db)
            sys.exit(1)

    if args.mode is not None:
        if args.mode == "nosplit":
            args.configuration.pick.chimera_split.execute = False
        else:
            args.configuration.pick.chimera_split.execute = True
            if args.mode == "split":
                args.configuration.pick.chimera_split.blast_check = False
            else:
                args.configuration.pick.chimera_split.blast_check = True
                args.configuration.pick.chimera_split.blast_params.leniency = args.mode.upper()

    if args.pad is not None:
        args.configuration.pick.alternative_splicing.pad = args.pad

    if args.min_clustering_cds_overlap is not None:
        args.configuration.pick.clustering.min_cds_overlap = args.min_clustering_cds_overlap

    if args.min_clustering_cdna_overlap is not None:
        args.configuration.pick.clustering.min_cdna_overlap = args.min_clustering_cdna_overlap
        if args.min_clustering_cds_overlap is None:
            args.configuration.pick.clustering.min_cds_overlap = args.min_clustering_cdna_overlap

    if args.pad_max_splices is not None:
        args.configuration.pick.alternative_splicing.ts_max_splices = True

    if args.pad_max_distance is not None:
        args.configuration.pick.alternative_splicing.ts_distance = True

    if args.intron_range is not None:
        args.configuration.pick.run_options.intron_range = tuple(sorted(args.intron_range))

    if args.max_intron_length is not None:
        args.configuration.prepare.max_intron_length = args.max_intron_length

    if args.cds_only is True:
        args.configuration.pick.clustering.cds_only = True

    if args.as_cds_only is True:
        args.configuration.pick.alternative_splicing.cds_only = True

    if args.gff:
        args.configuration.pick.files.input = args.gff
        if not os.path.exists(args.gff):
            raise ValueError("The input file {} does not exist. Please double check!".format(args.gff))

    prep_gtf = os.path.join(args.configuration.prepare.files.output_dir, args.configuration.prepare.files.out)
    if not os.path.exists(args.configuration.pick.files.input):
        if os.path.exists(prep_gtf):
            args.configuration.pick.files.input = prep_gtf
        elif os.path.exists(args.configuration.prepare.files.out):
            args.configuration.pick.files.input = args.configuration.prepare.files.out
        else:
            msg = "I tried to infer the input file from the prepare option, but failed. Please point me to the correct file \
through the command line or by correcting the configuration file."
            logger.critical(msg)
            raise ValueError(msg)

    if args.loci_out:
        args.configuration.pick.files.loci_out = args.loci_out if \
            args.loci_out.split('.')[-1] in ("gff3", "gff") else \
            "{0}.gff3".format(args.loci_out)

    if args.monoloci_out:
        args.configuration.pick.files.monoloci_out = args.monoloci_out if \
            args.monoloci_out.split('.')[-1] in ("gff3", "gff") else \
            "{0}.gff3".format(args.monoloci_out)

    if args.subloci_out:
        args.configuration.pick.files.subloci_out = args.subloci_out if \
            args.subloci_out.split('.')[-1] in ("gff3", "gff") else \
            "{0}.gff3".format(args.subloci_out)

    if args.log:
        args.configuration.pick.files.log = args.log

    if args.shm is True:
        args.configuration.pick.run_options.shm = True

    if args.only_reference_update is True:
        args.configuration.pick.run_options.only_reference_update = True
        args.configuration.pick.run_options.reference_update = True

    if args.reference_update is True:
        args.configuration.pick.run_options.reference_update = True

    if args.check_references is True:
        args.configuration.pick.run_options.check_references = True

    if args.report_all_orfs is True:
        args.configuration.pick.run_options.report_all_orfs = True

    if getattr(args, "fasta"):
        args.fasta.close()
        args.configuration.reference.genome = args.fasta.name

    if args.scoring_file is not None:
        if not os.path.exists(args.scoring_file) and os.path.isfile(args.scoring_file):
            raise ValueError("Invalid/inexistent scoring file: {}".format(args.scoring_file))
        args.configuration.pick.scoring_file = args.scoring_file

    if (args.configuration.pick.alternative_splicing.pad and
            not os.path.exists(args.configuration.reference.genome)):
        logger.critical("Transcript padding cannot function unless the genome file is specified. \
        Please either provide a valid genome file or disable the padding.")
        sys.exit(1)

    if args.keep_disrupted_cds is True:
        args.configuration.pick.alternative_splicing.keep_cds_disrupted_by_ri = True

    if args.exclude_retained_introns is True:
        args.configuration.pick.alternative_splicing.keep_retained_introns = False

    if args.codon_table not in (False, None, True):
        args.configuration.serialise.codon_table = str(args.codon_table)

    args.configuration = load_and_validate_config(args.configuration, logger=logger)
    return args


def pick(args):

    """
    This function launches the pick step, using the data derived from the Namespace.
    :param args: argparse Namespace with the configuration for the run.

    """

    logger = create_default_logger("pick_init")

    args.configuration.close()
    args.configuration = load_and_validate_config(args.configuration.name, logger=logger)

    try:
        args = check_log_settings(args)
    except Exception as exc:
        logger.error(exc)
        raise exc

    try:
        args = check_run_options(args, logger=logger)
    except Exception as exc:
        logger.error(exc)
        raise exc

    if args.regions is not None:
        regions = dict()
        if os.path.exists(args.regions):
            with open(args.regions) as f_regions:
                for line in f_regions:
                    chrom, start, end = to_region(line)
                    if chrom not in regions:
                        regions[chrom] = IntervalTree()
                    regions[chrom].add_interval(Interval(start, end))
        else:
            chrom, start, end = to_region(args.regions)
            regions[chrom] = IntervalTree.from_intervals([Interval(start, end)])
    else:
        regions = None

    from ..picking import Picker
    creator = Picker(args.configuration, commandline=" ".join(sys.argv), regions=regions)
    creator()
    sys.exit(0)


def pick_parser():
    """
    Parser for the picking step.
    """
    parser = argparse.ArgumentParser(description="Launcher of the Mikado pipeline.")
    parser.add_argument("--fasta", type=argparse.FileType(),
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
    parser.add_argument("--configuration", "--json-conf", dest="configuration",
                        type=argparse.FileType("r"), required=True,
                        help="JSON/YAML configuration file for Mikado.")
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

    parser.add_argument("--no_cds", action="store_true", default=False,
                        help="""Flag. If set, not CDS information will be printed out in the GFF output files.""")

    parser.add_argument("--flank", default=None, type=int,
                        help="""Flanking distance (in bps) to group non-overlapping transcripts into a single \
superlocus. Default: determined by the configuration file.""")
    parser.add_argument("--max-intron-length", default=None, type=int,
                        help="""Maximum intron length for a transcript. Default: inferred from the configuration \
file (default value there is 1,000,000 bps).""")
    parser.add_argument('--no-purge', action='store_true', default=False,
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
    parser.add_argument("-od", "--output-dir", dest="output_dir",
                        type=str, default=None,
                        help="Output directory. Default: current working directory")
    parser.add_argument("--single", action="store_true", default=False,
                        help="""Flag. If set, Creator will be launched with a single process, without involving the
multithreading apparatus. Useful for debugging purposes only.""")
    log_options = parser.add_argument_group("Log options")
    log_options.add_argument("-l", "--log", default=None,
                             help="""File to write the log to.
                             Default: decided by the configuration file.""")
    verbosity = log_options.add_mutually_exclusive_group()
    verbosity.add_argument("-v", "--verbose", action="store_true",
                           default=False, help="Flag. If set, the debug mode will be activated.")
    verbosity.add_argument("-nv", "--noverbose", action="store_true",
                           default=False, help="Flag. If set, the log will report only errors and critical events.")
    log_options.add_argument("-lv", "--log-level", dest="log_level",
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
    parser.add_argument("--seed", type=int, default=None,
                        help="Random seed number.")
    # parser.formatter_class = argparse.HelpFormatter
    parser.add_argument("gff", nargs="?", default=None)
    parser.set_defaults(func=pick)
    return parser
