#!/usr/bin/env python3

"""Stub of pre-configurer for Mikado"""
import yaml
import os
import re
from pkg_resources import resource_filename, resource_stream
import glob
import argparse
import sys
import dataclasses

from ._utils import _set_pick_mode
from ..configuration import DaijinConfiguration, MikadoConfiguration
from ..exceptions import InvalidConfiguration
from ..utilities import comma_split, percentage, merge_dictionaries
from ..utilities.namespace import Namespace
from ..configuration.configurator import load_and_validate_config
from ..configuration import print_config
import tempfile
from ..utilities.log_utils import create_null_logger, create_default_logger
try:
    from yaml import CSafeLoader as yLoader
except ImportError:
    from yaml import SafeLoader as yLoader
from .prepare import parse_prepare_options

__author__ = 'Luca Venturini'


def __add_daijin_specs(args, config):
    from ..configuration.daijin_configurator import create_cluster_config, create_daijin_config
    namespace = Namespace(default=False)
    namespace.r1 = []
    namespace.r2 = []
    namespace.samples = []
    namespace.strandedness = []
    namespace.asm_methods = []
    namespace.long_aln_methods = []
    namespace.aligners = []
    if args.mode is not None:
        namespace.modes = args.mode[:]
    else:
        namespace.modes = ["stringent"]
    if args.blast_targets:
        namespace.prot_db = args.blast_targets[:]
    else:
        namespace.prot_db = []
    namespace.cluster_config = None
    namespace.scheduler = ""
    namespace.flank = None
    namespace.intron_range = None
    namespace.genome = args.reference
    namespace.transcriptome = ""
    namespace.name = "Daijin"
    namespace.out_dir = args.out_dir if args.out_dir else "Daijin"
    namespace.threads = args.threads
    namespace.scoring = args.scoring
    namespace.new_scoring = getattr(args, "new_scoring", None)
    namespace.full = args.full
    namespace.no_files = args.no_files
    config = create_daijin_config(namespace, config, level="ERROR", piped=True)
    config.blastx.chunks = args.blast_chunks
    config.mikado.use_diamond = (not args.use_blast)
    config.mikado.use_prodigal = (not args.use_transdecoder)
    config.scheduler = args.scheduler
    create_cluster_config(config, args, create_null_logger())
    return config


def create_config(args):
    """
    Utility to create a default configuration file.
    :param args:
    :return:
    """

    if isinstance(args.seed, bool) or args.seed is None:
        raise OSError("Invalid seed: {}".format(args.seed))
    elif not (isinstance(args.seed, int) and 0 <= args.seed <= 2 ** 32 - 1):
        raise OSError("Invalid seed: {}".format(args.seed))

    args.daijin = args.daijin or [isinstance(args.mode, list) and len(args.mode) > 1]
    config = DaijinConfiguration() if args.daijin is True else MikadoConfiguration()

    if args.external is not None:
        other = dataclasses.asdict(load_and_validate_config(args.external, external=True))
        config = dataclasses.asdict(config)
        config = merge_dictionaries(config, other)
        config = load_and_validate_config(config)

    if args.daijin is True:
        config = __add_daijin_specs(args, config)

    config.pick.files.subloci_out = args.subloci_out if args.subloci_out else ""
    config.pick.files.monoloci_out = args.monoloci_out if args.monoloci_out else ""

    if isinstance(args.gff, str):
        args.gff = args.gff.split(",")
    elif not args.gff:
        args.gff = []
    config = parse_prepare_options(args, config)

    if args.random_seed is True:
        config.seed = None
    else:
        config.seed = args.seed

    config.serialise.files.junctions = args.junctions if args.junctions is not None else \
        config.serialise.files.junctions

    config.serialise.files.blast_targets = args.blast_targets if args.blast_targets is not None else \
        config.serialise.files.blast_targets

    if args.only_reference_update is True or args.reference_update is True:
        if len(config.prepare.files.reference) == 0:
            logger = create_default_logger("configure")
            logger.error(
                "No reference dataset provided! Please correct the issue or remove the \"--only-reference-update\" \
switch.")
            sys.exit(1)
        else:
            config.pick.run_options.only_reference_update = True

    config.pick.run_options.check_references = True if args.check_references is True else \
        config.pick.run_options.check_references

    config.pick.output_format.report_all_orfs = True if args.report_all_orfs is True else \
        config.pick.output_format.report_all_orfs

    config.pick.output_format.report_all_external_metrics = True if args.report_all_external_metrics else \
        config.pick.output_format.report_all_external_metrics

    if args.scoring is not None:
        if args.copy_scoring is not False:
            with open(args.copy_scoring, "wt") as out:
                with resource_stream("Mikado", os.path.join("configuration",
                                                            "scoring_files",
                                                            args.scoring)) as original:
                    for line in original:
                        print(line.decode().rstrip(), file=out)
            args.scoring = args.copy_scoring

        config.pick.scoring_file = args.scoring

    config.pick.clustering.cds_only = True if args.cds_only is True else config.pick.clustering.cds_only
    config.pick.alternative_splicing.cds_only = True if args.as_cds_only else config.pick.alternative_splicing.cds_only

    if args.daijin is False and args.mode is not None and len(args.mode) == 1:
        mode = args.mode.pop()
        config = _set_pick_mode(config, mode)

    if args.skip_split:
        if not all(_ in config.prepare.files.labels for _ in args.skip_split):
            raise InvalidConfiguration("Some of the labels to skip for splitting are invalid: {}".format(
                [_ for _ in args.skip_split if _ not in config.prepare.files.labels]
            ))
        config.pick.chimera_split.skip = list(set(config.pick.chimera_split.skip.extend(args.skip_split)))

    config.pick.alternative_splicing.pad = args.pad if args.pad is not None else config.pick.alternative_splicing.pad
    config.pick.clustering.min_cds_overlap = args.min_clustering_cds_overlap if args.min_clustering_cds_overlap is not \
        None else config.pick.clustering.min_cds_overlap

    if args.min_clustering_cdna_overlap is not None:
        config.pick.clustering.min_cdna_overlap = args.min_clustering_cdna_overlap
        if args.min_clustering_cds_overlap is None:
            config.pick.clustering.min_cds_overlap = args.min_clustering_cdna_overlap

    config.pick.run_options.intron_range = sorted(args.intron_range) if args.intron_range is not None else \
        config.pick.run_options.intron_range

    config.serialise.codon_table = str(args.codon_table) if args.codon_table not in (None, False, True) else \
        config.serialise.codon_table

    config.pick.alternative_splicing.keep_cds_disrupted_by_ri = True if args.keep_disrupted_cds is True else \
        config.pick.alternative_splicing.keep_cds_disrupted_by_ri

    config.pick.alternative_splicing.keep_retained_introns = False if args.exclude_retained_introns is True else \
        config.pick.alternative_splicing.keep_retained_introns

    if args.out_dir:
        config.prepare.files.output_dir = args.out_dir
        config.serialise.files.output_dir = args.out_dir
        config.pick.files.output_dir = args.out_dir

    config.check()

    # Check that the configuration file is correct
    with tempfile.NamedTemporaryFile("wt", suffix=".json", delete=True) as tempcheck:
        print_config(config, tempcheck, full=args.full, output_format="json")
        tempcheck.flush()
        try:
            load_and_validate_config(tempcheck.name)
        except InvalidConfiguration as exc:
            raise InvalidConfiguration("Created an invalid configuration file! Error:\n{}".format(exc))

    # Print out the final configuration file
    if args.json is True or args.out.name.endswith("json"):
        format_name = "json"
    elif args.yaml is True or args.out.name.endswith("yaml"):
        format_name = "yaml"
    else:
        format_name = "toml"

    print_config(config, args.out, output_format=format_name, no_files=args.no_files, full=args.full)


def configure_parser():
    """
    Parser for the configuration utility.
    :return: the parser.
    :rtype: argparse.ArgumentParser
    """

    scoring_folder = resource_filename("Mikado", os.path.join("configuration", "scoring_files"))
    trailing = re.compile(r"^{}".format(os.path.sep))
    fold_path = re.compile(scoring_folder)

    scoring_files = [trailing.sub(r"", re.sub(fold_path, r"", fname))
                     for fname in glob.iglob(os.path.join(scoring_folder, "**", "*yaml"), recursive=True)]

    parser = argparse.ArgumentParser(description="Configuration utility for Mikado")
    parser.add_argument("--full", action="store_true", default=False)
    seed_group = parser.add_mutually_exclusive_group()
    seed_group.add_argument("--seed", type=int, default=0, help="Random seed number. Default: 0.")
    seed_group.add_argument("--random-seed", action="store_true", default=False,
                            help="Generate a new random seed number (instead of the default of 0)")
    preparer = parser.add_argument_group("Options related to the prepare stage.")
    preparer.add_argument("--minimum-cdna-length", default=None, type=int, dest="minimum_cdna_length",
                          help="Minimum cDNA length for transcripts.")
    preparer.add_argument("--max-intron-length", default=None, type=int, dest="max_intron_length",
                          help="Maximum intron length for transcripts.")
    preparer.add_argument("--strip-faulty-cds", default=None, action="store_true",
                          help="Flag. If set, transcripts with an incorrect CDS will be retained but \
with their CDS stripped. Default behaviour: the whole transcript will be considered invalid and discarded.")
    scoring = parser.add_argument_group("Options related to the scoring system")
    scoring.add_argument("--scoring", type=str, default=None,
                         help="Scoring file to use. Mikado provides the following:\n{}".format(
                             ",\n".join(scoring_files)))
    scoring.add_argument("--copy-scoring", default=False,
                         type=str, dest="copy_scoring",
                         help="File into which to copy the selected scoring file, for modification.")
    picking = parser.add_argument_group("Options related to the picking")
    picking.add_argument("-i", "--intron-range",
                         dest="intron_range", type=int, nargs=2,
                         default=None,
                         help="""Range into which intron lengths should fall, as a couple of integers.
                             Transcripts with intron lengths outside of this range will be penalised.
                             Default: (60, 900)""")
    picking.add_argument("--subloci-out", default="", dest="subloci_out",
                         help="Name of the optional subloci output. By default, this will not be produced.")
    picking.add_argument("--monoloci-out", default="", dest="monoloci_out",
                         help="Name of the optional monoloci output. By default, this will not be produced.")
    picking.add_argument("--no-pad", dest="pad", default=None,
                         action="store_false", help="Disable transcript padding. On by default.")
    picking.add_argument("--reference-update", dest="reference_update", default=None,
                         action="store_true",
                         help="""Flag. If switched on, Mikado will prioritise transcripts marked as reference and will \
consider any other transcipt within loci only in reference to these reference transcripts. Novel loci will still be reported.""")
    picking.add_argument("--only-reference-update", dest="only_reference_update", default=None,
                         action="store_true",
                         help="""Flag. If switched on, Mikado will only keep loci where at least one of the transcripts \
    is marked as "reference". CAUTION: if no transcript has been marked as reference, \
    the output will be completely empty!""")
    picking.add_argument("-eri", "--exclude-retained-introns", default=None, action="store_true",
                         help="""Exclude all retained intron alternative splicing events from the final output. \
Default: False. Retained intron events that do not dirsupt the CDS are kept by Mikado in the final output.""")
    picking.add_argument("-kdc", "--keep-disrupted-cds", default=None, action="store_true",
                         help="""Keep in the final output transcripts whose CDS is most probably disrupted by a \
retained intron event. Default: False. Mikado will try to detect these instances and exclude them from the \
final output.""")
    picking.add_argument("--check-references", dest="check_references", default=None,
                         action="store_true",
                         help="""Flag. If switched on, Mikado will also check reference models against the general
    transcript requirements, and will also consider them as potential fragments. This is useful in the context of e.g.
    updating an *ab-initio* results with data from RNASeq, protein alignments, etc. 
    """)
    picking.add_argument("-mco", "--min-clustering-cdna-overlap", default=None, type=percentage,
                         help="Minimum cDNA overlap between two transcripts for them to be considered part of the same \
locus during the late picking stages. \
NOTE: if --min-cds-overlap is not specified, it will be set to this value! \
Default: 20%%.")
    picking.add_argument("-mcso", "--min-clustering-cds-overlap", default=None, type=percentage,
                         help="Minimum CDS overlap between two transcripts for them to be considered part of the same \
locus during the late picking stages. \
NOTE: if not specified, and --min-cdna-overlap is specified on the command line, min-cds-overlap will be set to this value! \
Default: 20%%.")
    picking.add_argument("--report-all-orfs", default=False, action="store_true",
                         help="Boolean switch. If set to true, all ORFs will be reported, not just the primary.")
    picking.add_argument("--report-all-external-metrics", default=None,
                         action="store_true",
                         help="Boolean switch. If activated, Mikado will report all available external metrics, not just \
those requested for in the scoring configuration. This might affect speed in Minos analyses.")
    picking.add_argument("--cds-only", dest="cds_only",
                         default=None, action="store_true",
                         help=""""Flag. If set, Mikado will only look for overlap in the coding features \
    when clustering transcripts (unless one transcript is non-coding, in which case  the whole transcript will \
    be considered). Please note that Mikado will only consider the **best** ORF for this. \
    Default: False, Mikado will consider transcripts in their entirety.""")
    picking.add_argument("--as-cds-only", dest="as_cds_only", default=None, action="store_true",
                         help="""Flag. If set, Mikado will only consider the CDS to determine whether a transcript
                            is a valid alternative splicing event in a locus.""")
    parser.add_argument("--strand-specific", default=False,
                        action="store_true",
                        help="""Boolean flag indicating whether all the assemblies are strand-specific.""")
    files = parser.add_mutually_exclusive_group()
    files.add_argument("--no-files", dest="no_files",
                       help="""Remove all files-specific options from the printed configuration file.
                       Invoking the "--gff" option will disable this flag.""",
                       default=False, action="store_true")
    files.add_argument("--gff", help="Input GFF/GTF file(s), separated by comma", type=str,
                       default="")
    files.add_argument("--list", type=argparse.FileType("r"),
                        help="""Tab-delimited file containing rows with the following format:
    <file>  <label> <strandedness(def. False)> <score(optional, def. 0)> <is_reference(optional, def. False)> <exclude_redundant(optional, def. True)> <strip_cds(optional, def. False)> <skip_split(optional, def. False)>
    "strandedness", "is_reference", "exclude_redundant", "strip_cds" and "skip_split" must be boolean values (True, False)
    "score" must be a valid floating number.
    """)
    parser.add_argument("--reference", "--genome", help="Fasta genomic reference.", default=None, dest="reference")
    serialisers = parser.add_argument_group(
        "Options related to the serialisation step")
    serialisers.add_argument("--junctions", type=comma_split, default=[])
    serialisers.add_argument("-bt", "--blast_targets", type=comma_split, default=None)
    parser.add_argument("--strand-specific-assemblies", type=str, default="",
                        dest="strand_specific_assemblies",
                        help="""List of strand-specific assemblies among the inputs.""")
    parser.add_argument("--labels", type=str, default="",
                        help="""Labels to attach to the IDs of the transcripts of the input files,
        separated by comma.""")
    parser.add_argument("--codon-table", dest="codon_table", default=None,
                        help="""Codon table to use. Default: 0 (ie Standard, NCBI #1, but only ATG is considered \
    a valid start codon.""")
    parser.add_argument("--external", help="""External configuration file to overwrite/add values from.
    Parameters specified on the command line will take precedence over those present in the configuration file.""")
    daijin = parser.add_argument_group("Options related to configuring a Daijin run.")
    daijin.add_argument("--daijin", action="store_true", default=False,
                        help="Flag. If set, the configuration file will be also valid for Daijin.")
    daijin.add_argument("-bc", "--blast-chunks", type=int, default=10, dest="blast_chunks",
                        help="Number of parallel DIAMOND/BLAST jobs to run. Default: %(default)s.")
    daijin.add_argument("--use-blast", action="store_true", default=False, dest="use_blast",
                        help="Flag. If switched on, Mikado will use BLAST instead of DIAMOND.")
    daijin.add_argument("--use-transdecoder", action="store_true", default=False, dest="use_transdecoder",
                        help="Flag. If switched on, Mikado will use TransDecoder instead of Prodigal.")
    daijin.add_argument("--mode", default=["stringent"], nargs="+",
                        choices=["nosplit", "stringent", "lenient", "permissive", "split"],
                        help="""Mode(s) in which Mikado will treat transcripts with multiple ORFs.
- nosplit: keep the transcripts whole.
- stringent: split multi-orf transcripts if two consecutive ORFs have both BLAST hits
             and none of those hits is against the same target.
- lenient: split multi-orf transcripts as in stringent, and additionally, also when
           either of the ORFs lacks a BLAST hit (but not both).
- permissive: like lenient, but also split when both ORFs lack BLAST hits
- split: split multi-orf transcripts regardless of what BLAST data is available.
If multiple modes are specified, Mikado will create a Daijin-compatible configuration file.""")
    daijin.add_argument("--scheduler", default="", choices=["local", "SLURM", "LSF", "PBS"],
                        help="Scheduler to use. Default: None - ie, either execute everything on the local machine or use DRMAA to submit and control jobs (recommended).")
    daijin.add_argument("--exe", default="daijin_exe.yaml",
                         help="Configuration file for the executables.")
    daijin.add_argument("-c", "--cluster_config",
                         type=str, default=None,
                         help="Cluster configuration file to write to.")
    parser.add_argument("-t", "--threads", default=1, type=int)
    parser.add_argument("--skip-split", dest="skip_split", default=[], nargs="+",
                        help="List of labels for which splitting will be disabled (eg long reads such as PacBio)")
    output_format = parser.add_mutually_exclusive_group()
    output_format.add_argument("-j", "--json", action="store_true", default=False,
                               help="Output will be in JSON (default: inferred by filename, with TOML as fallback).")
    output_format.add_argument("-y", "--yaml", action="store_true", default=False,
                               help="Output will be in YAML (default: inferred by filename, with TOML as fallback).")
    output_format.add_argument("--toml", action="store_true", default=False,
                               help="Output will be in TOML (default: inferred by filename, with TOML as fallback).")
    parser.add_argument("-od", "--out-dir", dest="out_dir", default=None,
                        help="Destination directory for the output.")
    parser.add_argument("out", nargs='?', default=sys.stdout, type=argparse.FileType('w'))
    parser.set_defaults(func=create_config)
    return parser
