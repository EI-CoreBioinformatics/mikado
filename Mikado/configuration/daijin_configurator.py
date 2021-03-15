import dataclasses
from typing import Union
from argparse import Namespace
import rapidjson as json
import os
import toml
import yaml
from pkg_resources import resource_stream
from .configurator import create_cluster_config, load_and_validate_config
from . import print_config
from .daijin_configuration import DaijinConfiguration
from .._transcripts.scoring_configuration import ScoringFile
from ..exceptions import InvalidConfiguration
from ..utilities.log_utils import create_default_logger
import sys
import pysam
try:
    from yaml import CSafeLoader as yLoader
except ImportError:
    from yaml import SafeLoader as yLoader


def _parse_sample_sheet(sample_sheet, config: DaijinConfiguration, logger) -> DaijinConfiguration:

    """Mini-function to parse the sample sheet for Daijin.
    :param sample_sheet: file name of the sample sheet in tabular format
    :type sample_sheet: str
    :param config: configuration object
    :param logger: logger instance
    :rtype DaijinConfiguration
    """

    config.short_reads.r1 = []
    config.short_reads.r2 = []
    config.short_reads.samples = []
    config.short_reads.strandedness = []

    config.long_reads.files = []
    config.long_reads.samples = []
    config.long_reads.strandedness = []

    with open(sample_sheet) as sample_sheet:
        for num, line in enumerate(sample_sheet):
            line = line.rstrip().split("\t")
            if not line:
                # Skip empty lines
                continue
            if len(line) < 3:
                logger.error("Invalid input line", num, "in the sample sheet:", "\t".join(line))
                sys.exit(1)
            r1, r2, sample = line[:3]
            if not r1:
                logger.error("Read 1 undefined for line {}, please correct.".format(num))
            elif not sample:
                logger.error("Sample undefined for line {}, please correct.".format(num))
            strandedness = line[3] if line[3:] else "fr-unstranded"
            if strandedness not in ("fr-unstranded", "fr-secondstrand", "fr-firststrand", "f", "r"):
                logger.error("Invalid strandedness at line {}:".format(num), strandedness)
                sys.exit(1)
            str_to_bool = {"False": False, "True": True}
            is_long_read = str_to_bool.get(line[4] if line[4:] else "False", False)
            if is_long_read and r2:
                logger.error(
                    "I found a long read with mates at line {}, this is not supported. Please double check. Line:\n{}".format(
                        num, "\t".join(line)))
                sys.exit(1)

            if is_long_read:
                config.long_reads.files.append(r1)
                config.long_reads.samples.append(sample)
                config.long_reads.strandedness.append(strandedness)
            else:
                config.short_reads.r1.append(r1)
                config.short_reads.r2.append(r2)
                config.short_reads.samples.append(sample)
                config.short_reads.strandedness.append(strandedness)
    return config


def _parse_reads_from_cli(args, config: DaijinConfiguration, logger) -> DaijinConfiguration:

    """Small function to infer the reads from the CLI.
    :param args: namespace coming from argparse
    :param config: initialised configuration object
    :param logger: logger instance
    """

    if len(args.r1) != len(args.r2):
        exc = InvalidConfiguration(
            """An invalid number of reads has been specified; there are {} left reads and {} right reads.
            Please correct the issue.""".format(len(args.r1), len(args.r2)))
        logger.exception(exc)
        sys.exit(1)
    elif len(args.r1) != len(args.samples):
        exc = InvalidConfiguration(
            """An invalid number of samples has been specified; there are {} left reads and {} samples.
            Please correct the issue.""".format(len(args.r1), len(args.samples)))
        logger.exception(exc)
        sys.exit(1)
    if len(args.strandedness) == 1 and len(args.r1) > 1:
        logger.warning(
            "Only one strand-specific setting has been specified even if there are multiple samples. \
            I will assume that all the samples have this strand-specificity.")
        args.strandedness *= len(args.r1)
    elif len(args.strandedness) == 0:
        logger.warning("No strand specific option specified, so I will assume all the samples are non-strand specific.")
        args.strandedness = ["fr-unstranded"] * len(args.r1)
    elif len(args.strandedness) != len(args.r1):
        exc = InvalidConfiguration(
            """An invalid number of strand-specific options has been specified; there are {} left reads and {} samples.
            Please correct the issue.""".format(len(args.r1), len(args.strandedness)))
        logger.exception(exc)
        sys.exit(1)

    config.short_reads.r1 = args.r1
    config.short_reads.r2 = args.r2
    config.short_reads.samples = args.samples
    config.short_reads.strandedness = args.strandedness
    return config


def create_daijin_config(args: Namespace, config=None, level="ERROR", piped=False) -> Union[DaijinConfiguration, None]:

    """
    Function to create the Daijin configuration object given the CLI arguments.
    :param args: argparse Namespace to use
    :param config: either None or a pre-existing DaijinConfiguration object to update
    :type config: (None|DaijinConfiguration)
    :param level: logging level to use.
    :type level: str
    :param piped: boolean flag. If true, the updated DaijinConfiguration will be returned. Otherwise, the function
    will write the configuration file to the output file specified in the args Namespace.
    :return: (None|DaijinConfiguration)
    """

    logger = create_default_logger("daijin_config", level=level)

    if not isinstance(config, DaijinConfiguration):
        config = DaijinConfiguration()

    if args.seed is not None and isinstance(args.seed, int) and args.seed not in (True, False):
        config.seed = args.seed
    if (args.genome is None or not os.path.exists(args.genome)) and not args.no_files:
        error = "The genome FASTA file {} does not exist! No files: {}".format(args.genome, args.no_files)
        logger.critical(error)
        raise ValueError(error)
    elif not args.no_files:
        config.reference.genome = args.genome
        logger.setLevel("INFO")
        index_present = os.path.exists(config.reference.genome + ".fai")
        if not index_present:
            logger.info("Indexing the genome")
        else:
            logger.debug("Loading the reference index")
        pysam.FastaFile(config.reference.genome)
        if not index_present:
            logger.info("Indexed the genome")
        else:
            logger.debug("Loaded the reference index")

        logger.setLevel(level)
        config.reference.transcriptome = args.transcriptome

    config.name = args.name
    if args.out_dir is None:
        args.out_dir = args.name
    config.out_dir = args.out_dir

    if args.sample_sheet:
        _parse_sample_sheet(args.sample_sheet, config, logger)
    else:
        _parse_reads_from_cli(args, config, logger)

    config.scheduler = args.scheduler
    create_cluster_config(config, args, logger)

    config.threads = args.threads

    config.mikado.modes = args.modes

    failed = False
    if config.short_reads.r1 and not args.aligners:
        logger.critical(
            "No short read aligner selected, but there are short read samples. Please select at least one alignment method.")
        failed = True
    if config.short_reads.r1 and not args.asm_methods:
        logger.critical(
            "No short read assembler selected, but there are short read samples. Please select at least one assembly method.")
        failed = True
    if config.long_reads.files and not args.long_aln_methods:
        logger.critical(
            "No long read aligner selected, but there are long read samples. Please select at least one assembly method.")
        failed = True

    if failed:
        sys.exit(1)

    # Methods that need to be used need to have at least one item in the list of execution parameters to be used.
    # Here we are telling Daijin that it will have to run these programs at least once with default parameters.
    for method in args.aligners:
        setattr(config.align_methods, method, [""])
    for method in args.asm_methods:
        setattr(config.asm_methods, method, [""])
    for method in args.long_aln_methods:
        setattr(config.long_read_align_methods, method, [""])

    # Set and eventually copy the scoring file.
    if args.scoring:
        if args.copy_scoring is not False:
            with open(args.copy_scoring, "wt") as out:
                with resource_stream("Mikado", os.path.join("configuration",
                                                            "scoring_files",
                                                            args.scoring)) as original:
                    for line in original:
                        print(line.decode(), file=out, end="")
            args.scoring = os.path.abspath(args.copy_scoring)
        config.pick.scoring_file = args.scoring
    elif args.new_scoring:
        if os.path.exists(args.new_scoring):
            # Check it's a valid scoring file
            with open(args.new_scoring) as _:
                if args.new_scoring.endswith("json"):
                    new_scoring = json.loads(_.read())
                else:
                    new_scoring = yaml.load(_, Loader=yLoader)

                _ = ScoringFile.Schema().load(new_scoring)
                _.check(minimal_orf_length=config.pick.orf_loading.minimal_orf_length)
        else:
            ns = ScoringFile()
            with open(args.new_scoring, "wt") as out:
                ns = dataclasses.asdict(ns)
                for key in ["as_requirements", "requirements", "not_fragmentary"]:
                    ns[key] = {"parameters": {}, "expression": []}
                if args.new_scoring.endswith("json"):
                    json.dump(ns, out)
                elif args.new_scoring.endswith("yaml"):
                    yaml.dump(ns, out)
                else:
                    toml.dumps(ns, out)
            config.pick.scoring_file = args.new_scoring

    if args.flank is not None:
        config.pick.clustering.flank = args.flank
        config.pick.fragments.max_distance = args.flank
    if args.intron_range is not None:
        args.intron_range = sorted(args.intron_range)
        config.pick.run_options.intron_range = args.intron_range

    config.blastx.prot_db = args.prot_db

    config.mikado.use_diamond = (not args.use_blast)
    config.mikado.use_prodigal = (not args.use_transdecoder)

    final_config = config.copy()

    final_config = load_and_validate_config(final_config)

    if args.exe:
        with open(args.exe, "wt") as out:
            for key, val in dataclasses.asdict(final_config.load).items():
                print("{}: \"{}\"".format(key, val), file=out)

    if piped is True:
        return final_config
    else:
        if args.json is True or (args.out != sys.stdout and args.out.name.endswith("json")):
            format_name = "json"
        elif args.yaml is True or (args.out != sys.stdout and args.out.name.endswith("yaml")):
            format_name = "yaml"
        else:
            format_name = "toml"

        print_config(final_config, args.out, output_format=format_name, no_files=args.no_files, full=args.full)
        args.out.close()
    return
