import rapidjson as json
import os
import io
import yaml
import tomlkit
import jsonschema
from pkg_resources import resource_stream, resource_filename
from .configurator import extend_with_default, merge_dictionaries, check_all_requirements, check_scoring, to_json
from .configurator import create_cluster_config
from . import print_config, check_has_requirements, print_toml_config
from ..exceptions import InvalidJson
from ..utilities.log_utils import create_default_logger
import sys
import jsonref
import pysam
try:
    from yaml import CSafeLoader as yLoader
except ImportError:
    from yaml import SafeLoader as yLoader


def _substitute_conf(schema):
    """Hack to solve my resolution JSON problems. It will change the $ref to the absolute location
    of this file."""
    base = resource_filename("Mikado", "configuration")

    for key in schema:
        if key == "$ref":
            schema[key] = "file://" + base + "/" + schema[key]
        elif isinstance(schema[key], dict):
            schema[key] = _substitute_conf(schema[key])
        else:
            continue
    return schema


def create_daijin_validator(simple=True):

    cname = resource_filename("Mikado.configuration", "configuration_blueprint.json")

    # We have to repeate twice the ending configuration (bug in jsonref?)
    baseuri = "file://" +  os.path.join(os.path.dirname(cname), os.path.basename(os.path.dirname(cname)))
    with io.TextIOWrapper(resource_stream("Mikado.configuration",
                                          "daijin_schema.json")) as blue:
        blue_print = jsonref.load(blue,
                                  jsonschema=True,
                                  base_uri=baseuri)

    # _substitute_conf(blue_print)
    validator = extend_with_default(jsonschema.Draft7Validator,
                                    simple=simple)
    validator = validator(blue_print)

    return validator


def check_config(config, logger=None):

    """
    Function to check that a configuration abides to the Daijin schema.

    :param config: the dictionary to validate
    :param logger: optional logger. If none is provided, one will be created
    :return:
    """

    if logger is None:
        logger = create_default_logger("daijin_validator")

    try:

        validator = create_daijin_validator()
        # validator = jsonschema.Draft4Validator(blue_print, resolver=resolver)
        validator.validate(config)
    except Exception as exc:
        logger.exception(exc)
        sys.exit(1)


def create_daijin_base_config(simple=True):

    validator = create_daijin_validator(simple=simple)
    conf = dict()
    validator.validate(conf)
    mikado_conf = to_json(None, simple=simple)
    mikado_conf = merge_dictionaries(mikado_conf, conf)

    return mikado_conf


def _parse_sample_sheet(sample_sheet, config, logger):

    """Mini-function to parse the sample sheet."""

    config["short_reads"]["r1"] = []
    config["short_reads"]["r2"] = []
    config["short_reads"]["samples"] = []
    config["short_reads"]["strandedness"] = []
    config["long_reads"]["files"] = []
    config["long_reads"]["samples"] = []
    config["long_reads"]["strandedness"] = []

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
                config["long_reads"]["files"].append(r1)
                config["long_reads"]["samples"].append(sample)
                config["long_reads"]["strandedness"].append(strandedness)
            else:
                config["short_reads"]["r1"].append(r1)
                config["short_reads"]["r2"].append(r2)
                config["short_reads"]["samples"].append(sample)
                config["short_reads"]["strandedness"].append(strandedness)
    return config


def _parse_reads_from_cli(args, config, logger):

    """Small function to infer the reads from the CLI."""

    if len(args.r1) != len(args.r2):
        exc = InvalidJson(
            """An invalid number of reads has been specified; there are {} left reads and {} right reads.
            Please correct the issue.""".format(len(args.r1), len(args.r2)))
        logger.exception(exc)
        sys.exit(1)
    elif len(args.r1) != len(args.samples):
        exc = InvalidJson(
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
        exc = InvalidJson(
            """An invalid number of strand-specific options has been specified; there are {} left reads and {} samples.
            Please correct the issue.""".format(len(args.r1), len(args.strandedness)))
        logger.exception(exc)
        sys.exit(1)

    config["short_reads"]["r1"] = args.r1
    config["short_reads"]["r2"] = args.r2
    config["short_reads"]["samples"] = args.samples
    config["short_reads"]["strandedness"] = args.strandedness
    return config


def create_daijin_config(args, level="ERROR", piped=False):

    logger = create_default_logger("daijin_config", level=level)

    config = create_daijin_base_config(simple=(not args.full))
    if args.seed is not None and isinstance(args.seed, int) and args.seed not in (True, False):
        config["seed"] = args.seed
    assert "reference" in config, config.keys()
    if not os.path.exists(args.genome):
        error = "The genome FASTA file {} does not exist!".format(args.genome)
        logger.critical(error)
        raise ValueError(error)

    config["reference"]["genome"] = args.genome
    logger.setLevel("INFO")
    logger.info("Indexing the genome")
    pysam.FastaFile(config["reference"]["genome"])
    logger.info("Indexed the genome")
    logger.setLevel(level)
    config["reference"]["transcriptome"] = args.transcriptome

    config["name"] = args.name
    if args.out_dir is None:
        args.out_dir = args.name
    config["out_dir"] = args.out_dir

    if args.sample_sheet:
        _parse_sample_sheet(args.sample_sheet, config, logger)
    else:
        _parse_reads_from_cli(args, config, logger)

    config["scheduler"] = args.scheduler
    create_cluster_config(config, args, logger)

    config["threads"] = args.threads

    config["mikado"]["modes"] = args.modes

    failed = False
    if config["short_reads"]["r1"] and not args.aligners:
        logger.critical(
            "No short read aligner selected, but there are short read samples. Please select at least one alignment method.")
        failed = True
    if config["short_reads"]["r1"] and not args.asm_methods:
        logger.critical(
            "No short read assembler selected, but there are short read samples. Please select at least one assembly method.")
        failed = True
    if config.get("long_reads", dict()).get("files", []) and not args.long_aln_methods:
        logger.critical(
            "No long read aligner selected, but there are long read samples. Please select at least one assembly method.")
        failed = True

    if failed:
        sys.exit(1)

    for method in args.aligners:
        config["align_methods"][method] = [""]
    for method in args.asm_methods:
        config["asm_methods"][method] = [""]
    for method in args.long_aln_methods:
        config["long_read_align_methods"][method] = [""]

    # Set and eventually copy the scoring file.
    if args.scoring is not None:
        if args.copy_scoring is not False:
            with open(args.copy_scoring, "wt") as out:
                with resource_stream("Mikado", os.path.join("configuration",
                                                            "scoring_files",
                                                            args.scoring)) as original:
                    for line in original:
                        print(line.decode(), file=out, end="")
            args.scoring = os.path.abspath(args.copy_scoring)
        config.get("pick", dict())["scoring_file"] = args.scoring
    elif args.new_scoring is not None:
        if os.path.exists(args.new_scoring):
            # Check it's a valid scoring file
            with open(args.new_scoring) as _:
                if args.new_scoring.endswith("json"):
                    new_scoring = json.loads(_.read())
                else:
                    new_scoring = yaml.load(_, Loader=yLoader)

                _ = check_all_requirements(new_scoring)
                _ = check_scoring(new_scoring)

        else:
            with io.TextIOWrapper(resource_stream("Mikado.configuration",
                                                  "scoring_blueprint.json")) as schema:
                _ = json.load(schema.read())

            ns = dict()
            with open(args.new_scoring, "wt") as out:
                ns["scoring"] = dict()
                for key in ["as_requirements", "requirements", "not_fragmentary"]:
                    ns["as_requirements"] = {"parameters": {}, "expression": []}
                if args.new_scoring.endswith("json"):
                    json.dump(ns, out)
                else:
                    yaml.dump(ns, out)
            config.get("pick", dict())["scoring_file"] = args.new_scoring

    if args.flank is not None:
        config.get("pick", {}).get("clustering", {})["flank"] = args.flank
        config.get("pick", {}).get("fragments", {})["max_distance"] = args.flank
    if args.intron_range is not None:
        args.intron_range = sorted(args.intron_range)
        config.get("pick", {}).get("run_options", {})["intron_range"] = args.intron_range

    config["blastx"]["prot_db"] = args.prot_db
    assert "prot_db" in config["blastx"]

    config["mikado"]["use_diamond"] = (not args.use_blast)
    config["mikado"]["use_prodigal"] = (not args.use_transdecoder)

    final_config = config.copy()
    check_config(final_config, logger)
    assert "prot_db" in final_config["blastx"]

    if args.exe:
        with open(args.exe, "wt") as out:
            for key, val in final_config["load"].items():
                if "Comment" in key:
                    continue
                else:
                    print("{}: \"{}\"".format(key, val), file=out)

    del final_config["load"]
    final_config.pop("as_requirements", None)
    final_config.pop("scoring", None)
    final_config.pop("requirements", None)
    final_config.pop("not_fragmentary", None)
    if any(key.startswith("_") for key in final_config):
        _to_remove = [_ for _ in final_config if isinstance(_, str) and _.startswith("_")]
        [final_config.pop(key) for key in _to_remove]

    def _rec_delete(d: dict, keep_comments=True):
        _to_remove = []
        if len(d) == 0:
            return d
        elif keep_comments is False and all(
                (isinstance(_, str) and "Comment" in _) for _ in d):
            d = dict()
            return d
        else:
            for key in d:
                if isinstance(d[key], dict):
                    d[key] = _rec_delete(d[key], ("Comment" in key))
                    if len(d[key]) == 0:
                        _to_remove.append(key)
        [d.pop(key) for key in _to_remove]
        return d

    final_config = _rec_delete(final_config)

    if piped is True:
        return final_config
    else:
        if args.json is True:
            json.dump(final_config, args.out)
        elif args.yaml is True:
            print_config(yaml.dump(final_config, default_flow_style=False), args.out)
        elif args.toml is True:
            print_toml_config(tomlkit.dumps(final_config), args.out)
        elif args.out != sys.stdout and args.out.name.endswith("json"):
            json.dump(final_config, args.out)
        elif args.out != sys.stdout and args.out.name.endswith("yaml"):
            print_config(yaml.dump(final_config, default_flow_style=False), args.out)
        else:
            print_toml_config(tomlkit.dumps(final_config), args.out)

        args.out.close()
    return
