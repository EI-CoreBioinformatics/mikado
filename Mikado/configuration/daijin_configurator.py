try:
    import ujson as json
except ImportError:
    import json
import os
import io
import yaml
import jsonschema
from pkg_resources import resource_stream, resource_filename
from .configurator import extend_with_default, merge_dictionaries
from ..subprograms.configure import print_config, check_has_requirements
from ..exceptions import InvalidJson
from ..utilities.log_utils import create_default_logger
import sys



def create_daijin_validator():

    resolver = jsonschema.RefResolver("file:///{}".format(
        resource_filename("Mikado.configuration", "configuration")), None)

    validator = extend_with_default(jsonschema.Draft4Validator,
                                    resolver=resolver,
                                    simple=True)

    with io.TextIOWrapper(resource_stream(__name__,
                                          "daijin_schema.json")) as blue:
        blue_print = json.load(blue)

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
        resolver = jsonschema.RefResolver("file:///{}".format(
            resource_filename("Mikado.configuration", "configuration")), None)
        with io.TextIOWrapper(resource_stream(__name__,
                                              "daijin_schema.json")) as blue:
            blue_print = json.load(blue)

        validator = jsonschema.Draft4Validator(blue_print, resolver=resolver)
        validator.validate(config)
    except Exception as exc:
        logger.exception(exc)
        sys.exit(1)


def create_daijin_base_config():

    validator = create_daijin_validator()
    conf = dict()
    validator.validate(conf)

    new_dict = dict()
    composite_keys = [(ckey[1:]) for ckey in
                      check_has_requirements(conf,
                                             validator.schema["properties"])]

    # Sort the composite keys by depth

    for ckey in sorted(composite_keys, key=len, reverse=True):
        assert isinstance(conf, dict), (conf, type(conf))
        defa = conf

        # Get to the latest position
        for key in ckey:
            # if key in ("Comment", "SimpleComment"):
            #     continue
            try:
                defa = defa[key]
            except KeyError:
                raise KeyError(key, defa)
            except TypeError:
                raise TypeError(key, defa)
        val = defa
        for k in reversed(ckey):
            val = {k: val}

        new_dict = merge_dictionaries(new_dict, val)

    return new_dict


def create_daijin_config(args, level="ERROR"):

    logger = create_default_logger("daijin_config", level=level)

    config = create_daijin_base_config()
    assert "reference" in config, config.keys()
    # print(config)
    config["reference"]["genome"] = args.genome
    config["reference"]["transcriptome"] = args.transcriptome

    config["name"] = args.name
    if args.out_dir is None:
        args.out_dir = args.name
    config["out_dir"] = args.out_dir

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
    config["scheduler"] = args.scheduler
    config["threads"] = args.threads

    config["mikado"]["modes"] = args.modes

    for method in args.asm_methods:
        config["asm_methods"][method] = [""]
    for method in args.aligners:
        config["align_methods"][method] = [""]

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
        config["mikado"]["pick"]["scoring_file"] = args.scoring

    if args.flank is not None:
        config["mikado"]["pick"]["clustering"]["flank"] = args.flank

    config["blastx"]["prot_db"] = args.prot_db
    assert "prot_db" in config["blastx"]

    config["mikado"]["use_diamond"] = args.use_diamond

    final_config = config.copy()
    check_config(config, logger)
    assert "prot_db" in config["blastx"]

    if args.cluster_config is not None:
        with open(args.cluster_config, "wb") as out:
            for line in resource_stream("Mikado", os.path.join("daijin", "hpc.yaml")):
                out.write(line)

    if args.out != sys.stdout and args.out.name.endswith("json"):
        json.dump(final_config, args.out)
    else:
        print_config(yaml.dump(final_config, default_flow_style=False), args.out)

    args.out.close()
    return
