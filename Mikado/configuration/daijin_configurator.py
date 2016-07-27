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


def create_daijin_config():

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


def create_config(args):

    logger = create_default_logger("daijin_config")

    config = create_daijin_config()
    assert "reference" in config, config.keys()
    # print(config)
    # config["reference"]["genome"] = args.genome
    # config["reference"]["transcriptome"] = args.transcriptome

    config["name"] = args.name

    if len(args.r1) != len(args.r2):
        exc=InvalidJson(
            """An invalid number of reads has been specified; there are {} left reads and {} right reads.
            Please correct the issue.""".format(len(args.r1), len(args.r2)))
        logger.exception(exc)
        sys.exit(1)

    elif len(args.r1) != len(args.samples):
        exc=InvalidJson(
            """An invalid number of samples has been specified; there are {} left reads and {} samples.
            Please correct the issue.""".format(len(args.r1), len(args.samples)))
        logger.exception(exc)
        sys.exit(1)
    if len(args.strand_specific) == 1 and len(args.r1) > 1:
        logger.warning(
            "Only one strand-specific setting has been specified even if there are multiple samples. I will assume that all the samples have this strand-specificity.")
        args.strand_specific *= len(args.r1)
    elif len(args.strand_specific) == 0:
        logger.warning("No strand specific option specified, so I will assume all the samples are non-strand specific.")
        args.strand_specific = ["fr-unstranded"] * len(args.r1)
    elif len(args.strand_specific) != len(args.r1):
        exc=InvalidJson(
            """An invalid number of strand-specific options has been specified; there are {} left reads and {} samples.
            Please correct the issue.""".format(len(args.r1), len(args.strand_specific)))
        logger.exception(exc)
        sys.exit(1)

    config["short_reads"]["r1"] = args.r1
    config["short_reads"]["r2"] = args.r2
    config["short_reads"]["samples"] = args.samples
    config["short_reads"]["strandedness"] = args.strand_specific


    for method in args.asm_methods:
        config["asm_methods"][method] = [""]
    for method in args.aligners:
        config["aligners"][method] = [""]

    final_config = config.copy()
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
        # sys.exit(1)

    print_config(yaml.dump(final_config, default_flow_style=False), sys.stdout)