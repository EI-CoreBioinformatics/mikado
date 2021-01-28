#!/usr/bin/env python3
# coding: utf-8

"""
This module defines the functionalities needed to verify the integrity and completeness
of Mikado configuration files. Missing values are replaced with default ones,
while existing values are checked for type and consistency.
"""
import dataclasses

import io
import json
import os.path
import re
import marshmallow
from multiprocessing import get_start_method
from logging import Logger
import jsonschema
import yaml
from pkg_resources import resource_stream, resource_filename
from ..exceptions import InvalidJson, UnrecognizedRescaler
from ..utilities.log_utils import create_default_logger
import random
import toml
from .configuration import MikadoConfiguration
from .daijin_configuration import DaijinConfiguration
from typing import Union
try:
    from yaml import CSafeLoader as yLoader
except ImportError:
    from yaml import SafeLoader as yLoader
from argparse import Namespace


__author__ = "Luca Venturini"


def create_cluster_config(config: Union[MikadoConfiguration, DaijinConfiguration], args: Namespace, logger: Logger):

    """
    Function to create a configuration for the scheduler for snakemake. It will contain information regarding
    queues, memory, threads etc. per program call.
    :param config: Configuration object
    :param args: namespace with the parameters to use for this run
    :param logger: logger for this run
    :return:
    """

    if (config.scheduler and config.scheduler != "local") or (not config.scheduler and args.cluster_config):
        if not args.queue:
            error = "A queue must be specified for the analysis when in HPC mode. Please relaunch."
            logger.error(error)
            exit(1)
        if args.cluster_config is not None:
            cluster_config = args.cluster_config
        else:
            cluster_config = "daijin_hpc.yaml"
        with open(cluster_config, "wt") as out, \
                resource_stream("Mikado", os.path.join("daijin", "hpc.yaml")) as original:
            for pos, line in enumerate(original):
                print(line.decode(), file=out, end="")
                if pos == 0:
                    print("    queue:", args.queue, file=out)


def extend_with_default(validator_class: jsonschema.Draft7Validator, resolver=None, simple=False):
    """
    Function to extend the normal validation classes for jsonschema
    so that they also set the default values provided inside the schema
    itself. Source:
    https://python-jsonschema.readthedocs.org/en/latest/faq/?highlight=default
    :param validator_class: the validator class to extend (e.g. Draft7Validator)

    :param simple: boolean flag. If set to True, only required properties will be extended.
    :type simple: bool

    :return:
    """
    validate_properties = validator_class.VALIDATORS["properties"]

    def set_default(instance, properties, simple_comment=False):
        """
        Recursive function that sets the default parameters inside "object"
        types for the dictionary instance. It also loads comments, if available.
        :param instance: a dictionary instance
        :type instance: dict

        :param properties: the "properties" dictionary inside the schema
        :type properties: dict

        :param simple_comment: a boolean flag to indicate whether
                               we want the simplified schema or not.
        :type simple_comment: bool

        :return: the prepared instance
        :rtype: dict
        """
        if not isinstance(instance, dict):
            return

        if simple_comment is True:
            if "required" in properties:
                required = properties["required"]
            else:
                required = []

        if "properties" in properties:
            properties = properties["properties"]

        to_remove = []

        for prop, subschema in properties.items():
            if instance is None:
                instance = dict()
            if "default" in subschema and ((not simple) or (simple and prop in required)):
                instance.setdefault(prop, subschema["default"])
            elif prop not in instance:
                if "type" not in subschema:
                    continue
                elif subschema["type"] == "object":
                    instance[prop] = dict()
                    instance[prop] = set_default(instance[prop],
                                                 subschema,
                                                 simple_comment=simple_comment)
                elif simple_comment is True and prop not in required:
                    to_remove.append(prop)
                    continue

        if to_remove:
            [instance.pop(prop, None) for prop in to_remove]

        return instance

    def set_defaults(validator, properties, instance, schema):
        """ Function that sets the default values from the schema
        into the dictionary.
        :param validator: the validator class
        :param properties: schema["properties"] if any
        :param instance: the dictionary loaded from the external file
        :param schema: the schema to be used as blueprint.
        :return:
        """
        for error in validate_properties(validator, properties, instance, schema):
            yield error
        # noinspection PyUnusedLocal
        instance = set_default(instance, properties, simple_comment=simple)

    return jsonschema.validators.extend(
        validator_class, {"properties": set_defaults, "resolver": resolver},
    )


def check_scoring(json_conf: Union[MikadoConfiguration, DaijinConfiguration]):
    """
    :param json_conf: configuration object to check.
    :type json_conf: (MikadoConfiguration|DaijinConfiguration)

    Function to check the "scoring" section of the configuration.
    Each scoring function will be checked for:
    - validity of the expression (it can be interpreted by Mikado)
    - validity of the parameter (it is a valid Metric)

    :return: json_conf
    :rtype: (MikadoConfiguration|DaijinConfiguration)

    """

    with io.TextIOWrapper(resource_stream("Mikado.configuration",
                                          "scoring_blueprint.json")) as schema:
        scoring_schema = json.loads(schema.read())

    parameters_found = set()
    parameters_not_found = []
    double_parameters = []
    invalid_raw = set()
    invalid_filter = set()
    from ..transcripts.transcript import Transcript
    available_metrics = Transcript.get_available_metrics()
    scoring = json_conf.scoring

    if scoring is None:
        raise InvalidJson("\"Scoring\" key not find in the scoring dictionary.")
    elif len(scoring) == 0:
        raise InvalidJson("Empty \"scoring\" dictionary.")

    jdict = scoring.copy()
    validator = extend_with_default(jsonschema.Draft7Validator)

    for parameter in jdict:
        if parameter not in available_metrics:
            # Leniency for external_scores
            if "." in parameter and len(parameter.split(".")) == 2 and (parameter.split(".")[0] == "external" or
                                                                        parameter.split(".")[0] == "attributes"):
                pass
            else:
                parameters_not_found.append(parameter)

        if parameter in parameters_found:
            double_parameters.append(parameter)

        if not jsonschema.Draft7Validator(scoring_schema).is_valid(
                jdict[parameter]):
            errors = [str(_) for _ in list(jsonschema.Draft7Validator(scoring_schema).iter_errors(
                jdict[parameter]))]
            raise InvalidJson("Invalid scoring for {}:\n{}".format(
                parameter, "\n".join(errors)))

        try:
            validator(scoring_schema).validate(jdict[parameter])
        except Exception as err:
            raise ValueError(parameter, err)

        if (not parameter.startswith("external") and not parameter.startswith("attributes") and
                    getattr(Transcript, parameter).usable_raw is False and
                    jdict[parameter]["use_raw"] is True):
            invalid_raw.add(parameter)

        if jdict[parameter]["rescaling"] == "target":
            if "value" not in jdict[parameter]:
                message = """Target rescaling requested for {0} but no target value specified.
                    Please specify it with the \"value\" keyword.\n{1}"""
                message = message.format(parameter, json_conf.scoring[parameter])
                raise UnrecognizedRescaler(message)
            elif jdict[parameter]["use_raw"] is True:
                invalid_raw.add(parameter)

        if "filter" in jdict[parameter] and "metric" in jdict[parameter]["filter"]:
            if jdict[parameter]["filter"]["metric"] not in available_metrics:
                parameters_not_found.append(jdict[parameter]["filter"]["metric"])

    if len(parameters_not_found) > 0 or len(double_parameters) > 0 or len(invalid_filter) > 0 or len(invalid_raw) > 0:
        err_message = ''
        if len(parameters_not_found) > 0:
            err_message = """The following parameters, present in the JSON file,
            are not available!\n\t{0}\n""".format(
                "\n\t".join(parameters_not_found))
        if len(double_parameters) > 0:
            err_message += """The following parameters have been specified more than once,
            please correct:\n\t{0}""".format("\n\t".join(list(double_parameters)))
        if len(invalid_filter) > 0:
            err_message += """The following parameters have an invalid filter,
            please correct:
            \t{0}""".format("\n\t".join(list(invalid_filter)))
        if len(invalid_raw) > 0:
            err_message += """The following parameters cannot be used as raw scores, either
            because they are not normalized on their own, or because a "target" rescaling has been asked for:
            \t{0}""".format("\n\t".join(list(invalid_raw)))
        raise InvalidJson(err_message)

    json_conf.scoring = jdict
    assert json_conf.scoring is not None

    return json_conf


def check_all_requirements(configuration: Union[MikadoConfiguration,DaijinConfiguration]):
    """
    Function to check all the "requirements" sections of the configuration.

    :param configuration: the configuration object to be checked.
    :type configuration: (MikadoConfiguration|DaijinConfiguration)

    :return: configuration
    :rtype (MikadoConfiguration|DaijinConfiguration)

    """

    with io.TextIOWrapper(resource_stream("Mikado.configuration",
                                          "requirements_blueprint.json")) as rs_blueprint:
        require_schema = json.loads(rs_blueprint.read())

    if not hasattr(configuration, "requirements") or len(configuration.requirements) == 0:
        raise InvalidJson("No minimal requirements specified!")

    for section_name in ["requirements", "as_requirements", "cds_requirements", "not_fragmentary"]:
        if not hasattr(configuration, section_name) or getattr(configuration, section_name) is None:
            section = dict()
            if section_name == "cds_requirements":
                section["parameters"] = {
                    "selected_cds_length": {"operator": "ge",
                                            "value": configuration.pick.orf_loading.minimal_orf_length}}
                section["expression"] = ["selected_cds_length"]
            else:
                section.update(configuration.requirements)
                section["expression"] = section["__expression"][:]
                del section["__expression"]
        else:
            section = getattr(configuration, section_name)

        # Check requirements will MODIFY IN PLACE the expression, so the copying
        # must happen before, not after.

        try:
            section = check_requirements(section, require_schema=require_schema, index=section_name)
        except (InvalidJson, SyntaxError) as exc:
            print(section["expression"])
            print(type(section["expression"]))
            raise exc
        setattr(configuration, section_name, section)

    return configuration


key_pattern = re.compile(r"([^ ()]+)")


def check_requirements(section: dict, require_schema: jsonschema.Draft7Validator, index: str):
    """
    Function to check the "requirements" section of the configuration.
    Each filtering function will be checked for:
    - validity of the expression (it can be interpreted by Mikado)
    - validity of the parameter (it is a valid Metric)

    :param section: configuration dictionary to check.
    :type section: dict

    :param require_schema: the requirements section of the JSON schema.
    :type require_schema: dict

    :param index: the key we are inspecting in the master dictionary
    :type index: str

    :return: section
    :rtype dict

    """

    # Check that the parameters are valid
    parameters_not_found = []
    from ..transcripts.transcript import Transcript
    available_metrics = Transcript.get_available_metrics()

    if "parameters" not in section:
        raise InvalidJson(
            "The {} field must have a \"parameters\" subfield!".format(index))
    for key in section["parameters"]:
        key_value = None
        dots = key.split(".")
        if dots[0] == "external" or dots[0] == "attributes":
            if len(dots) == 1 or len(dots) > 3:
                parameters_not_found.append(dots[0])
                continue
            elif len(dots) == 2:
                key_name = ".".join(dots)
            else:
                key_name = ".".join(dots[:-1])
                # print(key_name)
            key_value = dots[1]
        else:
            key_name = dots[0]
        if key_name not in available_metrics:
            if key_value is not None:
                pass
            else:
                parameters_not_found.append(key_name)
                continue
        if not jsonschema.Draft7Validator(require_schema["definitions"]["parameter"]).is_valid(
                section["parameters"][key]):
            errors = list(jsonschema.Draft7Validator(require_schema).iter_errors(
                section["parameters"][key]
            ))
            raise InvalidJson("Invalid parameter for {0} in {1}: \n{2}".format(
                key, index, errors))

        section["parameters"][key]["name"] = key_name
        if key_value:
            section["parameters"][key]["source"] = key_value

    if len(parameters_not_found) > 0:
        raise InvalidJson(
            "The following parameters, selected for filtering, are invalid:\n\t{0}".format(
                "\n\t".join(parameters_not_found)
            ))

    # Create automatically a filtering expression
    if "__expression" in section:
        section["expression"] = section["__expression"][:]

    if "expression" not in section:
        section["expression"] = " and ".join(
            list(section["parameters"].keys()))
        keys = section["parameters"].keys()
        newexpr = section["expression"][:]
        section["__expression"] = section["expression"][:]
    else:
        if not jsonschema.Draft7Validator(
                require_schema["definitions"]["expression"]).is_valid(
                    section["expression"]):
            raise InvalidJson("Invalid expression field")

        section["__expression"] = section["expression"][:]

        expr = " ".join(section["expression"])
        newexpr = expr[:]

        keys = set([key for key in key_pattern.findall(expr) if key not in ("and", "or", "not", "xor")])

        diff_params = set.difference(
            set(keys), set(section["parameters"].keys()))

        if len(diff_params) > 0:
            raise InvalidJson(
                "Expression and required parameters mismatch:\n\t{0}".format(
                    "\n\t".join(list(diff_params))))

    for key in keys:  # Create the final expression
        newexpr = re.sub(r"\b{}\b".format(key), "evaluated[\"{0}\"]".format(key), newexpr)

    # Test the expression
    try:
        compile(newexpr, "<json>", "eval")
    except SyntaxError:
        raise InvalidJson("Invalid expression for {}:\n{}".format(index, newexpr))

    section["expression"] = newexpr
    return section


def check_db(configuration: Union[MikadoConfiguration, DaijinConfiguration]) -> Union[MikadoConfiguration,
                                                                                      DaijinConfiguration]:

    """
    Function to check the validity of the database options.
    :param configuration: (MikadoConfiguration|DaijinConfiguration)
    :return: 
    """

    if configuration.db_settings.dbtype in ("mysql", "postgresql"):
        if configuration.db_settings.dbhost is None:
            raise InvalidJson(
                "No host specified for the {0} database!".format(
                    configuration.db_settings.dbtype))
        if configuration.db_settings.dbuser is None:
            raise InvalidJson(
                "No user specified for the {0} database!".format(
                    configuration.db_settings.dbtype))
        if configuration.db_settings.dbport == 0:
            if configuration.db_settings.dbtype == "mysql":
                configuration.db_settings.dbport = 3306
            else:
                configuration.db_settings.dbport = 5432
    else:
        if (configuration.db_settings.db is not None and
                not os.path.isabs(configuration.db_settings.db)):
            if os.path.exists(os.path.join(os.getcwd(),
                                           configuration.db_settings.db)):
                configuration.db_settings.db = os.path.join(
                    os.getcwd(), configuration.db_settings.db)
            elif os.path.exists(os.path.join(os.path.dirname(configuration.filename or ""), configuration.db_settings.db)):
                configuration.db_settings.db = os.path.join(os.path.dirname(configuration.filename or ""),
                                                        configuration.db_settings.db)
            else:
                configuration.db_settings.db = os.path.join(os.path.dirname(configuration.filename or ""),
                                                        configuration.db_settings.db)

    return configuration


def check_scoring_file(configuration: Union[MikadoConfiguration, DaijinConfiguration], logger: Logger) -> (
        Union[DaijinConfiguration, MikadoConfiguration], bool):

    """The purpose of this section is the following:
    - check that the scoring file exists somewhere different from the system folder. If it does, check whether it is
      a valid file.
    - If the local file is broken / wrong / inexistent, use the system one.
    - load the scoring file into the configuration object.
    
    :param configuration: configuration object to check
    :param logger: logger to use
    """

    overwritten = False
    if configuration.pick.scoring_file is None:
        if getattr(configuration, "_loaded_scoring", configuration.pick.scoring_file) != configuration.pick.scoring_file:
            logger.debug("Resetting the scoring to its previous value")
            configuration.pick.scoring_file = getattr(configuration, "_loaded_scoring", configuration.pick.scoring_file)

    elif getattr(configuration, "_loaded_scoring", configuration.pick.scoring_file) != configuration.pick.scoring_file:

            logger.debug("Overwriting the scoring configuration using '%s' as scoring file",
                         configuration.pick.scoring_file)
            [setattr(configuration, _, None) for _ in ("_loaded_scoring",
                "scoring", "requirements", "as_requirements", "not_fragmentary")]

    # FIXME: The scoring needs to be handled either within the Configuration object or made into a separate one
    elif all(getattr(configuration, _, None) is not None for _
             in ["scoring", "requirements", "as_requirements", "not_fragmentary"]):
        try:
            check_scoring(configuration)
            check_all_requirements(configuration)
            configuration._loaded_scoring = configuration.pick.scoring_file
            logger.debug("Verified everything is OK for the scoring, returning")
            return configuration, overwritten
        except InvalidJson:
            logger.debug("Invalid scoring for the jconf, resetting")
            [setattr(configuration, _, None) for _ in ("scoring", "requirements", "as_requirements", "not_fragmentary")]
            configuration._loaded_scoring = None
    else:
        logger.debug("Restarting")
        configuration._loaded_scoring = None

    overwritten = True

    options = [os.path.abspath(configuration.pick.scoring_file),
               os.path.abspath(os.path.join(os.path.dirname(configuration.filename or ""),
                                            configuration.pick.scoring_file)),
               os.path.abspath(os.path.join(resource_filename("Mikado.configuration", "scoring_files"),
                               configuration.pick.scoring_file))]

    found = False

    for option in options:
        if not os.path.exists(option):
            continue
        if option.endswith(("yaml", "json", "toml")):
            with open(option) as scoring_file:
                if option.endswith("yaml"):
                    scoring = yaml.load(scoring_file, Loader=yLoader)
                elif option.endswith("toml"):
                    scoring = toml.load(scoring_file)
                else:
                    scoring = json.loads(scoring_file.read())
                if not isinstance(scoring, dict):
                    continue
                for key, item in scoring.items():
                    assert hasattr(configuration, key), key
                    setattr(configuration, key, item)
                try:
                    configuration = check_scoring(configuration)
                    configuration = check_all_requirements(configuration)
                except InvalidJson:
                    continue
                for key in scoring:
                    setattr(configuration, key, scoring[key])
                found = True
                configuration.scoring_file = option
        if found is True:
            logger.info("Found the correct option: %s", option)
            configuration.scoring_file = option
            break

    return configuration, overwritten


def check_and_load_scoring(configuration: Union[DaijinConfiguration, MikadoConfiguration],
                           logger=None) -> Union[DaijinConfiguration, MikadoConfiguration]:

    """
    Wrapper for the various checks performed on the configuration file.
    :param configuration: the configuration object to finalise
    :type configuration: (MikadoConfiguration|DaijinConfiguration)
    :param logger: external logger instance
    :type logger: (None|Logger)
    :return configuration
    :rtype: (MikadoConfiguration|DaijinConfiguration)
    """

    if not isinstance(logger, Logger):
        logger = create_default_logger("check_json")

    try:

        configuration, overwritten = check_scoring_file(configuration, logger)
        configuration = check_db(configuration)
        if not configuration.multiprocessing_method:
            configuration.multiprocessing_method = get_start_method()
    except Exception as exc:
        logger.exception(exc)
        raise

    if overwritten is True and configuration.scoring is not None:
        logger.debug("Scoring parameters: {}".format("\n".join(["\n"] + [
            "{}: {}".format(_, configuration.scoring[_]) for _ in configuration.scoring.keys()])))

    seed = configuration.seed
    if seed is None:
        # seed = numpy.random.randint(0, 2**32 - 1)
        seed = random.randint(0, 2**32 - 1)
        logger.info("Random seed: {}", seed)
        configuration.seed = seed

    if seed is not None:
        # numpy.random.seed(seed % (2 ** 32 - 1))
        random.seed(seed % (2 ** 32 - 1))
    else:
        # numpy.random.seed(None)
        random.seed(None)

    return configuration


def load_and_validate_config(raw_configuration: Union[None, MikadoConfiguration, DaijinConfiguration,
                             str, dict], logger=None) -> Union[MikadoConfiguration, DaijinConfiguration]:
    """
    Function to serialise the JSON for configuration and check its consistency.

    :param raw_configuration: either the file name of the configuration or an initialised object to check and finalise.
    :type raw_configuration: (str | None | dict | MikadoConfiguration | DaijinConfiguration)

    :param logger: optional logger to be used.
    :type logger: Logger

    :rtype: (MikadoConfiguration|DaijinConfiguration)
    """

    if not isinstance(logger, Logger):
        logger = create_default_logger("to_json")

    try:
        if isinstance(raw_configuration, (MikadoConfiguration, DaijinConfiguration)):
            config = raw_configuration
        elif isinstance(raw_configuration, dict):
            try:
                config = MikadoConfiguration.Schema().load(raw_configuration)
            except marshmallow.exceptions.ValidationError:
                config = DaijinConfiguration.Schema().load(raw_configuration)
        elif raw_configuration is None or raw_configuration == '' or raw_configuration == dict():
            config = MikadoConfiguration()
        else:
            assert isinstance(raw_configuration, str), raw_configuration
            raw_configuration = os.path.abspath(raw_configuration)
            if not os.path.exists(raw_configuration) or os.stat(raw_configuration).st_size == 0:
                raise InvalidJson("JSON file {} not found!".format(raw_configuration))
            with open(raw_configuration) as json_file:
                # YAML *might* try to load up the file as the proper object
                if raw_configuration.endswith(".yaml"):
                    config = yaml.load(json_file, Loader=yaml.Loader)
                elif raw_configuration.endswith(".toml"):
                    config = toml.load(json_file)
                else:
                    config = json.loads(json_file.read())
            assert isinstance(config, (dict, MikadoConfiguration, DaijinConfiguration)), (config, type(config))
            if isinstance(config, dict):
                try:
                    config = MikadoConfiguration.Schema().load(config)
                except:
                    config = DaijinConfiguration.Schema().load(config)
            config.filename = raw_configuration

        assert isinstance(config, (MikadoConfiguration, DaijinConfiguration)), type(config)
        config = check_and_load_scoring(config, logger=logger)
    except Exception as exc:
        logger.exception(exc)
        raise InvalidJson("The configuration file passed is invalid. Please double check.")

    if config.seed == 0 or config.seed is None:
        config.seed = random.randint(1, 2 ** 32 - 1)
        logger.info("Random seed: {}", config.seed)

    random.seed(config.seed % (2 ** 32 - 1))

    return config
