#!/usr/bin/env python3
# coding: utf-8

"""
This module defines the functionalities needed to verify the integrity and completeness
of Mikado configuration files. Missing values are replaced with default ones,
while existing values are checked for type and consistency.
"""

import io
import json
import os.path
import pickle
import re
from multiprocessing import get_start_method
from logging import Logger
import jsonschema
import pkg_resources
import yaml
from pkg_resources import resource_stream, resource_filename
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from ..transcripts.transcript import Transcript
from ..exceptions import InvalidJson, UnrecognizedRescaler
from ..utilities import merge_dictionaries
from ..utilities.log_utils import create_default_logger

# from frozendict import frozendict

__author__ = "Luca Venturini"


def extend_with_default(validator_class, resolver=None, simple=False):
    """
    Function to extend the normal validation classes for jsonschema
    so that they also set the default values provided inside the schema
    itself. Source:
    https://python-jsonschema.readthedocs.org/en/latest/faq/?highlight=default
    :param validator_class: the validator class to extend (e.g. Draft4Validator)

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
        for prop, subschema in properties.items():
            if instance is None:
                instance = dict()
            if "$ref" in subschema:
                # Automatically resolve and load the reference
                assert resolver is not None
                properties[prop] = resolver.resolve(subschema["$ref"])[1]
                # subschema = resolver.resolve(subschema["$ref"])[1]
                subschema = properties[prop]
            if "default" in subschema:
                instance.setdefault(prop, subschema["default"])
            elif prop not in instance:
                if "type" not in subschema:
                    continue
                elif subschema["type"] == "object":
                    instance[prop] = dict()
                    if not simple_comment and "Comment" in subschema:
                        instance[prop].setdefault("Comment", subschema["Comment"])
                    elif simple_comment and "SimpleComment" in subschema:
                        instance[prop].setdefault("SimpleComment",
                                                  subschema["SimpleComment"])
                    instance[prop] = set_default(instance[prop],
                                                 subschema["properties"],
                                                 simple_comment=simple_comment)

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
        validator_class, {"properties": set_defaults},
    )


def check_scoring(json_conf):
    """
    :param json_conf: configuration dictionary to check.
    :type json_conf: dict

    Function to check the "scoring" section of the configuration.
    Each scoring function will be checked for:
    - validity of the expression (it can be interpreted by Mikado)
    - validity of the parameter (it is a valid Metric)

    :return: json_conf
    :rtype dict
    """

    with io.TextIOWrapper(resource_stream(__name__,
                                          "scoring_blueprint.json")) as schema:
        scoring_schema = json.load(schema)

    parameters_found = set()
    parameters_not_found = []
    double_parameters = []
    invalid_raw = set()
    invalid_filter = set()
    available_metrics = Transcript.get_available_metrics()
    if "scoring" not in json_conf or len(json_conf["scoring"].keys()) == 0:
        raise InvalidJson("No parameters specified for scoring!")

    validator = extend_with_default(jsonschema.Draft4Validator)

    for parameter in json_conf["scoring"]:
        if parameter not in available_metrics:
            # Leniency for external_scores
            if "." in parameter and len(parameter.split(".")) == 2 and parameter.split(".")[0] == "external":
                pass
            else:
                parameters_not_found.append(parameter)

        if parameter in parameters_found:
            double_parameters.append(parameter)

        if not jsonschema.Draft4Validator(scoring_schema).is_valid(
                json_conf["scoring"][parameter]):
            errors = [str(_) for _ in list(jsonschema.Draft4Validator(scoring_schema).iter_errors(
                json_conf["scoring"][parameter]))]
            raise InvalidJson("Invalid scoring for {}:\n{}".format(
                parameter, "\n".join(errors)))

        try:
            validator(scoring_schema).validate(json_conf["scoring"][parameter])
        except Exception as err:
            raise ValueError(parameter, err)

        if ("external" not in parameter and
                    getattr(Transcript, parameter).usable_raw is False and
                    json_conf["scoring"][parameter]["use_raw"] is True):
            invalid_raw.add(parameter)

        if json_conf["scoring"][parameter]["rescaling"] == "target":
            if "value" not in json_conf["scoring"][parameter]:
                message = """Target rescaling requested for {0} but no target value specified.
                    Please specify it with the \"value\" keyword.\n{1}"""
                message = message.format(parameter, json_conf["scoring"][parameter])
                raise UnrecognizedRescaler(message)
            elif json_conf["scoring"][parameter]["use_raw"] is True:
                invalid_raw.add(parameter)

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

    return json_conf


def check_all_requirements(json_conf):
    """
    Function to check all the "requirements" sections of the configuration.

    :param json_conf: the dictionary to be checked.
    :type json_conf: dict

    :return: json_conf
    :rtype dict
    """

    with io.TextIOWrapper(resource_stream(__name__,
                                          "requirements_blueprint.json")) as rs_blueprint:
        require_schema = json.load(rs_blueprint)

    if "requirements" not in json_conf:
        raise InvalidJson("No minimal requirements specified!")

    if "not_fragmentary" not in json_conf:
        json_conf["not_fragmentary"] = dict()
        json_conf["not_fragmentary"].update(json_conf["requirements"])
        # print("Original:", json_conf["not_fragmentary"]["expression"])
        # print("Original:", type(json_conf["not_fragmentary"]["expression"]))

    try:
        # print("Original:", json_conf["not_fragmentary"]["expression"])
        # print("Original:", type(json_conf["not_fragmentary"]["expression"]))
        json_conf = check_requirements(json_conf,
                                       require_schema,
                                       "not_fragmentary")
    except (InvalidJson, SyntaxError) as exc:
        print(json_conf["not_fragmentary"]["expression"])
        print(type(json_conf["not_fragmentary"]["expression"]))
        raise exc

    if "as_requirements" not in json_conf:
        json_conf["as_requirements"] = dict()
        json_conf["as_requirements"].update(json_conf["requirements"])

    # Check requirements will MODIFY IN PLACE the expression, so the copying
    # must happend before, not after.

    try:
        json_conf = check_requirements(json_conf,
                                       require_schema,
                                       "as_requirements")
    except (InvalidJson, SyntaxError) as exc:
        print(json_conf["as_requirements"]["expression"])
        print(type(json_conf["as_requirements"]["expression"]))
        raise exc

    try:
        json_conf = check_requirements(json_conf,
                                       require_schema,
                                       "requirements")
    except (InvalidJson, SyntaxError) as exc:
        print(json_conf["as_requirements"]["expression"])
        print(type(json_conf["as_requirements"]["expression"]))
        raise exc

    return json_conf


def check_requirements(json_conf, require_schema, index):
    """
    Function to check the "requirements" section of the configuration.
    Each filtering function will be checked for:
    - validity of the expression (it can be interpreted by Mikado)
    - validity of the parameter (it is a valid Metric)

    :param json_conf: configuration dictionary to check.
    :type json_conf: dict

    :param require_schema: the requirements section of the JSON schema.
    :type require_schema: dict

    :return: json_conf
    :rtype dict
    """

    # Check that the parameters are valid
    parameters_not_found = []
    available_metrics = Transcript.get_available_metrics()

    if "parameters" not in json_conf[index]:
        raise InvalidJson(
            "The {} field must have a \"parameters\" subfield!".format(index))
    for key in json_conf[index]["parameters"]:
        key_value = None
        dots = key.split(".")
        if dots[0] == "external":
            if len(dots) == 1 or len(dots) > 3:
                parameters_not_found.append(dots[0])
                continue
            elif len(dots) == 2:
                key_name = ".".join(dots)
            else:
                key_name = ".".join(dots[:-1])
                print(key_name)
            key_value = dots[1]
        else:
            key_name = dots[0]
        if key_name not in available_metrics:
            if key_value is not None:
                pass
            else:
                parameters_not_found.append(key_name)
                continue
        if not jsonschema.Draft4Validator(require_schema["definitions"]["parameter"]).is_valid(
                json_conf[index]["parameters"][key]):
            errors = list(jsonschema.Draft4Validator(require_schema).iter_errors(
                json_conf[index]["parameters"][key]
            ))
            raise InvalidJson("Invalid parameter for {0} in {1}: \n{2}".format(
                key, index, errors))

        json_conf[index]["parameters"][key]["name"] = key_name
        if key_value:
            json_conf[index]["parameters"][key]["source"] = key_value

    if len(parameters_not_found) > 0:
        raise InvalidJson(
            "The following parameters, selected for filtering, are invalid:\n\t{0}".format(
                "\n\t".join(parameters_not_found)
            ))

    # Create automatically a filtering expression
    if "__expression" in json_conf[index]:
        json_conf[index]["expression"] = json_conf[index]["__expression"][:]

    if "expression" not in json_conf[index]:
        json_conf[index]["expression"] = " and ".join(
            list(json_conf[index]["parameters"].keys()))
        keys = json_conf[index]["parameters"].keys()
        newexpr = json_conf[index]["expression"][:]
        json_conf[index]["__expression"] = json_conf[index]["expression"][:]
    else:
        if not jsonschema.Draft4Validator(
                require_schema["definitions"]["expression"]).is_valid(
                    json_conf[index]["expression"]):
            raise InvalidJson("Invalid expression field")

        json_conf[index]["__expression"] = json_conf[index]["expression"][:]

        expr = " ".join(json_conf[index]["expression"])
        newexpr = expr[:]

        keys = set([key for key in re.findall(
            "([^ ()]+)", expr) if key not in ("and", "or", "not", "xor")])

        diff_params = set.difference(
            set(keys), set(json_conf[index]["parameters"].keys()))

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

    json_conf[index]["expression"] = newexpr
    return json_conf


def check_db(json_conf):

    """
    Function to check the validity of the database options.
    :param json_conf:
    :return:
    """

    if json_conf["db_settings"]["dbtype"] in ("mysql", "postgresql"):
        if "dbhost" not in json_conf["db_settings"]:
            raise InvalidJson(
                "No host specified for the {0} database!".format(
                    json_conf["db_settings"]["dbtype"]))
        if "dbuser" not in json_conf["db_settings"]:
            raise InvalidJson(
                "No user specified for the {0} database!".format(
                    json_conf["db_settings"]["dbtype"]))
        if json_conf["db_settings"]["dbport"] == 0:
            if json_conf["db_settings"]["dbtype"] == "mysql":
                json_conf["db_settings"]["dbport"] = 3306
            else:
                json_conf["db_settings"]["dbport"] = 5432
    else:
        if (json_conf["db_settings"]["db"] is not None and
                not os.path.isabs(json_conf["db_settings"]["db"])):
            if os.path.exists(os.path.join(os.getcwd(),
                                           json_conf["db_settings"]["db"])):
                json_conf["db_settings"]["db"] = os.path.join(
                    os.getcwd(), json_conf["db_settings"]["db"])
            elif os.path.exists(os.path.join(
                os.path.dirname(json_conf["filename"]),
                    json_conf["db_settings"]["db"]
            )):
                json_conf["db_settings"]["db"] = os.path.join(
                    os.path.dirname(json_conf["filename"]),
                    json_conf["db_settings"]["db"])
            else:
                json_conf["db_settings"]["db"] = os.path.join(
                    os.path.dirname(json_conf["filename"]),
                    json_conf["db_settings"]["db"])

    return json_conf


def create_validator(simple=False):

    """Method to create a validator class (see extend_with_default).
    The simple keyword (boolean) is used to determine whether to keep
    only SimpleComment or full Comments from the schema.

    :type simple: bool

    :return validator
    :rtype: jsonschema.Draft4Validator
    """

    validator = extend_with_default(jsonschema.Draft4Validator,
                                    simple=simple)

    resolver = jsonschema.RefResolver("file:///{}".format(os.path.abspath(
        os.path.dirname(pkg_resources.resource_filename(__name__, __file__))
    )), None)

    with io.TextIOWrapper(resource_stream(__name__,
                                          "configuration_blueprint.json")) as blue:
        blue_print = json.load(blue)

    validator = validator(blue_print, resolver=resolver)

    return validator


def _check_scoring_file(json_conf, logger):

    overwritten = False

    if json_conf.get("__loaded_scoring", json_conf["pick"]["scoring_file"]) != json_conf["pick"]["scoring_file"]:
        logger.info("Overwriting the scoring configuration using '%s' as scoring file",
                    json_conf["pick"]["scoring_file"])
        overwritten = True
        [json_conf.pop(_, None) for _ in ("scoring", "requirements", "as_requirements", "not_fragmentary")]

    if os.path.exists(os.path.abspath(json_conf["pick"]["scoring_file"])):
        json_conf["pick"]["scoring_file"] = os.path.abspath(
            json_conf["pick"]["scoring_file"])
    elif os.path.exists(os.path.join(
            os.path.dirname(json_conf["filename"]),
            json_conf["pick"]["scoring_file"])):
        json_conf["pick"]["scoring_file"] = os.path.join(
            os.path.dirname(json_conf["filename"]),
            json_conf["pick"]["scoring_file"])
    elif os.path.exists(
            resource_filename(__name__, os.path.join("scoring_files",
                                                     json_conf["pick"]["scoring_file"]))):
        json_conf["pick"]["scoring_file"] = resource_filename(
            __name__,
            os.path.join("scoring_files",
                         json_conf["pick"]["scoring_file"]))
    else:
        raise InvalidJson(
            "Scoring file not found: {0}".format(
                json_conf["pick"]["scoring_file"]))
    return json_conf, overwritten


def check_json(json_conf, simple=False, external_dict=None, logger=None):

    """
    Wrapper for the various checks performed on the configuration file.

    :param json_conf: The dicitonary loaded from the configuration file.
    :type json_conf: dict

    :param simple: boolean flag indicating whether we desire
                   the simplified version of the configuration, or not.
    :type simple: bool

    :param external_dict: optional external dictionary with values to pass to the configuration.
    :type external_dict: (dict|None)

    :param logger: external logger instance
    :type logger: Logger

    :return json_conf
    :rtype: dict
    """

    # validator = extend_with_default(jsonschema.Draft4Validator)
    #
    # blue_print = os.path.join(
    #     os.path.dirname(os.path.abspath(__file__)),
    #     "configuration_blueprint.json"
    # )
    #
    # with open(blue_print) as blue:
    #     blue_print = json.load(blue)

    if not isinstance(logger, Logger):
        logger = create_default_logger("check_json")

    try:
        validator = create_validator(simple=simple)

        # config_folder = os.path.dirname(os.path.abspath(__file__))

        # This will check for consistency and add the default
        # values if they are missing
        validator.validate(json_conf)
        assert "files" in json_conf["pick"]

        overwritten = False

        json_conf, overwritten = _check_scoring_file(json_conf, logger)

        if json_conf["pick"]["scoring_file"].endswith(("yaml", "json")):
            with open(json_conf["pick"]["scoring_file"]) as scoring_file:
                if json_conf["pick"]["scoring_file"].endswith("yaml"):
                    scoring = yaml.load(scoring_file)
                else:
                    scoring = json.load(scoring_file)
            assert isinstance(json_conf, dict) and isinstance(scoring, dict),\
                (type(json_conf), type(scoring))
            json_conf = merge_dictionaries(json_conf, scoring)
            json_conf = check_all_requirements(json_conf)
            json_conf = check_scoring(json_conf)

        elif json_conf["pick"]["scoring_file"].endswith(("model", "pickle")):
            with open(json_conf["pick"]["scoring_file"], "rb") as forest:
                scoring = pickle.load(forest)
                assert isinstance(scoring, dict)
                assert "scoring" in scoring and isinstance(scoring["scoring"], (RandomForestRegressor, RandomForestClassifier))
                del scoring["scoring"]
                json_conf = merge_dictionaries(json_conf, scoring)
                json_conf = check_all_requirements(json_conf)
        else:
            raise InvalidJson(
                "Invalid scoring file: {0}".format(
                    json_conf["pick"]["scoring_file"]))

        json_conf["__loaded_scoring"] = json_conf["pick"]["scoring_file"]

        if external_dict is not None:
            if not isinstance(external_dict, dict):
                raise TypeError("Passed an invalid external dictionary, type {}".format(
                    type(external_dict)))
            json_conf = merge_dictionaries(json_conf, external_dict)

        json_conf = check_db(json_conf)
        # json_conf = check_blast(json_conf, json_file)
        validator.validate(json_conf)
        # json_conf["prepare"]["canonical"] = tuple([tuple(_) for _
        #                                            in json_conf["prepare"]["canonical"]])
        if not json_conf["multiprocessing_method"]:
            json_conf["multiprocessing_method"] = get_start_method()
    except Exception as exc:
        logger.exception(exc)
        raise

    if overwritten is True:
        logger.debug("Scoring parameters: {}".format("\n".join(["\n"] + [
            "{}: {}".format(_, json_conf["scoring"][_]) for _ in json_conf["scoring"].keys()])))

    return json_conf


def to_json(string, simple=False, logger=None):
    """
    Function to serialise the JSON for configuration and check its consistency.

    :param string: the configuration file name.
    :type string: (str | None | dict)

    :param simple: boolean flag indicating whether we desire
                   the simplified version of the configuration, or not.
    :type simple: bool

    :param logger: optional logger to be used.
    :type logger: Logger

    :rtype: dict
    """

    if not isinstance(logger, Logger):
        logger = create_default_logger("to_json")

    try:
        if string is None or string == '' or string == dict():
            json_dict = dict()
            string = os.path.join(os.path.abspath(os.getcwd()), "mikado.json")
        else:
            string = os.path.abspath(string)
            if not os.path.exists(string) or os.stat(string).st_size == 0:
                raise InvalidJson("JSON file {} not found!".format(string))
            with open(string) as json_file:
                if string.endswith(".yaml"):
                    json_dict = yaml.load(json_file)
                else:
                    json_dict = json.load(json_file)
        json_dict["filename"] = string
        # json_dict = frozendict(check_json(json_dict, simple=simple))
        json_dict = check_json(json_dict, simple=simple, logger=logger)
    except Exception as exc:
        logger.exception(exc)
        raise

    return json_dict
