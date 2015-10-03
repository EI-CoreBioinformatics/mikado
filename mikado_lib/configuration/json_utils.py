#!/usr/bin/env python3
# coding: utf-8

"""
This module defines the functionalities needed to verify the integrity and completeness
of Mikado configuration files. Missing values are replaced with default ones,
while existing values are checked for type and consistency.
"""

import os.path
import re
from distutils import spawn
import yaml
from mikado_lib.exceptions import InvalidJson
from mikado_lib.loci_objects.transcript import Transcript
import mikado_lib.exceptions
import json
import sys
import subprocess
import jsonschema


def extend_with_default(validator_class):
    """
    Function to extend the normal validation classes for jsonschema
    so that they also set the default values provided inside the schema
    itself. Source:
    https://python-jsonschema.readthedocs.org/en/latest/faq/?highlight=default
    :param validator_class: the validator class to extend (e.g. Draft4Validator)
    :return:
    """
    validate_properties = validator_class.VALIDATORS["properties"]

    def set_default(instance, properties):
        """
        Recursive function that sets the default parameters inside "object"
        types for the dictionary instance. It also loads comments, if available.
        :param instance: a dictionary instance
        :param properties: the "properties" dictionary inside the schema
        :return:
        """
        for prop, subschema in properties.items():
            if "default" in subschema:
                instance.setdefault(prop, subschema["default"])
            elif prop not in instance:
                if "type" not in subschema: continue
                elif subschema["type"] == "object":
                    instance[prop] = dict()
                    if "comment" in subschema:
                        instance[prop].setdefault("comment", subschema["comment"])
                    instance[prop] = set_default(instance[prop], subschema["properties"])
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
        for error in validate_properties(
            validator, properties, instance, schema,
        ):
            yield error
        instance = set_default(instance, properties)

    return jsonschema.validators.extend(
        validator_class, {"properties": set_defaults},
    )


def merge_dictionaries(dict_a, dict_b, path=None):
    """Recursive function to merge two dictionaries.

    :param dict_a: first dictionary
    :type dict_a: dict

    :param dict_b: second dictionary
    :type dict_b: dict

    :param path: list to be updated during recursion to indicate that we are in a secondary node
    :type path: list(str)

    Source: http://stackoverflow.com/questions/7204805/dictionaries-of-dictionaries-merge
    """

    if path is None:
        path = []
    for key in dict_b:
        if key in dict_a:
            if isinstance(dict_a[key], dict) and isinstance(dict_b[key], dict):
                merge_dictionaries(
                    dict_a[key],
                    dict_b[key], path + [str(key)])
            else:
                pass  # same leaf value
        else:
            dict_a[key] = dict_b[key]
    return dict_a


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

    scoring_schema = json.load(open(os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "scoring_blueprint.json"
    )))

    parameters_found = set()
    parameters_not_found = []
    double_parameters = []
    invalid_filter = set()
    available_metrics = Transcript.get_available_metrics()
    if "scoring" not in json_conf or len(json_conf["scoring"].keys()) == 0:
        raise mikado_lib.exceptions.InvalidJson("No parameters specified for scoring!")

    validator = extend_with_default(jsonschema.Draft4Validator)

    for parameter in json_conf["scoring"]:
        if parameter not in available_metrics:
            parameters_not_found.append(parameter)
        if parameter in parameters_found:
            double_parameters.append(parameter)
        if not jsonschema.Draft4Validator(scoring_schema).is_valid(json_conf["scoring"][parameter]):
            errors = list(jsonschema.Draft4Validator(scoring_schema).iter_errors(
                json_conf["scoring"][parameter]))
            raise mikado_lib.exceptions.InvalidJson("Invalid scoring for {}:\n{}".format(
                parameter,"\n".join(errors)))
        try:
            validator(scoring_schema).validate(json_conf["scoring"][parameter])
        except Exception as err:
            raise ValueError(parameter, err)
        if json_conf["scoring"][parameter]["rescaling"] == "target":
            if "value" not in json_conf["scoring"][parameter]:
                message = """Target rescaling requested for {0} but no target value specified.
                    Please specify it with the \"value\" keyword.\n{1}""".format(parameter,
                                                                                 json_conf["scoring"][parameter])
                raise mikado_lib.exceptions.UnrecognizedRescaler(message)

    if len(parameters_not_found) > 0 or len(double_parameters) > 0 or len(invalid_filter) > 0:
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
        raise mikado_lib.exceptions.InvalidJson(err_message)

    return json_conf


def check_requirements(json_conf):
    """
    :param json_conf: configuration dictionary to check.
    :type json_conf: dict

    Function to check the "requirements" section of the configuration.
    Each filtering function will be checked for:
    - validity of the expression (it can be interpreted by Mikado)
    - validity of the parameter (it is a valid Metric)

    :return: json_conf
    :rtype dict
    """

    available_metrics = Transcript.get_available_metrics()
    parameters_not_found = []
    require_schema = json.load(open(os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "requirements_blueprint.json"
    )))

    if "requirements" in json_conf:
        # Check that the parameters are valid
        if "parameters" not in json_conf["requirements"]:
            raise mikado_lib.exceptions.InvalidJson(
                "The requirements field must have a \"parameters\" subfield!")
        for key in json_conf["requirements"]["parameters"]:
            key_name = key.split(".")[0]
            if key_name not in available_metrics:
                parameters_not_found.append(key_name)
                continue
            if not jsonschema.Draft4Validator(require_schema["definitions"]["parameter"]).is_valid(
                    json_conf["requirements"]["parameters"][key]):
                errors = list(jsonschema.Draft4Validator(require_schema).iter_errors(
                    json_conf["requirements"]["parameters"][key]
                ))
                raise mikado_lib.exceptions.InvalidJson("Invalid parameter for {0}: \n{1}".format(
                    key, errors))

            json_conf["requirements"]["parameters"][key]["name"] = key_name

        if len(parameters_not_found) > 0:
            raise mikado_lib.exceptions.InvalidJson(
                "The following parameters, selected for filtering, are invalid:\n\t{0}".format(
                    "\n\t".join(parameters_not_found)
                ))

        # Create automatically a filtering expression
        if "expression" not in json_conf["requirements"]:
            json_conf["requirements"]["expression"] = " and ".join(
                list(json_conf["requirements"]["parameters"].keys()))
            keys = json_conf["requirements"]["parameters"].keys()
            newexpr = json_conf["requirements"]["expression"][:]
        else:

            if not jsonschema.Draft4Validator(require_schema["definitions"]["expression"]).is_valid(
                json_conf["requirements"]["expression"]
            ):
                raise mikado_lib.exceptions.InvalidJson("Invalid expression field")
            expr = " ".join(json_conf["requirements"]["expression"])
            newexpr = expr[:]

            keys = list(key for key in re.findall(
                "([^ ()]+)", expr) if key not in ("and", "or", "not", "xor"))

            diff_params = set.difference(
                set(keys), set(json_conf["requirements"]["parameters"].keys()))

            if len(diff_params) > 0:
                raise mikado_lib.exceptions.InvalidJson(
                    "Expression and required parameters mismatch:\n\t{0}".format(
                        "\n\t".join(list(diff_params))))

        for key in keys:  # Create the final expression
            newexpr = re.sub(key, "evaluated[\"{0}\"]".format(key), newexpr)
        json_conf["requirements"]["expression"] = newexpr

    if "soft_requirements" not in json_conf:
        json_conf["soft_requirements"] = dict()
        json_conf["soft_requirements"]["intron_range"] = (0, sys.maxsize)
    if "intron_range" not in json_conf["soft_requirements"]:
        raise InvalidJson("No intron range found!")
        # json_conf["soft_requirements"]["intron_range"] = (0, sys.maxsize)
    else:
        try:
            minimal, maximal = (
                int(_) for _ in json_conf["soft_requirements"]["intron_range"])
            json_conf["soft_requirements"]["intron_range"] = (minimal, maximal)
        except Exception:
            raise InvalidJson("Invalid intron range: {0}".format(
                json_conf["soft_requirements"]["intron_range"]
            ))

    return json_conf


def check_blast(json_conf, json_file):
    """
    :param json_conf: configuration dictionary to check.
    :type json_conf: dict

    :param json_file: the original JSON file name (necessary to derive its parent folder)
    :type json_file: str

    Function to check the optional "blast" section of the configuration.
    It will be created (with "execute" set to False) if it is missing.

    :return: json_conf
    :rtype dict
    """

    if "blast" not in json_conf:
        json_conf["blast"] = dict()
        json_conf["blast"]["execute"] = False
        json_conf["blast"]["max_target_seqs"] = sys.maxsize
        return json_conf

    if json_conf["blast"]["execute"] is False:
        return json_conf

    if "program" not in json_conf["blast"]:
        raise mikado_lib.exceptions.InvalidJson("No BLAST program specified.")
    elif os.path.basename(json_conf["blast"]["program"]) not in ("blastn", "blastx", "tblastx"):
        raise mikado_lib.exceptions.InvalidJson("""Invalid BLAST program specified: {0}.
        Supported options: blastn, blastx, tblastx.""")
    if os.path.dirname(json_conf["blast"]["program"]) == "":
        program = spawn.find_executable(json_conf["blast"]["program"])
    else:
        try:
            program = os.path.abspath(json_conf["blast"]["program"])
        except OSError:
            program = None
    if program is None:
        raise mikado_lib.exceptions.InvalidJson(
            "The selected BLAST program {0} has not been found on this system!".format(
                json_conf["blast"]["program"]))
    json_conf["blast"]["program"] = program

    if "evalue" not in json_conf["blast"]:
        json_conf["blast"]["evalue"] = 10
    else:
        evalue = json_conf["blast"]["evalue"]
        if not isinstance(evalue, (float, int)) or evalue < 0:
            raise mikado_lib.exceptions.InvalidJson(
                "Invalid evalue: {0}".format(evalue))
    if "max_target_seqs" in json_conf["blast"]:
        assert isinstance(json_conf["blast"]["max_target_seqs"], int)
    else:
        json_conf["blast"]["max_target_seqs"] = sys.maxsize
    if "database" not in json_conf["blast"]:
        raise mikado_lib.exceptions.InvalidJson("No BLAST database provided!")
    json_conf["blast"]["database"] = os.path.abspath(json_conf["blast"]["database"])
    if not os.path.exists(json_conf["blast"]["database"]):
        database = os.path.join(
            os.path.dirname(json_file),
            os.path.basename(json_conf["blast"]["database"])
        )
        if not os.path.exists(database):
            if os.path.exists("{0}.gz".format(database)):
                retcode = subprocess.call("gzip -dc {0}.gz > {0}".format(database), shell=True)
                if retcode != 0:
                    raise mikado_lib.exceptions.InvalidJson("Failed to decompress the BLAST database!")
            else:
                raise mikado_lib.exceptions.InvalidJson(
                    """I need a valid BLAST database! This file does not exist:
                    {0}""".format(json_conf["blast"]["database"]))
        else:
            json_conf["blast"]["database"] = os.path.abspath(database)
    else:
        json_conf["blast"]["database"] = os.path.abspath(
            json_conf["blast"]["database"])

    makeblastdb_cmd = os.path.join(
        os.path.dirname(json_conf["blast"]["program"]), "makeblastdb")
    assert os.path.exists(makeblastdb_cmd)
    retcode = 0
    program = os.path.basename(json_conf["blast"]["program"])
    if program == "blastx" and not os.path.exists("{0}.pog".format(
            json_conf["blast"]["database"])):
        retcode = subprocess.call(
            "{0} -in {1} -dbtype prot -parse_seqids".format(
                makeblastdb_cmd, json_conf["blast"]["database"]), shell=True)
    elif program in ("blastn", "tblastx") and not os.path.exists("{0}.nog".format(
            json_conf["blast"]["database"])):
        retcode = subprocess.call(
            "{0} -in {1} -dbtype nucl -parse_seqids".format(
                makeblastdb_cmd, json_conf["blast"]["database"]), shell=True)
    if retcode != 0:
        raise OSError("BLAST indexing failed.")

    return json_conf


def check_db(json_conf):

    """
    Function to check the validity of the database options.
    :param json_conf:
    :return:
    """

    db_settings = json_conf["db_settings"]

    if json_conf["db_settings"]["dbtype"] in ("mysql", "postgresql"):
        if "dbhost" not in json_conf["db_settings"]:
            raise mikado_lib.exceptions.InvalidJson(
                "No host specified for the {0} database!".format(json_conf["db_settings"]["dbtype"]))
        if "dbuser" not in json_conf["db_settings"]:
            raise mikado_lib.exceptions.InvalidJson(
                "No user specified for the {0} database!".format(json_conf["db_settings"]["dbtype"]))
        if json_conf["db_settings"]["dbport"] == 0:
            if json_conf["db_settings"]["dbtype"] == "mysql":
                json_conf["db_settings"]["dbport"] = 3306
            else:
                json_conf["db_settings"]["dbport"] = 5432

    return json_conf


def check_json(json_conf, json_file):

    """
    :param json_conf: The dicitonary loaded from the configuration file.
    :type json_conf: dict

    :param json_file: The original configuration file name.
    :type json_file: str


    Wrapper for the various checks performed on the configuration file.

    :return json_conf
    :rtype dict
    """

    validator = extend_with_default(jsonschema.Draft4Validator)

    blue_print = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "configuration_blueprint.json"
    )
    blue_print = json.load(open(blue_print))

    config_folder = os.path.dirname(os.path.abspath(__file__))

    # This will check for consistency and add the default
    # values if they are missing
    validator(blue_print).validate(json_conf)

    if "scoring_file" in json_conf:
        if os.path.exists(os.path.abspath(json_conf["scoring_file"])):
            json_conf["scoring_file"] = os.path.abspath(json_conf["scoring_file"])
        elif os.path.exists(
                    os.path.join(os.path.dirname(
                        json_conf["filename"]), json_conf["scoring_file"])):
                json_conf["scoring_file"] = os.path.join(
                    os.path.dirname(json_conf["filename"]),
                    json_conf["scoring_file"])
        elif os.path.exists(os.path.join(config_folder, json_conf["scoring_file"])):
                json_conf["scoring_file"] = os.path.join(config_folder, json_conf["scoring_file"])
        else:
            raise mikado_lib.exceptions.InvalidJson(
                "Scoring file not found: {0}".format(json_conf["scoring_file"]))

        with open(json_conf["scoring_file"]) as scoring_file:
            if json_conf["scoring_file"].endswith("yaml"):
                scoring = yaml.load(scoring_file)
            else:
                scoring = json.load(scoring_file)
        assert isinstance(json_conf, dict) and isinstance(scoring, dict),\
            (type(json_conf), type(scoring))

        json_conf = merge_dictionaries(json_conf, scoring)

    if "input" not in json_conf:
        json_conf["input"] = None
    else:
        assert os.path.exists(json_conf["input"]) and os.path.isfile(json_conf["input"])

    for prefix in ["", "mono", "sub"]:
        key = "{0}loci_out".format(prefix)
        if key not in json_conf:
            json_conf[key] = None

    json_conf = check_db(json_conf)
    json_conf = check_blast(json_conf, json_file)
    json_conf = check_requirements(json_conf)
    json_conf = check_scoring(json_conf)
    validator(blue_print).validate(json_conf)
    return json_conf


def to_json(string):
    """
    :param string: the configuration file name.
    :type string: str

    Function to serialise the JSON for configuration and check its consistency."""

    string = os.path.abspath(string)
    with open(string) as json_file:
        if string.endswith(".yaml"):
            json_dict = yaml.load(json_file)
        else:
            json_dict = json.load(json_file)
    json_dict["filename"] = string
    json_dict = check_json(json_dict, json_file.name)
    return json_dict
