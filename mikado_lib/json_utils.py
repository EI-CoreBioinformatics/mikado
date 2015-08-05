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
import subprocess


def check_log(json_conf):
    """
    :param json_conf: configuration dictionary to check.
    :type json_conf: dict

    This function checks that the logging configuration section of the JSON/YAML is properly formatted
    and adds missing values if necessary.

    :return: json_conf
    :rtype dict
    """

    if "log_settings" not in json_conf:
        json_conf["log_settings"] = dict()
    if "log" not in json_conf["log_settings"]:
        json_conf["log_settings"]["log"] = None
    if "log_level" not in json_conf["log_settings"]:
        json_conf["log_settings"]["log_level"] = "WARN"
    else:
        valid_levels = ["INFO", "WARN", "ERROR", "CRITICAL", "DEBUG"]
        if json_conf["log_settings"]["log_level"] not in valid_levels:
            raise mikado_lib.exceptions.InvalidJson("Invalid log level: {0}\mValid levels: {1}".format(
                json_conf["log_level"],
                "\n\t".join(valid_levels)
            ))
    return json_conf


def merge_dictionaries(a, b, path=None):
    """Recursive function to merge two dictionaries.

    :param a: first dictionary
    :type a: dict

    :param b: second dictionary
    :type b: dict

    :param path: list to be updated during recursion to indicate that we are in a secondary node
    :type path: list(str)

    Source: http://stackoverflow.com/questions/7204805/dictionaries-of-dictionaries-merge
    """

    if path is None:
        path = []
    for key in b:
        if key in a:
            if isinstance(a[key], dict) and isinstance(b[key], dict):
                merge_dictionaries(a[key], b[key], path + [str(key)])
            else:
                pass  # same leaf value
        else:
            a[key] = b[key]
    return a


def check_alternative_splicing(json_conf):
    """
    :param json_conf: configuration dictionary to check.
    :type json_conf: dict

    This function checks the options for the alternative splicing reporting.
    If no option is present in the configuration file, the basic "alternative_splicing": {"report" = False }
    will be added to the json_conf.

    Possible options:

    - min_cds_overlap:    Amount of CDS recall on the primary transcript. It can be expressed as values b/w 01 and 1,
                          or 1 to 100. Default: 0.
    - max_utr_length      Maximum UTR length of the AS isoforms.
    - max_isoforms        Maximum number of isoforms per Locus.

    :return: json_conf
    :rtype dict
    """

    if "alternative_splicing" not in json_conf:
        json_conf["alternative_splicing"] = dict()

    if "report" not in json_conf["alternative_splicing"]:
        json_conf["alternative_splicing"]["report"] = False
    else:
        assert type(json_conf["alternative_splicing"]["report"]) is bool
        if "min_cds_overlap" not in json_conf["alternative_splicing"]:
            json_conf["alternative_splicing"]["min_cds_overlap"] = 0
        else:
            assert type(json_conf["alternative_splicing"]["min_cds_overlap"]) in (float, int)
            if json_conf["alternative_splicing"]["min_cds_overlap"] < 0 or json_conf["alternative_splicing"][
                    "min_cds_overlap"] > 100:
                raise InvalidJson("Invalid percentage value for min_cds_overlap: {0}".format(
                    json_conf["alternative_splicing"]["min_cds_overlap"]))
            if 1 < json_conf["alternative_splicing"]["min_cds_overlap"] <= 100:
                json_conf["alternative_splicing"]["min_cds_overlap"] /= 100
        if "max_isoforms" not in json_conf["alternative_splicing"]:
            json_conf["alternative_splicing"]["max_isoforms"] = 10000
        else:
            assert type(json_conf["alternative_splicing"]["max_isoforms"]) is int and \
                json_conf["alternative_splicing"]["max_isoforms"] >= 1
        if "max_utr_length" not in json_conf["alternative_splicing"]:
            json_conf["alternative_splicing"]["max_utr_length"] = float("Inf")
        else:
            assert type(json_conf["alternative_splicing"]["max_utr_length"]) is int and \
                json_conf["alternative_splicing"]["max_utr_length"] >= 0
        if "max_fiveutr_length" not in json_conf["alternative_splicing"]:
            json_conf["alternative_splicing"]["max_fiveutr_length"] = float("Inf")
        else:
            assert type(json_conf["alternative_splicing"]["max_fiveutr_length"]) is int and \
                json_conf["alternative_splicing"]["max_fiveutr_length"] >= 0
        if "max_threeutr_length" not in json_conf["alternative_splicing"]:
            json_conf["alternative_splicing"]["max_threeutr_length"] = float("Inf")
        else:
            assert type(json_conf["alternative_splicing"]["max_threeutr_length"]) is int and \
                json_conf["alternative_splicing"]["max_threeutr_length"] >= 0
        if "keep_retained_introns" not in json_conf["alternative_splicing"]:
            json_conf["alternative_splicing"]["keep_retained_introns"] = True
        else:
            assert type(json_conf["alternative_splicing"]["keep_retained_introns"]) is bool

    return json_conf


def check_chimera_split(json_conf):
    """
    :param json_conf: configuration dictionary to check.
    :type json_conf: dict

    Function to check the "chimera_split" section of the configuration.

    :return: json_conf
    :rtype dict
    """

    assert ("execute" in json_conf["chimera_split"] and type(json_conf["chimera_split"]["execute"]) is bool)
    if json_conf["chimera_split"]["execute"] is True:
        assert ("blast_check" in json_conf["chimera_split"] and type(json_conf["chimera_split"]["blast_check"]) is bool)
        if json_conf["chimera_split"]["blast_check"] is True:
            assert "blast_params" in json_conf["chimera_split"]
            if "evalue" in json_conf["chimera_split"]["blast_params"]:
                assert type(json_conf["chimera_split"]["blast_params"]["evalue"]) in (float, int) or \
                    json_conf["chimera_split"]["blast_params"]["evalue"] is None
            else:
                json_conf["chimera_split"]["blast_params"]["evalue"] = None
            if "hsp_evalue" in json_conf["chimera_split"]["blast_params"]:
                assert type(json_conf["chimera_split"]["blast_params"]["hsp_evalue"]) in (float, int) or \
                    json_conf["chimera_split"]["blast_params"]["hsp_evalue"] is None
                if json_conf["chimera_split"]["blast_params"]["evalue"] is None:
                    json_conf["chimera_split"]["blast_params"]["evalue"] = json_conf["chimera_split"][
                        "blast_params"]["hsp_evalue"]
                elif json_conf["chimera_split"]["blast_params"]["evalue"] > \
                        json_conf["chimera_split"]["blast_params"]["hsp_evalue"]:
                    raise mikado_lib.exceptions.InvalidJson(
                        "Maximum HSP evalues cannot be higher than global e-values.")
            else:
                if json_conf["chimera_split"]["blast_params"]["evalue"] is None:
                    json_conf["chimera_split"]["blast_params"]["hsp_evalue"] = None
                else:
                    json_conf["chimera_split"]["blast_params"]["hsp_evalue"] = \
                        json_conf["chimera_split"]["blast_params"]["evalue"]
            if "max_target_seqs" in json_conf["chimera_split"]["blast_params"]:
                assert type(json_conf["chimera_split"]["blast_params"]["max_target_seqs"]) is int or \
                    json_conf["chimera_split"]["blast_params"]["max_target_seqs"] is None
            else:
                json_conf["chimera_split"]["blast_params"]["max_target_seqs"] = None
            if "minimal_hsp_overlap" in json_conf["chimera_split"]["blast_params"]:
                assert type(json_conf["chimera_split"]["blast_params"]["minimal_hsp_overlap"]) is float and \
                    0 <= json_conf["chimera_split"]["blast_params"]["minimal_hsp_overlap"] <= 1
            else:
                json_conf["chimera_split"]["blast_params"]["minimal_hsp_overlap"] = 0
            if "leniency" not in json_conf["chimera_split"]["blast_params"]:
                json_conf["chimera_split"]["blast_params"]["leniency"] = "STRINGENT"
            else:
                assert json_conf["chimera_split"]["blast_params"]["leniency"] in ("STRINGENT", "PERMISSIVE", "LENIENT")

    return json_conf


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

    parameters_found = set()
    parameters_not_found = []
    double_parameters = []
    invalid_filter = set()
    available_metrics = Transcript.get_available_metrics()
    if "scoring" not in json_conf or len(json_conf["scoring"].keys()) == 0:
        raise mikado_lib.exceptions.InvalidJson("No parameters specified for scoring!")

    for parameter in json_conf["scoring"]:
        if parameter not in available_metrics:
            parameters_not_found.append(parameter)
        if parameter in parameters_found:
            double_parameters.add(parameter)
        if "filter" in json_conf["scoring"][parameter]:
            conf = json_conf["scoring"][parameter]["filter"]
            if "operator" not in conf or "value" not in conf:
                invalid_filter.add(parameter)
            elif conf["operator"] not in ("gt", "ge", "eq", "lt", "le", "ne", "in", "not in"):
                invalid_filter.add(parameter)
        if "rescaling" not in json_conf["scoring"][parameter]:
            raise mikado_lib.exceptions.UnrecognizedRescaler(
                "No rescaling specified for {0}. Must be one among \"max\",\"min\", and \"target\".".format(parameter))
        elif json_conf["scoring"][parameter]["rescaling"] not in ("max", "min", "target"):
            raise mikado_lib.exceptions.UnrecognizedRescaler(
                "Invalid rescaling specified for {0}. Must be one among \"max\",\"min\", and \"target\".".format(
                    parameter))
        elif json_conf["scoring"][parameter]["rescaling"] == "target":
            if "value" not in json_conf["scoring"][parameter]:
                message = "Target rescaling requested for {0}".format(parameter)
                message += " but no target value specified. Please specify it with the \"value\" keyword."
                raise mikado_lib.exceptions.UnrecognizedRescaler(message)

            json_conf["scoring"][parameter]["value"] = float(json_conf["scoring"][parameter]["value"])

        if "multiplier" not in json_conf["scoring"][parameter]:
            json_conf["scoring"][parameter]["multiplier"] = 1
        else:
            if type(json_conf["scoring"][parameter]["multiplier"]) not in (float, int) or \
                    json_conf["scoring"][parameter]["multiplier"] == 0:
                raise mikado_lib.exceptions.InvalidJson(
                    "Invalid multiplier: {0}".format(json_conf["scoring"][parameter]["multiplier"]))
            json_conf["scoring"][parameter]["multiplier"] = float(json_conf["scoring"][parameter]["multiplier"])

    if len(parameters_not_found) > 0 or len(double_parameters) > 0 or len(invalid_filter) > 0:
        err_message = ''
        if len(parameters_not_found) > 0:
            err_message = "The following parameters, present in the JSON file, are not available!\n\t{0}\n".format(
                "\n\t".join(parameters_not_found))
        if len(double_parameters) > 0:
            err_message += "The following parameters have been specified more than once, please correct:"
            err_message += "\n\t{0}".format("\n\t".join(list(double_parameters)))
        if len(invalid_filter) > 0:
            err_message += "The following parameters have an invalid filter, please correct:\n\t{0}".format(
                "\n\t".join(list(invalid_filter)))
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

    if "requirements" in json_conf:
        # Check that the parameters are valid
        if "parameters" not in json_conf["requirements"]:
            raise mikado_lib.exceptions.InvalidJson("The requirements field must have a \"parameters\" subfield!")
        for key in json_conf["requirements"]["parameters"]:
            key_name = key.split(".")[0]
            if key_name not in available_metrics:
                parameters_not_found.append(key_name)
            if "operator" not in json_conf["requirements"]["parameters"][key]:
                raise mikado_lib.exceptions.InvalidJson("No operator provided for requirement {0}".format(key))
            elif "value" not in json_conf["requirements"]["parameters"][key]:
                raise mikado_lib.exceptions.InvalidJson("No value provided for requirement {0}".format(key))
            elif json_conf["requirements"]["parameters"][key]["operator"] not in (
                    "gt", "ge", "eq", "lt", "le", "ne", "in", "not in"):
                raise mikado_lib.exceptions.UnrecognizedOperator(
                    "Unrecognized operator: {0}".format(json_conf["requirements"]["parameters"][key]["operator"]))
            json_conf["requirements"]["parameters"][key]["name"] = key_name
        if len(parameters_not_found) > 0:
            raise mikado_lib.exceptions.InvalidJson(
                "The following parameters, selected for filtering, are invalid:\n\t{0}".format(
                    "\n\t".join(parameters_not_found)
                ))
        if "expression" not in json_conf["requirements"]:  # Create automatically a filtering expression
            json_conf["requirements"]["expression"] = " and ".join(list(json_conf["requirements"]["parameters"].keys()))
            keys = json_conf["requirements"]["parameters"].keys()
            newexpr = json_conf["requirements"]["expression"][:]
        else:
            # Parse the filtering expression, verify that it is syntactically correct
            if type(json_conf["requirements"]["expression"]) is list:
                json_conf["requirements"]["expression"] = " ".join(json_conf["requirements"]["expression"])
            newexpr = json_conf["requirements"]["expression"][:]
            keys = list(filter(lambda x: x not in ("and", "or", "not", "xor"),
                               re.findall("([^ ()]+)", json_conf["requirements"]["expression"])))
            diff_params = set.difference(set(keys), set(json_conf["requirements"]["parameters"].keys()))
            if len(diff_params) > 0:
                raise mikado_lib.exceptions.InvalidJson(
                    "Expression and required parameters mismatch:\n\t{0}".format("\n\t".join(list(diff_params))))
        for key in keys:  # Create the final expression
            newexpr = re.sub(key, "evaluated[\"{0}\"]".format(key), newexpr)
        json_conf["requirements"]["expression"] = newexpr

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

    '''Function to check the BLAST section of the JSON and eventually perform the database indexing.'''

    if "blast" not in json_conf:
        json_conf["blast"] = dict()
        json_conf["blast"]["execute"] = False
        return json_conf

    if json_conf["blast"]["execute"] is False:
        return json_conf

    if "program" not in json_conf["blast"]:
        raise mikado_lib.exceptions.InvalidJson("No BLAST program specified.")
    elif os.path.basename(json_conf["blast"]["program"]) not in ("blastn", "blastx", "tblastx"):
        raise mikado_lib.exceptions.InvalidJson("""Invalid BLAST program specified: {0}.
        Supported options: blastn, blastx, tblastx.""")
    if os.path.dirname(json_conf["blast"]["program"]) == "":
        program = spawn.find_executable(json_conf["blast"]["program"])  # @UndefinedVariable
    else:
        try:
            program = os.path.abspath(json_conf["blast"]["program"])
        except OSError:
            program = None
    if program is None:
        raise mikado_lib.exceptions.InvalidJson(
            "The selected BLAST program {0} has not been found on this system!".format(json_conf["blast"]["program"]))
    json_conf["blast"]["program"] = program

    if "evalue" not in json_conf["blast"]:
        json_conf["blast"]["evalue"] = 10
    else:
        if type(json_conf["blast"]["evalue"]) not in (float, int) or 0 > json_conf["blast"]["evalue"]:
            raise mikado_lib.exceptions.InvalidJson("Invalid evalue: {0}".format(json_conf["blast"]["evalue"]))
    if "max_target_seqs" in json_conf["blast"]:
        assert type(json_conf["blast"]["max_target_seqs"]) is int
    if "database" not in json_conf["blast"]:
        raise mikado_lib.exceptions.InvalidJson("No BLAST database provided!")
    json_conf["blast"]["database"] = os.path.abspath(json_conf["blast"]["database"])
    if not os.path.exists(json_conf["blast"]["database"]):
        db = os.path.join(
            os.path.dirname(json_file),
            os.path.basename(json_conf["blast"]["database"])
        )
        if not os.path.exists(db):
            raise mikado_lib.exceptions.InvalidJson(
                "I need a valid BLAST database! This file does not exist:\n{0}".format(json_conf["blast"]["database"]))
        else:
            json_conf["blast"]["database"] = os.path.abspath(db)
    else:
        json_conf["blast"]["database"] = os.path.abspath(json_conf["blast"]["database"])
    makeblastdb_cmd = os.path.join(os.path.dirname(json_conf["blast"]["program"]), "makeblastdb")
    assert os.path.exists(makeblastdb_cmd)
    retcode = 0
    if os.path.basename(json_conf["blast"]["program"]) == "blastx" and not os.path.exists(
            "{0}.pog".format(json_conf["blast"]["database"])):
        retcode = subprocess.call(
            "{0} -in {1} -dbtype prot -parse_seqids".format(makeblastdb_cmd, json_conf["blast"]["database"]),
            shell=True)
    elif os.path.basename(json_conf["blast"]["program"]) in ("blastn", "tblastx") and not os.path.exists(
            "{0}.nog".format(json_conf["blast"]["database"])):
        retcode = subprocess.call(
            "{0} -in {1} -dbtype nucl -parse_seqids".format(makeblastdb_cmd, json_conf["blast"]["database"]),
            shell=True)
    if retcode != 0:
        raise OSError("BLAST indexing failed.")

    return json_conf


def check_orf_loading(json_conf):
    """
    :param json_conf: configuration dictionary to check.
    :type json_conf: dict

    Function to check the options related to ORF loading.

    :return: json_conf
    :rtype dict
    """

    if "orf_loading" not in json_conf:
        json_conf["orf_loading"] = dict()
        json_conf["orf_loading"]["strand_specific"] = False
        json_conf["orf_loading"]["minimal_secondary_orf_length"] = 0
    else:
        if "strand_specific" not in json_conf:
            json_conf["orf_loading"]["strand_specific"] = False
        else:
            if not type(json_conf["orf_loading"]["strand_specific"]) is bool:
                raise mikado_lib.exceptions.InvalidJson(
                    "Invalid strand_specific value: {0}".format(json_conf["orf_loading"]["strand_specific"]))
        if "minimal_secondary_orf_length" not in json_conf["orf_loading"]:
            json_conf["orf_loading"]["minimal_secondary_orf_length"] = 0
        else:
            if not type(json_conf["orf_loading"]["minimal_secondary_orf_length"]) is int:
                raise mikado_lib.exceptions.InvalidJson("Invalid minimal_secondary_orf_length value: {0}".format(
                    json_conf["orf_loading"]["minimal_secondary_orf_length"]))

    return json_conf


def check_run_options(json_conf):
    """
    :param json_conf: configuration dictionary to check.
    :type json_conf: dict

    Function to check general run options such as number of threads to use.

    :return: json_conf
    :rtype dict
    """

    if "run_options" not in json_conf:
        json_conf["run_options"] = dict()
    if "purge" not in json_conf["run_options"]:
        json_conf["run_options"]["purge"] = False
    else:
        assert type(json_conf["run_options"]["purge"]) is bool
    if "preload" not in json_conf["run_options"]:
        json_conf["run_options"]["preload"] = False
    else:
        assert type(json_conf["run_options"]["preload"]) is bool
    if "exclude_cds" not in json_conf["run_options"]:
        json_conf["run_options"]["exclude_cds"] = False
    else:
        assert type(json_conf["run_options"]["exclude_cds"]) is bool
    if "remove_overlapping_fragments" not in json_conf["run_options"]:
        json_conf["run_options"]["remove_overlapping_fragments"] = False
    else:
        assert type(json_conf["run_options"]["remove_overlapping_fragments"]) is bool
    if "fragments_maximal_cds" not in json_conf["run_options"]:
        json_conf["run_options"]["fragments_maximal_cds"] = 100
    else:
        assert type(json_conf["run_options"]["fragments_maximal_cds"]) is int and \
            json_conf["run_options"]["fragments_maximal_cds"] >= 0

    if "threads" not in json_conf["run_options"]:
        json_conf["run_options"]["threads"] = 1
    else:
        assert type(json_conf["run_options"]["threads"]) is int

    if "subloci_from_cds_only" not in json_conf["run_options"]:
        json_conf["run_options"]["subloci_from_cds_only"] = False
    else:
        assert type(json_conf["run_options"]["subloci_from_cds_only"]) is bool

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

    if "scoring_file" in json_conf:
        if os.path.exists(os.path.abspath(json_conf["scoring_file"])):
            json_conf["scoring_file"] = os.path.abspath(json_conf["scoring_file"])
        else:
            if os.path.exists(os.path.join(os.path.dirname(json_conf["filename"]), json_conf["scoring_file"])):
                json_conf["scoring_file"] = os.path.join(os.path.dirname(json_conf["filename"]),
                                                         json_conf["scoring_file"])
            else:
                raise mikado_lib.exceptions.InvalidJson("Scoring file not found: {0}".format(json_conf["scoring_file"]))

        with open(json_conf["scoring_file"]) as scoring_file:
            if json_conf["scoring_file"].endswith("yaml"):
                scoring = yaml.load(scoring_file)
            else:
                scoring = json.load(scoring_file)
        assert type(json_conf) is dict and type(scoring) is dict, (type(json_conf), type(scoring))
        json_conf = merge_dictionaries(json_conf, scoring)

    if "db" not in json_conf:
        raise mikado_lib.exceptions.InvalidJson("No database specified.")
    if "dbtype" not in json_conf:
        raise mikado_lib.exceptions.InvalidJson("DB type not specified.")
    if json_conf["dbtype"] not in ("sqlite", "mysql", "postgresql"):
        raise mikado_lib.exceptions.InvalidJson(
            "Invalid DB type: {0}. At the moment we support sqlite, mysql, postgresql".format(json_conf["dbtype"]))

    if json_conf["dbtype"] in ("mysql", "postgresql"):
        if "dbhost" not in json_conf:
            raise mikado_lib.exceptions.InvalidJson(
                "No host specified for the {0} database!".format(json_conf["dbtype"]))
        if "dbuser" not in json_conf:
            raise mikado_lib.exceptions.InvalidJson(
                "No user specified for the {0} database!".format(json_conf["dbtype"]))
        if "dbpasswd" not in json_conf or json_conf['dbpasswd'] is None:
            json_conf["dbpasswd"] = ''
        if "dbport" in json_conf and json_conf["dbport"] is not None:
            if not type(json_conf["dbport"]) is int:
                raise mikado_lib.exceptions.InvalidJson(
                    "Invalid type for dbport: {0}".format(type(json_conf["dbport"])))

        else:
            # Default ports
            if json_conf["dbtype"] == "mysql":
                json_conf["dbport"] = 3306
            else:
                json_conf["dbport"] = 5432

    if "input" not in json_conf:
        json_conf["input"] = None
    else:
        assert os.path.exists(json_conf["input"]) and os.path.isfile(json_conf["input"])

    if "source" not in json_conf:
        json_conf["source"] = "Mikado"

    for prefix in ["", "mono", "sub"]:
        key = "{0}loci_out".format(prefix)
        if key not in json_conf:
            json_conf[key] = None

    json_conf = check_blast(json_conf, json_file)
    json_conf = check_requirements(json_conf)
    json_conf = check_scoring(json_conf)
    json_conf = check_orf_loading(json_conf)
    json_conf = check_chimera_split(json_conf)
    json_conf = check_run_options(json_conf)
    json_conf = check_log(json_conf)
    json_conf = check_alternative_splicing(json_conf)
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
