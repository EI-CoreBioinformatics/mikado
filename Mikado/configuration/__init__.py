"""
This module defines the functions needed to check the sanity of the configuration file,
plus the JSON schemas for the configuration and scoring files.
"""
from typing import Union

import os
from marshmallow import fields
from .configuration import MikadoConfiguration
from .daijin_configuration import DaijinConfiguration
from . import configurator
import itertools
import textwrap
import rapidjson as json
import toml
import dataclasses
import yaml


__author__ = 'Luca Venturini'


def print_config(config: Union[MikadoConfiguration, DaijinConfiguration], out,
                 output_format="yaml", no_files=False, full=True):
    """
    Function to print out the configuration in TOML format, adding the descriptions as comments preceded by #.
    :param config: configuration
    :type config: (MikadoConfiguration|DaijinConfiguration)

    :param out: output handle
    :type out: [io.TextIOWrapper|io.TextIO]

    :param output_format: one of yaml, json or toml (case-insensitive)
    :type output_format: str

    :param no_files: boolean. If on, sections pertaining to files will be deleted from the output.
    :type no_files: bool
    """

    if not isinstance(output_format, str) or output_format.lower() not in ("yaml", "json", "toml"):
        raise ValueError("Unknown format: {}. I can only accept yaml, json or toml as options.")

    output_format = output_format.lower()
    config_dict = dataclasses.asdict(config)
    for key in ["scoring", "cds_requirements", "requirements", "not_fragmentary", "as_requirements"]:
        # Necessary otherwise we will be deleting fields from the *original* object!
        config_dict.pop(key, None)
    if no_files is True:
        for stage in ["pick", "prepare", "serialise"]:
            if "files" in config_dict[stage]:
                del config_dict[stage]["files"]
        del config_dict["reference"]
        del config_dict["db_settings"]

    if not full:
        config_dict = filter_config(config_dict, config)

    if output_format == "toml":
        print_toml_config(config_dict, config, out)
    elif output_format == "yaml":
        print_yaml_config(config_dict, config, out)
    elif output_format == "json":
        # Necessary otherwise we will be deleting fields from the *original* object!
        print(json.dumps(config_dict, indent=4, sort_keys=True), file=out)


def select_attribute_for_output(final_config_level: dict, attr_parent, attr_name, attr_value):
    is_nested = isinstance(attr_parent.Schema._declared_fields[attr_name], fields.Nested)
    is_required = attr_parent.Schema._declared_fields[attr_name].required

    # This is not going to happen if we are looking at a nested property
    if attr_parent.Schema._declared_fields[attr_name].default == attr_value:
        if is_required:
            final_config_level[attr_name] = attr_value
        else:
            return
    elif not is_nested:
        callable_and_equal_to_default = (callable(attr_parent.Schema._declared_fields[attr_name].default) and
                                         attr_parent.Schema._declared_fields[attr_name].default() == attr_value)
        if callable_and_equal_to_default and not is_required:
            return
        else:
            final_config_level[attr_name] = attr_value
    else:
        final_config_level[attr_name] = dict()
        for key, value in dataclasses.asdict(attr_value).items():
            select_attribute_for_output(
                final_config_level[attr_name], attr_value, key, getattr(attr_value, key))


def clear_dict(data: dict) -> dict:
    keys = list(data.keys())
    for key in keys:
        if isinstance(data[key], dict) and len(data[key]) > 0:
            data[key] = clear_dict(data[key])
        if isinstance(data[key], dict) and len(data[key]) == 0:
            del data[key]
        else:
            pass
    return data


def filter_config(config_dict: dict, config: Union[DaijinConfiguration, MikadoConfiguration]):
    for key in ["scoring", "cds_requirements", "requirements", "not_fragmentary", "as_requirements"]:
        config_dict.pop(key, None)
    filtered_config = dict()

    # Recursive method here
    for key, value in config_dict.items():
        orig_type = type(value)
        select_attribute_for_output(filtered_config, config, key, getattr(config, key))
    clear_dict(filtered_config)
    return filtered_config


def print_toml_config(config_dict: dict, config: Union[DaijinConfiguration, MikadoConfiguration], out):

    """Function to print out the configuration in TOML format, adding the descriptions as comments preceded by #.
    :param config_dict: the configuration dictionary, optionally filtered
    :type config_dict: dict

    :param config: configuration object
    :type config: (MikadoConfiguration|DaijinConfiguration)

    :param out: output handle
    :type out: io.TextIOWrapper

    :param no_files: boolean. If on, sections pertaining to files will be deleted from the output.
    :type no_files: bool
    """

    output = toml.dumps(config_dict)

    lines = []
    level = config

    for line in output.split("\n"):
        if line.startswith("["):
            keys = line.rstrip().replace("[", "").replace("]", "").split(".")
            level = config
            for key in keys[:-1]:
                level = getattr(level, key)
            metadata = level.Schema._declared_fields[keys[-1]].metadata
            comment = []
            description = metadata.get('description', None)
            if description:
                comment += ["# " + _ for _ in textwrap.wrap(description)]
            dataclassObject = getattr(level, keys[-1])
            if hasattr(dataclassObject, "Schema"):
                level = getattr(level, keys[-1])
            lines.append(line)
            lines.extend(comment)
        else:
            comment = []
            if "=" in line:
                key = line.split("=")[0].strip()
                meta = level.Schema._declared_fields.get(key)

                description = None
                if meta:
                    description = meta.metadata.get("description")
                if description:
                    _comment = textwrap.wrap(description)
                    if _comment:
                        _comment[0] = key + ": " + _comment[0]
                        comment += ["# " + _ for _ in _comment]
                lines.extend(comment)
            # if level.Schema._declared_fields[key].required or full:
            lines.append(line.rstrip())

    if config.__doc__:
        print(*["# " + _ for _ in textwrap.wrap(config.__doc__.strip())], sep="\n", file=out)
        print("#", file=out)

    print(*lines, sep="\n", file=out)


def print_yaml_config(config_dict: dict, config: Union[DaijinConfiguration, MikadoConfiguration], out):
    """
    Function to print out the configuration in YAML format, adding the descriptions as comments preceded by #.
    :param config_dict: the configuration dictionary, optionally filtered
    :type config_dict: dict

    :param config: configuration object
    :type config: (MikadoConfiguration|DaijinConfiguration)

    :param out: output handle
    :type out: io.TextIOWrapper

    :param no_files: boolean. If on, sections pertaining to files will be deleted from the output.
    :type no_files: bool
    """

    output = yaml.dump(config_dict, default_flow_style=False)
    lines = []
    nesting = []
    comment = []

    for line in output.split("\n"):
        if line and line.startswith(line.lstrip()[0]):  # We are back at the head level
            nesting = []
        if line.endswith(":"):  # New level
            spaces = sum(1 for _ in itertools.takewhile(str.isspace, line))
            key = line.strip().rstrip(":")
            while nesting and nesting[-1][1] >= spaces:
                nesting = nesting[:-1]
            level = config
            for oldkey, _ in nesting:
                _ = getattr(level, oldkey)
                level = _
            try:
                field_dc = level.Schema._declared_fields[key]
            except KeyError:
                raise KeyError(key, level)
            except AttributeError:
                raise AttributeError(key, level)
            comment = []
            description = field_dc.metadata.get("description", None)
            if description:
                comment += ["# " + _ for _ in textwrap.wrap(description)]
            lines.append(line)
            nesting.append((key, spaces))
        else:
            spaces = sum(1 for _ in itertools.takewhile(str.isspace, line))
            comment = [" " * spaces + _ for _ in comment]
            lines.extend(comment)
            comment = []
            key = line.split(":")[0].strip()
            level = config
            if ":" in line:
                for stop, (nest_level, _) in enumerate(nesting):
                    _ = getattr(level, nest_level)
                    if not hasattr(_, "Schema"):
                        break
                    elif key in _.Schema._declared_fields:
                        level = _
                        break
                    level = _
                try:
                    meta = level.Schema._declared_fields[key].metadata
                    description = meta.get("description", None)
                except AttributeError:
                    raise AttributeError(key, level)
                except KeyError:
                    description = None

                if description:
                    _comment = textwrap.wrap(description)
                    if _comment:
                        _comment[0] = key + ": " + _comment[0]
                        comment += [" " * spaces + "# " + _ for _ in _comment]

            lines.extend(comment)
            lines.append(line.rstrip())
            comment = []

    if config.__doc__:
        print(*["# " + _ for _ in textwrap.wrap(config.__doc__.strip())], sep="\n", file=out)
        print("#", file=out)

    print(*lines, sep="\n", file=out)


def parse_list_file(cfg: Union[MikadoConfiguration,DaijinConfiguration],
                    list_file: str) -> Union[MikadoConfiguration,DaijinConfiguration]:
    configuration = {
        "pick": {
            "chimera_split": {
                "skip": []
            }
        },
        "prepare": {
            "files": {
                "gff": [],
                "labels": [],
                "strand_specific_assemblies": [],
                "source_score": {},
                "reference": [],
                "exclude_redundant": [],
                "strip_cds": [],
                "skip": []
            }
        }
    }

    if isinstance(list_file, str):
        list_file = open(list_file)

    for counter, line in enumerate(list_file):
        line = line.strip()
        if not line or line.startswith("#"):
            continue  # Skip empty lines or headers
        fields = line.rstrip().split("\t")
        try:
            gff_name, label, stranded = fields[:3]
        except (IndexError, ValueError):
            raise IndexError(
                "The input list file must contain, for each line, at a minimum file name, unique label and "
                "strandedness. Line no. {} does not contain the expected fields.".format(counter))
        if not os.path.exists(gff_name):
            raise ValueError("Non-existent file name at line {}: {}".format(counter, gff_name))
        if label in configuration["prepare"]["files"]["labels"]:
            raise ValueError("Non-unique label specified at line {}: {}".format(counter, label))
        if stranded.lower() not in ("true", "false"):
            raise ValueError("Malformed line (no. {}) for the list: {}".format(counter, line))
        if gff_name in configuration["prepare"]["files"]["gff"]:
            raise ValueError("Repeated input file at line no. {}: {}".format(counter, line))
        elif label != '' and label in configuration["prepare"]["files"]["labels"]:
            raise ValueError("Repeated label at line {}: {}".format(counter, line))
        configuration["prepare"]["files"]["gff"].append(gff_name)
        configuration["prepare"]["files"]["labels"].append(label)
        if stranded.capitalize() == "True":
            configuration["prepare"]["files"]["strand_specific_assemblies"].append(gff_name)
        elif stranded.capitalize() == "False":
            pass
        else:
            raise ValueError("Invalid strandedness value at line no. {} (only True/False accepted): {}".format(
                counter, stranded
            ))
        if len(fields) >= 4:
            score = fields[3]
            if score == "":
                score = 0
            else:
                try:
                    score = float(fields[3])
                except ValueError:
                    raise ValueError("Invalid score specified for {} at line no. {} (only int/floats accepted): {}".format(
                        gff_name, counter, fields[3]
                    ))
            configuration["prepare"]["files"]["source_score"][label] = score
        else:
            configuration["prepare"]["files"]["source_score"][label] = 0

        for arr, pos, default in [("reference", 4, False), ("exclude_redundant", 5, False),
                                  ("strip_cds", 6, False), ("skip_split", 7, False)]:
            try:
                val = fields[pos]
                if val.lower() in ("false", "true"):
                    val = eval(val.capitalize())
                else:
                    raise ValueError(
                        "Malformed line (no. {}). The last fields should be either True or False. At line "
                        "no. {}, the {} value (position {}) is instead {}.".format(counter, counter,
                                                                                   arr, pos, val))
            except IndexError:
                val = default
            if arr == "skip_split":
                configuration["pick"]["chimera_split"]["skip"].append(val)
            else:
                configuration["prepare"]["files"][arr].append(val)

    assert "exclude_redundant" in configuration["prepare"]["files"]

    cfg.prepare.files.gff = configuration["prepare"]["files"]["gff"]
    cfg.prepare.files.labels = configuration["prepare"]["files"]["labels"]
    cfg.prepare.files.strand_specific_assemblies = configuration["prepare"]["files"]["strand_specific_assemblies"]
    cfg.prepare.files.source_score = configuration["prepare"]["files"]["source_score"]
    cfg.prepare.files.reference = configuration["prepare"]["files"]["reference"]
    cfg.prepare.files.exclude_redundant = configuration["prepare"]["files"]["exclude_redundant"]
    cfg.prepare.files.strip_cds = configuration["prepare"]["files"]["strip_cds"]
    cfg.pick.chimera_split.skip = configuration["pick"]["chimera_split"]["skip"]

    return cfg
