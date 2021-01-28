"""
This module defines the functions needed to check the sanity of the configuration file,
plus the JSON schemas for the configuration and scoring files.
"""
from typing import Union

from .configuration import MikadoConfiguration
from .daijin_configuration import DaijinConfiguration
from . import configurator
import itertools
import textwrap
try:
    import rapidjson as json
except (ImportError,ModuleNotFoundError):
    import json
import toml
import dataclasses
import yaml


__author__ = 'Luca Venturini'


def print_toml_config(config, out, no_files=False):

    """Function to print out the configuration in TOML format, adding the descriptions as comments preceded by #.
    :param config: configuration
    :type config: (MikadoConfiguration|DaijinConfiguration)

    :param out: output handle
    :type out: io.TextIOWrapper

    :param no_files: boolean. If on, sections pertaining to files will be deleted from the output.
    :type no_files: bool
    """

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

    output = toml.dumps(config_dict)

    lines = []
    level = config

    for line in output.split("\n"):
        if line.startswith("["):
            keys = line.rstrip().replace("[", "").replace("]", "").split(".")
            level = config
            for key in keys[:-1]:
                level = getattr(level, key)
            metadata = level.__dataclass_fields__[keys[-1]].metadata
            comment = []
            description = metadata.get("metadata", dict()).get("description", None)
            if description:
                comment += ["# " + _ for _ in textwrap.wrap(description)]
            if hasattr(getattr(level, keys[-1]), "__dataclass_fields__"):
                level = getattr(level, keys[-1])
            lines.append(line)
            lines.extend(comment)
        else:
            comment = []
            if "=" in line:
                key = line.split("=")[0].strip()
                try:
                    meta = level.__dataclass_fields__[key].metadata
                    description = meta.get("metadata", dict()).get("description", None)
                except AttributeError:
                    raise AttributeError(key, level)
                except KeyError:
                    description = None
                if description:
                    _comment = textwrap.wrap(description)
                    if _comment:
                        _comment[0] = key + ": " + _comment[0]
                        comment += ["# " + _ for _ in _comment]
                lines.extend(comment)
            lines.append(line.rstrip())

    if config.__doc__:
        print(*["# " + _ for _ in textwrap.wrap(config.__doc__.strip())], sep="\n", file=out)
        print("#", file=out)

    print(*lines, sep="\n", file=out)


def print_config(config: Union[MikadoConfiguration, DaijinConfiguration], out, format="yaml", no_files=False):
    """
    Function to print out the configuration in TOML format, adding the descriptions as comments preceded by #.
    :param config: configuration
    :type config: (MikadoConfiguration|DaijinConfiguration)

    :param out: output handle
    :type out: io.TextIOWrapper

    :param format: one of yaml, json or toml (case-insensitive)
    :type format: str

    :param no_files: boolean. If on, sections pertaining to files will be deleted from the output.
    :type no_files: bool
    """

    if not isinstance(format, str) or format.lower() not in ("yaml", "json", "toml"):
        raise ValueError("Unknown format: {}. I can only accept yaml, json or toml as options.")

    format = format.lower()

    if format == "toml":
        print_toml_config(config, out, no_files=no_files)
    elif format == "yaml":
        print_yaml_config(config, out, no_files=no_files)
    elif format == "json":
        config_dict = dataclasses.asdict(config)
        # Necessary otherwise we will be deleting fields from the *original* object!
        for key in ["scoring", "cds_requirements", "requirements", "not_fragmentary", "as_requirements"]:
            config_dict.pop(key, None)

        if no_files is True:
            for stage in ["pick", "prepare", "serialise"]:
                if "files" in config_dict[stage]:
                    del config_dict[stage]["files"]
            del config_dict["reference"]
            del config_dict["db_settings"]

        print(json.dumps(config_dict, indent=4, sort_keys=True), file=out)


def print_yaml_config(config, out, no_files=False):
    """
    Function to print out the configuration in YAML format, adding the descriptions as comments preceded by #.
    :param config: configuration
    :type config: (MikadoConfiguration|DaijinConfiguration)

    :param out: output handle
    :type out: io.TextIOWrapper

    :param no_files: boolean. If on, sections pertaining to files will be deleted from the output.
    :type no_files: bool
    """

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

    output = yaml.dump(config_dict, default_flow_style=False)

    # TODO currently this does not print all comments. E.g. use_diamond does not have the correct explanation.
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
                if not hasattr(_, "__dataclass_fields__"):
                    break
                level = _
            try:
                metadata = level.__dataclass_fields__[key].metadata
            except KeyError:
                raise KeyError(key, level)
            except AttributeError:
                raise AttributeError(key, level)
            comment = []
            description = metadata.get("metadata", dict()).get("description", None)
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
                    if not hasattr(_, "__dataclass_fields__"):
                        break
                    elif key in _.__dataclass_fields__:
                        level = _
                        break
                    level = _
                try:
                    meta = level.__dataclass_fields__[key].metadata
                    description = meta.get("metadata", dict()).get("description", None)
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
            comment = []
            lines.append(line.rstrip())

    if config.__doc__:
        print(*["# " + _ for _ in textwrap.wrap(config.__doc__.strip())], sep="\n", file=out)
        print("#", file=out)

    print(*lines, sep="\n", file=out)
