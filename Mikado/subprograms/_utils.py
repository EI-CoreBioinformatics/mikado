import argparse
import dataclasses
from typing import Union
import io
from ..configuration import MikadoConfiguration, DaijinConfiguration
import logging
import os

from ..exceptions import InvalidConfiguration
from ..utilities.log_utils import create_logger_from_conf


def _set_pick_mode(conf: Union[MikadoConfiguration, DaijinConfiguration], mode: str) -> Union[
        MikadoConfiguration, DaijinConfiguration]:
    if mode is not None:
        mode = mode.lower()
        if mode == "nosplit":
            conf.pick.chimera_split.execute = False
            conf.pick.chimera_split.blast_check = False
        else:
            conf.pick.chimera_split.execute = True
            if mode == "split":
                conf.pick.chimera_split.blast_check = False
            else:
                conf.pick.chimera_split.blast_check = True
                conf.pick.chimera_split.blast_params.leniency = mode.upper()

    errors = conf.pick.chimera_split.Schema().validate(dataclasses.asdict(conf.pick.chimera_split))
    if len(errors) > 0:  # Dict of errors
        raise InvalidConfiguration(f"Invalid mode selected for Mikado: {mode}\nErrors: {errors}\nDict: "
                                   f"{dataclasses.asdict(conf.pick.chimera_split)}")

    return conf


def check_log_settings_and_create_logger(conf: Union[MikadoConfiguration, DaijinConfiguration],
                                         log: Union[str, io.TextIOWrapper, io.BufferedWriter, bytes, None],
                                         log_level: Union[str, None],
                                         section=None) -> (
        Union[MikadoConfiguration, DaijinConfiguration], logging.Logger):
    """
    Quick method to check the consistency of log settings from the namespace and create a logger.
    :param conf: the configuration object
    :param log: Handle to use for the logger. Use "stderr" for using a streaming logger. If log is None:
    - if the section (see below) if also None, use stderr
    - else, try to infer the log from the parameter in the section
    :param log_level: level of the log to use
    :param section: section of the configuration to use, if any. One between None, pick, serialise, prepare
    :return: MikadoConfiguration, logger
    """

    if not isinstance(conf, MikadoConfiguration):
        raise InvalidConfiguration(f"Invalid object type; it should be MikadoConfiguration, it is {type(conf)}")
    out_dir = getattr(conf, section).files.output_dir if section is not None else ""
    # Set to standard error
    if log == "stderr" or (log is None and section is None):
        conf.log_settings.log = None
        if section == "pick":
            conf.pick.files.log = None
        elif section == "serialise":
            conf.serialise.files.log = None
        elif section == "prepare":
            conf.serialise.files.log = None
    # Set to the specified log output
    elif log is not None:
        if isinstance(log, (io.TextIOWrapper, io.BufferedWriter)):
            log.close()
            conf.log_settings.log = log.name
        elif isinstance(log, bytes):
            conf.log_settings.log = log.decode()
        elif isinstance(log, str):
            conf.log_settings.log = log
        else:
            raise TypeError("Invalid log type: {}".format(log))
    # infer from the section log, if present
    elif log is None and section is not None:
        other = getattr(conf, section).files.log
        if other is not None:
            conf.log_settings.log = other

    if conf.log_settings.log is not None:
        if not out_dir or (os.path.dirname(os.path.abspath(conf.log_settings.log)) == os.path.abspath(out_dir)):
            pass
        elif os.path.dirname(conf.log_settings.log) != "":
            conf.log_settings.log = os.path.abspath(conf.log_settings.log)
        elif os.path.dirname(os.path.abspath(conf.log_settings.log)) != os.path.abspath(out_dir):
            conf.log_settings.log = os.path.join(out_dir, conf.log_settings.log)

    if section == "pick":
        conf.pick.files.log = conf.log_settings.log
    elif section == "serialise":
        conf.serialise.files.log = conf.log_settings.log
    elif section == "prepare":
        conf.prepare.files.log = conf.log_settings.log

    try:
        conf.log_settings.log_level = log_level.upper() if log_level is not None else conf.log_settings.log_level
    except (TypeError, AttributeError):
        raise TypeError(f"Invalid log_level: {log_level}, type {type(log_level)}")
    # Check the log section
    errors = conf.log_settings.Schema().validate(dataclasses.asdict(conf.log_settings))
    if len(errors) > 0:
        raise InvalidConfiguration(
            f"Invalid logging selected for Mikado: log {log}, level {log_level}, section {section}\n"
            f"Errors: {errors}\nDict: {dataclasses.asdict(conf.pick.chimera_split)}")

    logger = create_logger_from_conf(conf, name=section if section else "mikado")
    return conf, logger
