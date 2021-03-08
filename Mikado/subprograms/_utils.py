import argparse
from typing import Union
import io
from ..configuration import MikadoConfiguration, DaijinConfiguration
import logging
import os
from ..utilities.log_utils import create_logger_from_conf


def _set_pick_mode(conf: Union[MikadoConfiguration, DaijinConfiguration], mode: str) -> Union[
        MikadoConfiguration, DaijinConfiguration]:
    if mode is not None:
        if mode == "nosplit":
            conf.pick.chimera_split.execute = False
        else:
            conf.pick.chimera_split.execute = True
            if mode == "split":
                conf.pick.chimera_split.blast_check = False
            else:
                conf.pick.chimera_split.blast_check = True
                conf.pick.chimera_split.blast_params.leniency = mode.upper()
    return conf


def check_log_settings_and_create_logger(conf: Union[MikadoConfiguration, DaijinConfiguration],
                                         args: argparse.Namespace, level=None) -> (
        Union[MikadoConfiguration, DaijinConfiguration], argparse.Namespace, logging.Logger):
    """
    Quick method to check the consistency of log settings from the namespace and create a logger.
    :param conf: the configuration object
    :param args: a Namespace
    :return: args
    """

    assert isinstance(conf, (MikadoConfiguration, DaijinConfiguration))
    out_dir = getattr(conf, level).files.output_dir if level is not None else ""
    conf = conf.copy()
    if args.log == "stderr":
        conf.log_settings.log = None
    elif args.log is not None:
        if isinstance(args.log, (io.TextIOWrapper, io.BufferedWriter)):
            args.log.close()
            conf.log_settings.log = args.log.name
        elif isinstance(args.log, bytes):
            conf.log_settings.log = args.log.decode()
        elif isinstance(args.log, str):
            conf.log_settings.log = args.log
        else:
            raise TypeError("Invalid log type: {}".format(args.log))
    elif conf.log_settings.log is None and level is not None:
        other = getattr(conf, level).files.log
        if other is not None:
            conf.log_settings.log = other

    if conf.log_settings.log is not None:
        if not out_dir or (os.path.dirname(os.path.abspath(conf.log_settings.log)) == os.path.abspath(out_dir)):
            pass
        elif os.path.dirname(os.path.abspath(conf.log_settings.log)) != os.path.abspath(out_dir):
            conf.log_settings.log = os.path.join(out_dir, conf.log_settings.log)

    if level == "pick":
        conf.pick.files.log = conf.log_settings.log
    elif level == "serialise":
        conf.serialise.files.log = conf.log_settings.log
    elif level == "prepare":
        conf.prepare.files.log = conf.log_settings.log

    conf.log_settings.log_level = args.log_level if args.log_level is not None else conf.log_settings.log_level
    logger = create_logger_from_conf(conf, name=level if level else "mikado")
    return conf, args, logger
