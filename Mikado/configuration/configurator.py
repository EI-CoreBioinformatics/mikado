#!/usr/bin/env python3
# coding: utf-8

"""
This module defines the functionalities needed to verify the integrity and completeness
of Mikado configuration files. Missing values are replaced with default ones,
while existing values are checked for type and consistency.
"""
import dataclasses
import json
import os.path
import pprint
import marshmallow
from multiprocessing import get_start_method
from logging import Logger
import yaml
from pkg_resources import resource_stream
from ..exceptions import InvalidConfiguration
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
        with open(cluster_config, "wb") as out, resource_stream("Mikado.daijin", "hpc.yaml") as original:
            for pos, line in enumerate(original):
                out.write(line)
                if pos == 0:
                    out.write("    queue: {}\n".format(args.queue).encode())


def check_db(configuration: Union[MikadoConfiguration, DaijinConfiguration]) -> Union[MikadoConfiguration,
                                                                                      DaijinConfiguration]:

    """
    Function to check the validity of the database options.
    :param configuration: (MikadoConfiguration|DaijinConfiguration)
    :return: 
    """

    if configuration.db_settings.dbtype in ("mysql", "postgresql"):
        if configuration.db_settings.dbhost is None:
            raise InvalidConfiguration(
                "No host specified for the {0} database!".format(
                    configuration.db_settings.dbtype))
        if configuration.db_settings.dbuser is None:
            raise InvalidConfiguration(
                "No user specified for the {0} database!".format(
                    configuration.db_settings.dbtype))
        if configuration.db_settings.dbport == 0:
            if configuration.db_settings.dbtype == "mysql":
                configuration.db_settings.dbport = 3306
            else:
                configuration.db_settings.dbport = 5432
    else:
        if configuration.db_settings.db is not None and not os.path.isabs(configuration.db_settings.db):
            if os.path.exists(os.path.join(os.getcwd(), configuration.db_settings.db)):
                configuration.db_settings.db = os.path.join(os.getcwd(), configuration.db_settings.db)
            elif os.path.exists(
                    os.path.join(os.path.dirname(configuration.filename or ""), configuration.db_settings.db)):
                configuration.db_settings.db = os.path.join(
                    os.path.dirname(configuration.filename or ""), configuration.db_settings.db)
            elif os.path.exists(os.path.join(os.path.dirname(configuration.filename or ""),
                                             configuration.db_settings.db)):
                configuration.db_settings.db = os.path.join(os.path.dirname(configuration.filename or ""),
                                                            configuration.db_settings.db)
            else:
                pass

    return configuration


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
        configuration.load_scoring(logger=logger)
        configuration.check(logger=logger)
        configuration = check_db(configuration)
        if not configuration.multiprocessing_method:
            configuration.multiprocessing_method = get_start_method()
    except InvalidConfiguration as exc:
        logger.exception(exc)
        raise

    assert configuration.seed is not None
    random.seed(configuration.seed % (2 ** 32 - 1))
    return configuration


def load_and_validate_config(raw_configuration: Union[None, MikadoConfiguration, DaijinConfiguration,
                             str, dict], logger=None, external=False) -> Union[MikadoConfiguration, DaijinConfiguration]:
    """
    Function to serialise the JSON for configuration and check its consistency.

    :param raw_configuration: either the file name of the configuration or an initialised object to check and finalise.
    :type raw_configuration: (str | None | dict | MikadoConfiguration | DaijinConfiguration)

    :param external: boolean. If True, presume the file might be an old/incomplete configuration to pass in values.
    Therefore, accept also *partial* configuration files.
    :type external: bool

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
                config = MikadoConfiguration.Schema().load(raw_configuration, partial=external)
            except marshmallow.exceptions.ValidationError:
                config = DaijinConfiguration.Schema().load(raw_configuration, partial=external)
        elif raw_configuration is None or raw_configuration == '' or raw_configuration == dict():
            config = MikadoConfiguration()
        else:
            assert isinstance(raw_configuration, str), raw_configuration
            raw_configuration = os.path.abspath(raw_configuration)
            if not os.path.exists(raw_configuration) or os.stat(raw_configuration).st_size == 0:
                raise InvalidConfiguration("JSON file {} not found!".format(raw_configuration))
            with open(raw_configuration) as json_file:
                # YAML *might* try to load up the file as the proper object
                if raw_configuration.endswith(".yaml"):
                    config = yaml.load(json_file, Loader=yaml.Loader)
                elif raw_configuration.endswith(".toml"):
                    config = toml.load(json_file)
                elif raw_configuration.endswith(".json"):
                    config = json.loads(json_file.read())
                else:
                    config = toml.load(json_file)
            assert isinstance(config, dict), (config, type(config))
            config["filename"] = raw_configuration
            try:
                config = MikadoConfiguration.Schema().load(config, partial=external)
            except marshmallow.exceptions.ValidationError as mikado_exc:
                try:
                    config = DaijinConfiguration.Schema().load(config, partial=external)
                except marshmallow.exceptions.ValidationError as exc:
                    logger.critical("The configuration file is invalid. Validation errors:\n%s\n%s\n\n",
                                    pprint.pformat(mikado_exc.messages),
                                    pprint.pformat(exc.messages))
                    logger.critical(exc)
                    raise

            config.filename = raw_configuration
            if config.scoring is None:
                check_and_load_scoring(config, logger=logger)
            config.scoring.check(config.pick.orf_loading.minimal_orf_length)

        assert isinstance(config, (MikadoConfiguration, DaijinConfiguration)), type(config)
        config = check_and_load_scoring(config, logger=logger)
    except KeyboardInterrupt:
        raise
    except Exception as exc:
        logger.exception("Loading the configuration file failed with error:\n%s\n\n\n", exc)
        raise InvalidConfiguration("The configuration file passed is invalid. Please double check.")

    random.seed(config.seed % (2 ** 32 - 1))

    return config
