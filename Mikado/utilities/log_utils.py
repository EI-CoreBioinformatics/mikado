"""
Module which contains all functions related to logging.
"""


import logging
import logging.handlers

__author__ = 'Luca Venturini'


formatter = logging.Formatter(
        "{asctime} - {name} - {filename}:{lineno} - {levelname} - {funcName} \
- {processName} - {message}",
        style="{"
        )


null_logger = logging.getLogger("null")
null_handler = logging.NullHandler()
null_handler.setFormatter(formatter)
null_logger.setLevel(logging.CRITICAL)
null_logger.addHandler(null_handler)


def create_null_logger(*args, **kwargs):
    """Static method to create a default logging instance for the loci.
    The default is a null handler (no log).

    :param instance: the instance used to derive a name for the logger. It must be either a string
    or a class instance with a __name__ attribute."""

    if len(args) > 0:
        null_special_logger = logging.getLogger(args[0])
        null_special_handler = logging.NullHandler()
        null_special_handler.setFormatter(formatter)
        null_special_logger.handlers = [null_special_handler]
        if "level" in kwargs:
            null_special_logger.setLevel(kwargs["level"])
        else:
            null_special_logger.setLevel(logging.CRITICAL)
        return null_special_logger

    return null_logger


def check_logger(logger):
    """Quick function to verify that a logger is really a logger,
    otherwise it raises a ValueError.

    :param logger: the logger instance
    :type logger: logging.Logger
    """

    if isinstance(logger, logging.Logger):
        return logger
    else:
        raise ValueError("{0} is not a logger but rather {1}".format(
            logger, type(logger)
        ))


def create_default_logger(name, level="WARN"):
    """Default logger
    :param name: string used to give a name to the logger.
    :type name: str

    :param level: level of the logger. Default: WARN
    """

    logger = logging.getLogger(name)
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logger.setLevel(level)
    logger.addHandler(handler)
    logger.propagate = False
    return logger


def create_queue_logger(instance, prefix=""):
    """
    Create a queue logger for a specific class, which *must* have
     a "logging_queue" property redirecting to a queue-like object.
     If the instance possesses a "log_level" attribute, the log level
     will be set to its value; otherwise, the logger will be configured
     with a default level of "WARNING".

    :param instance:
    :param prefix:
    :return:
    """

    if not hasattr(instance, "logging_queue"):
        raise AttributeError(
            "A queue logger can be built only using a class with the \"logging_queue\" property!")
    instance._log_handler = logging.handlers.QueueHandler(instance.logging_queue)
    if hasattr(instance, "log_level"):
        instance._log_handler.setLevel(instance.log_level)
    else:
        instance._log_handler.setLevel("WARNING")

    if prefix:
        instance.logger = logging.getLogger("{0}.{1}".format(prefix,
                                                             instance.name if
                                                             hasattr(instance, "name") else
                                                             "default"))
    else:
        instance.logger = logging.getLogger(instance.name)
    instance.logger.addHandler(instance._log_handler)
    instance.logger.setLevel(instance._log_handler.level)
    instance.logger.propagate = False
    return


def create_logger_from_conf(conf, name="mikado", mode="a"):

    logger = logging.getLogger(name)
    handle = conf["log_settings"].get("log", None)
    if handle is None:
        handler = logging.StreamHandler()
    else:
        handler = logging.FileHandler(conf["log_settings"]['log'], mode=mode)

    handler.setFormatter(formatter)
    logger.setLevel(conf["log_settings"]['log_level'])
    logger.addHandler(handler)
    logger.propagate = False
    return logger
