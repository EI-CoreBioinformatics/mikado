"""
Module which contains all functions related to logging.
"""


import logging

__author__ = 'Luca Venturini'


def create_null_logger(instance):
    """Static method to create a default logging instance for the loci.
    The default is a null handler (no log)
    """

    formatter = logging.Formatter(
        "{asctime} - {levelname} - {lineno} - {funcName} - {processName} - {message}",
        style="{"
        )

    if isinstance(instance, str):
        name = instance
    else:
        name = instance.__name__

    logger = logging.getLogger("{0}_logger".format(name))
    handler = logging.NullHandler()
    handler.setFormatter(formatter)
    logger.setLevel(logging.WARN)
    logger.addHandler(handler)
    return logger
