#!/usr/bin/env python3
# coding: utf-8

"""
Launcher for the Mikado compare utility.
Launcher for the Mikado compare utility.
"""

import logging
import os
import sys
from logging import handlers as log_handlers
from ..utilities.log_utils import create_default_logger, formatter
import multiprocessing as mp


__author__ = 'Luca Venturini'


def setup_logger(args):

    """
    Function to setup the logger for the compare function.
    :param args:
    :param manager:
    :return:
    """

    args.log_queue = mp.Queue(-1)
    args.queue_handler = log_handlers.QueueHandler(args.log_queue)

    if args.log is not None:
        _log_folder = os.path.dirname(args.log)
        if _log_folder and not os.path.exists(_log_folder):
            os.makedirs(_log_folder)
        handle = open(args.log, mode="wt")
        handle.close()
        handler = logging.FileHandler(handle.name)
        logger = logging.getLogger("main_compare")
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.propagate = False
    else:
        logger = create_default_logger("main_compare")
        handler = logger.handlers[0]

    if args.verbose is False:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.DEBUG)

    logger.propagate = False
    log_queue_listener = log_handlers.QueueListener(args.log_queue, logger)
    log_queue_listener.propagate = False
    log_queue_listener.start()

    queue_logger = logging.getLogger("main_queue")
    for handler in queue_logger.handlers:
        queue_logger.removeHandler(handler)
    if args.verbose is False:
        queue_logger.setLevel(logging.INFO)
    else:
        queue_logger.setLevel(logging.DEBUG)
    main_queue_handler = log_handlers.QueueHandler(args.log_queue)
    queue_logger.propagate = False
    queue_logger.addHandler(main_queue_handler)

    return args, handler, logger, log_queue_listener, queue_logger


def _shutdown(args, index_name, logger, handler, queue_logger, log_queue_listener):
    queue_logger.info("Finished")
    log_queue_listener.enqueue_sentinel()
    log_queue_listener.stop()
    args.queue_handler.close()
    [_.close() for _ in logger.handlers]
    handler.flush()
    handler.close()
    args.reference.close()
    if hasattr(args.prediction, "close"):
        args.prediction.close()
    if args.no_save_index is True:
        os.remove(index_name)
    return


def compare(args):
    """
    This function performs the comparison between two different files.

    :param args: the argparse Namespace
    """

    # Flags for the parsing
    # pylint: disable=no-member
    # multiprocessing.set_start_method(method="spawn", force=True)
    # pylint: enable=no-member

    if isinstance(args.out, str):
        _out_folder = os.path.dirname(args.out)
        if _out_folder and not os.path.exists(_out_folder):
            os.makedirs(_out_folder)

    args, handler, logger, log_queue_listener, queue_logger = setup_logger(args)
    queue_logger.info("Start")
    args.commandline = " ".join(sys.argv)
    queue_logger.info("Command line: %s", args.commandline)
    logger.handlers[0].flush()

    from .reference_preparation import prepare_index
    index_name = prepare_index(args, queue_logger)
    assert os.path.exists(index_name)

    if args.index is True:
        _shutdown(args, index_name, logger, handler, queue_logger, log_queue_listener)
        return

    # Create the output folder if it does not exist.
    if os.path.dirname(args.out) and os.path.dirname(args.out) != os.path.dirname(os.path.abspath(".")):
        dirname = os.path.dirname(args.out)
        if os.path.exists(dirname):
            assert os.path.isdir(dirname)
        else:
            os.makedirs(dirname)

    try:
        if hasattr(args, "internal") and args.internal is True:
            # raise NotImplementedError()
            from .prediction_parsers import parse_self
            from .reference_preparation.gene_dict import GeneDict
            parse_self(args, GeneDict(index_name), queue_logger)
        else:
            from .prediction_parsers import parse_prediction
            parse_prediction(args, index_name, queue_logger)
    except Exception as err:
        queue_logger.exception(err)
        log_queue_listener.enqueue_sentinel()
        handler.close()
        log_queue_listener.stop()
        args.queue_handler.close()
        raise

    _shutdown(args, index_name, logger, handler, queue_logger, log_queue_listener)

    return
