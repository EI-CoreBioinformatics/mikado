#!/usr/bin/env python3
# coding: utf-8

"""
Script to extract features from a GTF with certain coordinates.
"""

import sys
import argparse
from ...parsers.blast_utils import XMLMerger, BlastOpener
import logging

__author__ = 'Luca Venturini'


def launch(args):

    """
    Simple launcher script.

    :param args: the argparse Namespace
    """

    if len(args.xml) == 1:
        parser = BlastOpener(args.xml[0])
    else:
        if args.verbose is True:
            level = logging.DEBUG
        else:
            level = logging.WARN
        parser = XMLMerger(args.xml, log=args.log, log_level=level)

    for line in parser:
        print(line, file=args.out)
    parser.close()

def merger_parser():
    """
    Parser for the command line
    :return:
    """

    parser = argparse.ArgumentParser(
        "Quick utility to extract features from a GTF with certain coordinates.")
    parser.add_argument("-v", "--verbose", action="store_true", default=False)
    parser.add_argument("-l", "--log", default=None)
    parser.add_argument("--out", type=argparse.FileType("w"), nargs="?", default=sys.stdout)
    parser.add_argument("xml", type=str, nargs="+")
    parser.set_defaults(func=launch)
    return parser
