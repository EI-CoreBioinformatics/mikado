#!/usr/bin/env python3
# coding: utf-8

""" Script to calculate statistics about an annotation file.
It can take both GTF and GFF files as input."""

import sys
import argparse


__author__ = "Luca Venturini"


def launch(args):

    """
    Very simple launcher function, calls Calculator from this module.

    :param args: the argparse Namespace.
    """

    from ...scales.calculator import Calculator
    from ...parsers import to_gff
    args.gff = to_gff(args.gff)
    calculator = Calculator(args)
    calculator()


def stats_parser():

    """
    Argument parser.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('--only-coding', dest="only_coding", action="store_true", default=False)
    parser.add_argument('--tab-stats',
                        dest="tab_stats",
                        default=None,
                        type=argparse.FileType('w'),
                        help="Optional tabular file to write statistics for each transcript.")
    parser.add_argument("-v", "--verbose", action="store_true", default=False)
    parser.add_argument('gff', type=str, help="GFF file to parse.")
    parser.add_argument('out', type=argparse.FileType('w'), default=sys.stdout, nargs='?')
    parser.set_defaults(func=launch)
    return parser
