#!/usr/bin/env python3
# coding: utf-8

"""Quick script to automate the generation of metrics definition from the files."""
import sys

from mikado_lib.loci_objects.transcript import Transcript
import re
import argparse

__author__ = 'Luca Venturini'


def launch(args):

    """
    Main caller.
    """

    metrics = Transcript.get_available_metrics()
    print(file=args.out)
    for metric in ["tid", "parent", "score"]+sorted(filter(lambda x: x not in ("tid", "parent", "score", metrics), metrics)):
        docstring=getattr(Transcript,metric).__doc__
        if docstring is None:
            docstring=''
        print( "- *{0}*:".format(metric), re.sub(" +", " ", re.sub("\n", " ", docstring)),
               sep="\t",
               file=args.out)

    print(file=args.out)


def metric_parser():

    parser = argparse.ArgumentParser("Script to generate the available metrics")
    parser.add_argument("-o", "--out", type=argparse.FileType("w"), default=sys.stdout)
    parser.set_defaults(func=launch)
    return parser