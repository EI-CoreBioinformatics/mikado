#!/usr/bin/env python3
# coding: utf-8

"""Quick script to automate the generation of metrics definition from the files."""

import sys
from ...loci import Transcript
import re
import argparse
import tabulate
import textwrap
from itertools import zip_longest


__author__ = 'Luca Venturini'


def launch(args):

    """
    Main caller.

    :param args: the argparse Namespace.
    """

    metric_names = Transcript.get_available_metrics()
    print(file=args.out)
    metrics = ["tid", "parent", "score"]
    metrics.extend(
        sorted(metric for metric in metric_names
               if metric not in ("tid", "parent", "score")))

    rows = []

    for metric in metrics:
        docstring = getattr(Transcript, metric).__doc__
        if not hasattr(getattr(Transcript, metric), "category"):
            category = "Descriptive"  # Hack for tid etc
        else:
            category = getattr(Transcript, metric).category
        if docstring is None:
            docstring = ''
        else:
            docstring = textwrap.wrap(docstring, 57)

        met_rows = zip_longest([metric], docstring, [category])
        rows.extend(met_rows)
        # print("- *{0}*:".format(metric), re.sub(" +", " ", re.sub("\n", " ", docstring)),
        #       sep="\t", file=args.out)

    separator = None
    out_of_header = False
    for row in tabulate.tabulate(rows,
                            headers=["Metric name", "Description", "Data type", "Category"],
                            tablefmt="grid").split("\n"):
        if row[:2] == "+-":
            separator = row
            if not out_of_header:
                print(row)
            continue
        if row[:2] == "+=":
            out_of_header = True
            print(row, file=args.out)
            continue
        elif out_of_header is False:
            print(row, file=args.out)
            continue
        elif row[:2] != "+-":
            if row.rstrip().split("|")[1].strip() != "":
                print(separator, file=args.out)
            print(row, file=args.out)
    print(separator, file=args.out)

    # print(file=args.out)


def metric_parser():

    """
    Command line parser for the utility.
    """

    parser = argparse.ArgumentParser("Script to generate the available metrics")
    parser.add_argument("-o", "--out", type=argparse.FileType("w"), default=sys.stdout)
    parser.set_defaults(func=launch)
    return parser
