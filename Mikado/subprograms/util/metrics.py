#!/usr/bin/env python3
# coding: utf-8

"""Quick script to automate the generation of metrics definition from the files."""

import sys
from ...transcripts import Transcript
import argparse
import tabulate
import textwrap
from itertools import zip_longest
from inspect import getdoc


__author__ = 'Luca Venturini'


def launch(args):

    """
    Main caller.

    :param args: the argparse Namespace.
    """

    metric_names = Transcript.get_available_metrics()
    print(file=args.out)
    metrics = ["tid", "parent", "score", "external_scores"]
    metrics.extend(
        sorted(metric for metric in metric_names
               if metric not in ("tid", "parent", "score")))

    if len(args.metric) > 0:
        if not all(metric in metrics for metric in args.metric):
            print("Invalid metrics selected: {}".format(
                ", ".join(sorted(metric for metric in args.metric if metric not in metrics))))
        metrics = args.metric
    elif len(args.category) > 0:
        metrics = [metric for metric in Transcript.get_available_metrics() if
                   getattr(getattr(Transcript, metric), "category", "Descriptive") in args.category]

    rows = []

    for metric in metrics:
        docstring = getdoc(getattr(Transcript, metric))
        category = getattr(getattr(Transcript, metric), "category", "Descriptive")
        usable_raw = getattr(getattr(Transcript, metric), "usable_raw", False)
        rtype = getattr(getattr(Transcript, metric), "rtype", "str")

        if metric == "external_scores":
            usable_raw = True
            rtype = "Namespace"
            category = "External"

        if docstring is None:
            docstring = ''
        else:
            docstring = textwrap.wrap(docstring, 57)

        met_rows = zip_longest([metric], docstring, [category], [rtype], [str(usable_raw)])
        rows.extend(met_rows)

    __table_format = tabulate._table_formats[args.format]

    # TableFormat(lineabove=Line("+", "-", "+", "+"),
    #             linebelowheader=Line("+", "=", "+", "+"),
    #             linebetweenrows=Line("+", "-", "+", "+"),
    #             linebelow=Line("+", "-", "+", "+"),
    #             headerrow=DataRow("|", "|", "|"),
    #             datarow=DataRow("|", "|", "|"),
    #             padding=1, with_header_hide=None),
    # Line = namedtuple("Line", ["begin", "hline", "sep", "end"])
    #
    # DataRow = namedtuple("DataRow", ["begin", "sep", "end"])

    if args.format not in ("grid", "fancy_grid"):

        print(tabulate.tabulate(rows,
                                headers=["Metric name", "Description", "Category", "Data type", "Usable raw"],
                                tablefmt=args.format))
    else:
        out_of_header = False
        separator = None

        for row in tabulate.tabulate(rows,
                                     headers=["Metric name", "Description", "Category", "Data type", "Usable raw"],
                                     tablefmt=args.format).split("\n"):
            if row[:2] == __table_format.lineabove.begin + __table_format.lineabove.hline:
                separator = row
                if not out_of_header:
                    print(row)
                continue
            if row[:2] == __table_format.linebelowheader.begin + __table_format.linebelowheader.hline:
                out_of_header = True
                print(row, file=args.out)
                continue
            elif out_of_header is False:
                print(row, file=args.out)
            elif row[:2] == __table_format.linebetweenrows[0] + __table_format.linebetweenrows[1]:
                continue
            elif row[0] == __table_format.datarow.begin:
                if row.strip().split(__table_format.datarow.sep)[1].strip() != "":
                    print(separator, file=args.out)
                print(row, file=args.out)
        print(separator, file=args.out)
        print(file=args.out)
    return


def metric_parser():

    """
    Command line parser for the utility.
    """

    parser = argparse.ArgumentParser("Script to generate the available metrics")
    parser.add_argument("-f", "--format",
                        help="Format of the table to be printed out.",
                        choices=tabulate.tabulate_formats, default="rst")
    parser.add_argument("-o", "--out",
                        help="Optional output file",
                        type=argparse.FileType("w"), default=sys.stdout)
    parser.add_argument("-c", "--category",
                        help="Available categories to select from.",
                        default=[], nargs="+",
                        choices=sorted(set(
                            [_ for _ in [getattr(getattr(Transcript, metric), "category", "Descriptive") for metric in
                             Transcript.get_available_metrics()] if _ is not None] + ["Descriptive"])))
    parser.add_argument("metric", nargs="*")
    parser.set_defaults(func=launch)
    return parser
