#!/usr/bin/env python3
# coding: utf-8

"""
Script to extract features from a GTF with certain coordinates.
"""

import sys
import argparse
from mikado_lib.parsers.GTF import GTF
from mikado_lib.loci_objects.transcript import Transcript


__author__ = 'Luca Venturini'


def launch(args):

    """
    Simple launcher script.
    """

    if args.start >= args.end:
        raise ValueError("Start greater than end: {0}\t{1}".format(args.start, args.end))

    transcript = None
    with GTF(args.gtf) as gtf:
        for row in gtf:
            if row.chrom != args.chrom:
                continue
            else:
                if row.is_transcript is True:
                    if transcript is not None and \
                            transcript.start >= args.start and transcript.end <= args.end:
                        print(transcript.__str__(to_gtf=True), file=args.out)
                    if args.assume_sorted is True and row.end > args.end:
                        break
                    transcript = Transcript(row)
                else:
                    transcript.add_exon(row)

    if transcript is not None and transcript.start >= args.start and transcript.end <= args.end:
        print(transcript.__str__(to_gtf=True), file=args.out)


def awk_parser():
    """
    Parser for the command line
    :return:
    """

    parser = argparse.ArgumentParser("Quick utility to extract features from a GTF with certain coordinates.")
    parser.add_argument("--chrom", required=True)
    parser.add_argument("-as", "--assume-sorted", dest="assume_sorted",
                        default=False, action="store_true")
    parser.add_argument("--start", type=int, default=float("-inf"))
    parser.add_argument("--end", type=int, default=float("inf"))
    parser.add_argument("gtf", type=argparse.FileType())
    parser.add_argument("out", type=argparse.FileType("w"), nargs="?", default=sys.stdout)
    parser.set_defaults(func=launch)
    return parser
