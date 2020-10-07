#!/usr/bin/env python3
# coding: utf-8

"""
Script to extract features from a GTF with certain coordinates.
"""

import argparse
import sys
from ...utilities import to_region


__author__ = 'Luca Venturini'


def launch(args):

    """
    Simple launcher script.

    :param args: the argparse Namespace
    """
    from ...transcripts.transcript import Transcript
    from ...parsers.GTF import GTF

    if hasattr(args, "region") and args.region is not None:
        try:
            args.chrom, args.start, args.end = args.region
        except ValueError as exc:
            raise ValueError("{0} {1}".format(exc, args.region))

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
                        print(transcript.format("gtf"), file=args.out)
                        transcript = None
                    if args.assume_sorted is True and row.start > args.end:
                        break
                    transcript = Transcript(row)
                else:
                    transcript.add_exon(row)

    if transcript is not None and transcript.start >= args.start and transcript.end <= args.end:
        print(transcript.format("gtf"), file=args.out)


def awk_parser():
    """
    Parser for the command line
    :return:
    """

    parser = argparse.ArgumentParser(
        "Utility to extract features from a GTF within given coordinates.")
    region = parser.add_mutually_exclusive_group(required=True)
    region.add_argument("-r", "--region", type=to_region,
                        help="Region defined as a string like <chrom>:<start>..<end>")
    region.add_argument("--chrom", type=str)
    parser.add_argument("-as", "--assume-sorted", dest="assume_sorted",
                        default=False, action="store_true")
    parser.add_argument("--start", type=int, default=float("-inf"))
    parser.add_argument("--end", type=int, default=float("inf"))

    parser.add_argument("gtf", type=argparse.FileType())
    parser.add_argument("out", type=argparse.FileType("w"), nargs="?", default=sys.stdout)
    parser.set_defaults(func=launch)
    return parser
