#!/usr/bin/env python3
# coding: utf-8

"""
Script to add a transcript feature to an input GTF.
"""

import sys
from Mikado.parsers.GTF import GTF
from Mikado.transcripts import Transcript
from copy import deepcopy
import operator
import argparse
from collections import defaultdict


class Obj(object):
    """ Simple container. Just a basic object."""
    pass


def main():
    """
    Main script function.
    """

    parser = argparse.ArgumentParser("Script to add a transcript feature to e.g. Cufflinks GTFs")
    parser.add_argument("gtf", type=argparse.FileType(),
                        help="Input GTF")
    parser.add_argument("out", default=sys.stdout, nargs="?",
                        type=argparse.FileType("w"),
                        help="Output file. Default: stdout.")
    args = parser.parse_args()

    args.gtf.close()

    transcript_lines = defaultdict(list)

    [transcript_lines[_.transcript].append(_) for _ in GTF(args.gtf.name) if _.header is False and _.is_exon is True]
    args.gtf.close()
    transcripts = list()

    for tid in transcript_lines:
        transcript = Transcript(transcript_lines[tid][0])
        transcript.add_exons(transcript_lines[tid])
        transcripts.append(transcript)

    for transcript in sorted(transcripts):
        print(transcript.format("gtf"), file=args.out)

    if args.out is not sys.stdout:
        args.out.close()

if __name__ == '__main__':
    main()
