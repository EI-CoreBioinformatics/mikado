#!/usr/bin/env python3
# coding: utf-8

"""
Script to add a transcript feature to an input GTF.
"""

import sys
from Mikado.parsers.GTF import GTF, GtfLine
from Mikado.transcripts import Transcript
from Mikado.utilities import overlap
from copy import deepcopy
import operator
import argparse
from collections import defaultdict
from typing import List, Generator


class Obj(object):
    """ Simple container. Just a basic object."""
    pass


def create_transcript(tid: str, lines: List[GtfLine], args: argparse.Namespace) -> Generator[Transcript]:

    """"""

    chroms = defaultdict(list)
    for line in lines:
        chroms[line.chrom].append(line)

    if len(chroms) == 1:
        # Everything as it should.
        pass
    else:
        # Recursively
        for chrom in chroms:
            for transcript in create_transcript(tid + "." + chrom, chroms[chrom], args):
                yield transcript

    # Now we are sure that we only have one chromosome
    exons = sorted([line for line in lines if line.is_exon],
                   key=operator.attrgetter("chrom", "start", "end"))

    if len(exons) == 1:
        transcript = Transcript(exons[0])
        transcript.finalize()
        yield transcript
    else:
        new_exons = []
        identifier = ord("A")
        for pos in range(1, len(exons)):
            exon = exons[pos]
            prev = exons[pos - 1]
            if overlap((prev.start, prev.end), (exon.start, exon.end)) > 0:
                # Merge the two exons
                exons[pos].start = prev.start
            elif exon.start - prev.end + 1 < args.min_intron:
                if args.split is False:
                    exons[pos].start = prev.start
                else:
                    # we have to split
                    pass
            elif exon.start - prev.end + 1 > args.max_intron:
                # we have to split
                pass
            else:
                new_exons.append(prev)





def main():
    """
    Main script function.
    """

    parser = argparse.ArgumentParser("Script to add a transcript feature to e.g. Cufflinks GTFs")
    parser.add_argument("-mai", "--max-intron", dest="max_intron",
                        help="Maximum intron length before splitting a transcript into different pieces.")
    parser.add_argument("-mi", "--min-intron", dest="min_intron",
                        help="""Minimum intron length; intron lengths lower than this will cause two consecutive exons
                        to be merged.""")
    parser.add_argument("--split-small-introns", dest="split", action="store_true", default=False,
                        help="""Flag. If set, transcripts with very small introns will end up 
                        split into two (or more) transcripts rather than having their exons merged.""")
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

    for tid, lines in transcript_lines.items():
        transcripts.extend(*create_transcript(tid, lines, args))

    for transcript in sorted(transcripts):
        print(transcript.format("gtf"), file=args.out)

    if args.out is not sys.stdout:
        args.out.close()

if __name__ == '__main__':
    main()
