#!/usr/bin/env python3
# coding: utf-8

"""
Script to add a transcript feature to an input GTF.
"""

import sys
from Mikado.parsers.GTF import GTF, GtfLine
from Mikado.transcripts import Transcript
from Mikado.utilities import overlap
from sys import maxsize
import operator
import argparse
from collections import defaultdict
from typing import List, Generator
from collections import deque


class Obj(object):
    """ Simple container. Just a basic object."""
    pass


def create_transcript(tid: str, parent: str, lines: List[GtfLine], args: argparse.Namespace):

    """"""

    chroms = defaultdict(list)
    for line in lines:
        chroms[line.chrom].append(line)

    if len(chroms) > 1:
        # Recursively
        for chrom in chroms:
            newtid = tid + "." + chrom
            newparent = parent + "." + chrom
            for transcript in create_transcript(newtid, newparent, chroms[chrom], args):
                assert transcript.id == newtid, (newtid, transcript.id)
                assert transcript.parent[0] == newparent
                yield transcript
    else:
        # Now we are sure that we only have one chromosome
        exons = sorted([line for line in lines if line.is_exon],
                       key=operator.attrgetter("chrom", "start", "end"))

        if len(exons) == 1:
            transcript = Transcript(exons[0])
            transcript.id = tid
            transcript.parent = parent
            transcript.finalize()
            yield transcript
        else:
            new_exons = deque()
            identifier = ord("A") - 1
            current = exons[0]

            for exon in exons[1:]:
                if ((overlap((exon.start, exon.end), (current.start, current.end)) > 0) or
                        (exon.start - current.end + 1 <= args.min_intron and args.split is False)):
                    # Merge the two exons
                    current.end = exon.end
                elif ((exon.start - current.end + 1 <= args.min_intron and args.split is True) or
                       exon.start - current.end + 1 > args.max_intron):
                    # TODO: split
                    new_exons.append(current)
                    transcript = Transcript(new_exons.popleft())
                    transcript.add_exons(new_exons)
                    transcript.finalize()
                    identifier += 1
                    transcript.parent = parent + "." + chr(identifier)
                    transcript.id = tid + "." + chr(identifier)
                    yield transcript
                    current = exon
                    new_exons = deque()
                else:
                    new_exons.append(current)
                    current = exon

            new_exons.append(current)
            transcript = Transcript(new_exons.popleft())
            transcript.add_exons(new_exons)

            if identifier == ord("A") - 1:
                transcript.id = tid
                transcript.parent = parent
            else:
                identifier += 1
                transcript.id = tid + "." + chr(identifier)
                transcript.parent = parent + "." + chr(identifier)

            transcript.finalize()
            yield transcript


def main():
    """
    Main script function.
    """

    parser = argparse.ArgumentParser("Script to add a transcript feature to e.g. Cufflinks GTFs")
    parser.add_argument("-mai", "--max-intron", dest="max_intron", default=maxsize, type=int,
                        help="Maximum intron length before splitting a transcript into different pieces.")
    parser.add_argument("-mi", "--min-intron", dest="min_intron", default=0, type=int,
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
        parent = lines[0].gene
        for transcript in create_transcript(tid, parent, lines, args):
            transcripts.append(transcript)

    for transcript in sorted(transcripts):
        print(transcript.format("gtf"), file=args.out)

    if args.out is not sys.stdout:
        args.out.close()


if __name__ == '__main__':
    main()
