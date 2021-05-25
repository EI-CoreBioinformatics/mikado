#!/usr/bin/env python3

import sys
import argparse
from collections import OrderedDict


__doc__ = """Script to remove UTRs from a GFF/GTF."""

from Mikado.parsers import parser_factory
from Mikado.transcripts import Transcript, Gene


def strip_utr(gene: Gene):

    for transcript in gene:
        assert isinstance(transcript, Transcript)
        if len(transcript.combined_cds) > 0:
            transcript.exons = transcript.combined_cds.copy()
            transcript.start = min([_[0] for _ in transcript.combined_cds])
            transcript.end = max([_[1] for _ in transcript.combined_cds])
            transcript.combined_utr = []
        transcript.finalize()

    gene.finalize()
    return gene


def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("-f", "--format", default=None, choices=["gff3", "gtf"])
    parser.add_argument("gff", type=parser_factory)
    parser.add_argument("out", nargs="?", default=sys.stdout, type=argparse.FileType("wt"))
    args = parser.parse_args()

    is_gff = (args.gff.file_format == "gff3")
    if args.format is None:
        args.format = args.gff.file_format

    tid2gid = dict()
    genes = OrderedDict()

    for row in args.gff:
        if row.header is True:
            continue
        elif row.is_gene is True:
            genes[row.id] = Gene(row)
        elif row.is_transcript is True:
            assert len(row.parent) == 1
            parent = row.parent[0]
            tid2gid[row.id] = parent
            genes[parent].add(Transcript(row))
        elif row.is_exon is True:
            if row.gene is None:
                gene = tid2gid[row.parent[0]]
            else:
                gene = row.gene
            genes[gene].add_exon(row)

    for gid, gene in genes.items():
        print(strip_utr(gene).format(args.format), file=args.out)
        continue

    return


main()
