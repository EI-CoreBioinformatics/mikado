#!/usr/bin/env python3

import argparse
import sys
import mikado.parsers.GTF
from mikado.loci_objects.reference_gene import Gene
from mikado.loci_objects.transcript import Transcript
from collections import defaultdict

__author__ = 'Luca Venturini'


def main():

    """
    Main utility function
    :return:
    """

    parser = argparse.ArgumentParser("Script to convert from GTF to GFF3")
    parser.add_argument("-s", "--source", default=None, type=str,
                        help="Optional new source for the GFF")
    parser.add_argument("gtf", type=mikado.parsers.GTF.GTF)
    parser.add_argument("out", default=sys.stdout, help="Output file, optional",
                        nargs="?", type=argparse.FileType("w"))
    args = parser.parse_args()

    gene_dict = defaultdict(dict)

    for row in args.gtf:
        if row.header is True:
            continue
        if row.is_transcript is True:
            gene_dict[
                row.attributes["gene_id"]][
                row.attributes["transcript_id"]] = Transcript(row)
            if args.source is not None:
                gene_dict[
                    row.attributes["gene_id"]][
                    row.attributes["transcript_id"]].source = args.source
        else:
            gene_dict[
                row.attributes["gene_id"]][
                row.attributes["transcript_id"]].add_exon(row)

    print("##gff-version 3", file=args.out)
    for gene_name in gene_dict:
        transcripts = list(gene_dict[gene_name].values())
        gene = Gene(transcripts[0])
        gene.id = gene_name
        if args.source is not None:
            gene.source = args.source
        if len(transcripts) > 1:
            for transcript in transcripts[1:]:
                gene.add(transcript)
        print(gene.format("gff3"), file=args.out)
        print("###", file=args.out)

if __name__ == '__main__':
    main()
