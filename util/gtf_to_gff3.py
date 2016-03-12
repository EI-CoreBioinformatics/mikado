#!/usr/bin/env python3

import argparse
import sys
from Mikado.parsers.GTF import GTF, GtfLine
from Mikado.loci.reference_gene import Gene
from Mikado.loci.transcript import Transcript

__author__ = 'Luca Venturini'


def main():

    """
    Main utility function
    :return:
    """

    parser = argparse.ArgumentParser("Script to convert from GTF to GFF3")
    parser.add_argument("-s", "--source", default=None, type=str,
                        help="Optional new source for the GFF")
    parser.add_argument("gtf", type=GTF)
    parser.add_argument("out", default=sys.stdout, help="Output file, optional",
                        nargs="?", type=argparse.FileType("w"))
    args = parser.parse_args()

    gene_dict = dict()

    for row in args.gtf:
        if row.header is True:
            continue
        if row.attributes["gene_id"] not in gene_dict:
            gene_dict[row.attributes["gene_id"]] = dict()
        if row.attributes["transcript_id"] not in gene_dict[row.attributes["gene_id"]]:
            gene_dict[row.attributes["gene_id"]][row.attributes["transcript_id"]] = []
        if row.is_transcript is True:
            assert len(gene_dict[row.attributes["gene_id"]][row.attributes["transcript_id"]]) == 0
        gene_dict[row.attributes["gene_id"]][
            row.attributes["transcript_id"]].append(row._line)

    print("##gff-version 3", file=args.out)
    for gene_name in gene_dict:
        transcripts = []
        for tid in gene_dict[gene_name]:
            lines = [GtfLine(line) for line in gene_dict[gene_name][tid]]
            transcript_line = [_ for _ in lines if _.is_transcript is True]
            if len(transcript_line) == 1:
                transcript = Transcript(transcript_line[0])
            elif len(transcript_line) == 0:
                transcript = Transcript()
                e_line = lines[0]
                transcript.chrom = e_line.chrom
                transcript.feature = "transcript"
                transcript.strand = e_line.strand
                transcript.id = tid
                transcript.parent = gene_name
                transcript.start = min(_.start for _ in lines)
                transcript.end = max(_.end for _ in lines)
            else:
                raise ValueError("Multiple transcript lines found for {0}!".format(tid))
            transcript.add_exons([_ for _ in lines if _.is_exon is True])
            transcript.finalize()
            if args.source is not None:
                transcript.source = args.source
            transcripts.append(transcript)

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
