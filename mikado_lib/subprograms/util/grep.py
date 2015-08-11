#!/usr/bin/env python3
# coding: utf-8

"""Script to parse and retrieve given features from a GTF/GFF file."""

import sys
import argparse
from mikado_lib.parsers import GFF, GTF
from mikado_lib.subprograms import to_gff

__author__ = 'Luca Venturini'


def grep_gff(args):

    gene_ids, mrna_ids = set(), set()
    for line in args.ids:
        if args.genes is False:
            mrna_id, gene_id = line.rstrip().split()[:2]
            mrna_ids.add(mrna_id)
            gene_ids.add(gene_id)
        else:
            gene_id = line.rstrip()
            gene_ids.add(gene_id)

    curr_gene = None
    curr_transcripts = dict()

    print("##gff-version 3", file=args.out)

    for record in args.gff:
        if record.is_transcript is True:  # Potential gene line
            if args.reverse is False and (
                    record.id in mrna_ids or (args.genes is True and any([p in gene_ids for p in record.parent]))):
                curr_transcripts[record.id] = [record]
            elif args.reverse is True and ((args.genes is False and record.id not in mrna_ids) or (
                    args.genes is True and not any([p in gene_ids for p in record.parent]))):
                curr_transcripts[record.id] = [record]
        elif record.is_exon is True:
            for parent in record.parent:
                if parent in mrna_ids:
                    curr_transcripts[parent].append(record)
        else:
            if curr_gene is not None and len(curr_transcripts) > 0:
                print(curr_gene, file=args.out)
                for tid in curr_transcripts:
                    print(curr_transcripts[tid][0], file=args.out)
                    for rec in curr_transcripts[tid][1:]:
                        rec.parent = [tid]
                        print(rec, file=args.out)

            curr_gene = None
            curr_transcripts = dict()

            if record.id is None:
                continue
            elif args.reverse is True and record.id not in gene_ids:
                curr_gene = record
            elif args.reverse is False and record.id in gene_ids:
                curr_gene = record

    if curr_gene is not None and len(curr_transcripts) > 0:
        print(curr_gene, file=args.out)
        for tid in curr_transcripts:
            print(curr_transcripts[tid][0], file=args.out)
            for rec in curr_transcripts[tid][1:]:
                rec.parent = [tid]
                print(rec, file=args.out)


def grep_gtf(args):
    gene_ids, mrna_ids = set(), set()
    for line in args.ids:
        if args.genes is False:
            mrna_id, gene_id = line.rstrip().split()[:2]
            mrna_ids.add(mrna_id)
            gene_ids.add(gene_id)
        else:
            gene_id = line.rstrip()
            gene_ids.add(gene_id)

    for record in args.gtf:
        if not record:
            continue
        if args.genes is False:
            if record.transcript in mrna_ids:
                if not args.reverse:
                    print(record, file=args.out)
                elif args.reverse:
                    print(record, file=args.out)
        else:
            if record.gene in gene_ids:
                if not args.reverse:
                    print(record, file=args.out)
                elif args.reverse:
                    print(record, file=args.out)


def launch(args):
    if type(args.gff) is GTF.GTF:
        grep_gtf(args)
    else:
        grep_gff(args)


def grep_parser():
    parser = argparse.ArgumentParser('Script to parse and retrieve given features from a GFF file.')
    parser.add_argument('-v', action='store_true', dest='reverse',
                        help="Exclude from the gff all the records in the id file.")
    parser.add_argument("--genes", action="store_true",
                        help="""Flag. If set, the program expects as ids only a list of genes,
                        and will exclude/include all the transcripts children of the selected genes.""")
    parser.add_argument('ids', type=argparse.FileType(), help="ID file (format: mrna_id, gene_id - tab separated)")
    parser.add_argument('gff', type=to_gff, help="The GFF file to parse.")
    parser.add_argument('out', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help="Optional output file")
    parser.set_defaults(func=launch)
    return parser
