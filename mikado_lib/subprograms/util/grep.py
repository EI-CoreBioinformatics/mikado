#!/usr/bin/env python3
# coding: utf-8

"""Script to parse and retrieve given features from a GTF/GFF file."""

import sys
import argparse
import functools
from mikado_lib.parsers import GFF, GTF
from mikado_lib.subprograms import to_gff

__author__ = 'Luca Venturini'


def print_gff_gene(curr_gene, curr_transcripts, args):
    """
    :param curr_gene: gene record
    :param curr_transcripts: dictionary of transcripts
    :param args: argparse
    :return: None
    """

    if curr_gene is not None and len(curr_transcripts) > 0:
        starts, ends = [], []
        lines = []
        for tid in curr_transcripts:
            lines.append(str(curr_transcripts[tid][0]))
            starts.append(curr_transcripts[tid][0].start)
            ends.append(curr_transcripts[tid][0].end)
            for rec in curr_transcripts[tid][1:]:
                rec.parent = [tid]
                lines.append(str(rec))
        curr_gene.start = min(starts)
        curr_gene.end = max(ends)
        print(curr_gene, file=args.out)
        print(*lines, sep="\n", file=args.out)
        print("###", file=args.out)


def verify_storability(record, mrna_ids, gene_ids, args):
    """
    This function verifies whether a GFF transcript has to be kept
    or not.
    :param record: the record to be evaluated
    :param gene_ids: a set of gene IDs
    :param mrna_ids: a set of mRNA IDs
    :param args: the namespace args
    :return:
    """

    bool_flag = False

    if args.reverse is False:
        if record.id in mrna_ids:
            bool_flag = True
        elif args.genes is True and any([p in gene_ids for p in record.parent]):
            bool_flag = True
            mrna_ids.add(record.id)
    elif args.reverse is True:
        if args.genes is False and record.id not in mrna_ids:
            bool_flag = True
        elif args.genes is True and not any([p in gene_ids for p in record.parent]):
            bool_flag = True
            mrna_ids.add(record.id)

    return bool_flag, mrna_ids


def grep_gff(args, gene_ids, mrna_ids):
    """
    Grep-like main function for *GFF3* files.
    :param args:
    :return:
    """

    curr_gene = None
    curr_transcripts = dict()

    print("##gff-version 3", file=args.out)

    evaluator = functools.partial(verify_storability,
                                  **{"gene_ids": gene_ids,
                                     "args": args})

    for record in args.gff:
        if record.is_transcript is True:
            evaluated, mrna_ids = evaluator(record, mrna_ids)
            if evaluated is True:
                curr_transcripts[record.id] = [record]
        elif record.is_exon is True:
            for parent in record.parent:
                if parent in mrna_ids and args.reverse is False:
                    curr_transcripts[parent].append(record)
                elif parent not in mrna_ids and args.reverse is True:
                    curr_transcripts[parent].append(record)
        elif record.is_derived is True:
            for derivation in record.derived_from:
                if args.reverse is True and derivation not in mrna_ids:
                    curr_transcripts[derivation].append(record)
                elif args.reverse is False and derivation in mrna_ids:
                    curr_transcripts[derivation].append(record)
        else:
            print_gff_gene(curr_gene, curr_transcripts, args)
            curr_gene = None
            curr_transcripts = dict()

            if record.id is None:
                continue
            curr_gene = record

    print_gff_gene(curr_gene, curr_transcripts, args)


def grep_gtf(args, gene_ids, mrna_ids):
    """
    Grep-like main function for *GTF* files.
    :param args:
    :return:
    """

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
    """
    Function which chooses whether to use the grep_gtf
    or the grep_gff function for the analysis.
    :param args:
    :return:
    """

    gene_ids, mrna_ids = set(), set()
    for line in args.ids:
        if args.genes is False:
            mrna_id, gene_id = line.rstrip().split()[:2]
            mrna_ids.add(mrna_id)
            gene_ids.add(gene_id)
        else:
            gene_id = line.rstrip()
            gene_ids.add(gene_id)

    if isinstance(args.gff, GTF.GTF):
        grep_gtf(args, gene_ids, mrna_ids)
    elif isinstance(args.gff, GFF.GFF3):
        grep_gff(args, gene_ids, mrna_ids)
    else:
        raise TypeError(type(args.gff))


def grep_parser():
    """
    Simple command line parser for the utility.
    :return: args
    """

    parser = argparse.ArgumentParser(
        'Script to parse and retrieve given features from a GFF file.')
    parser.add_argument('-v', action='store_true', dest='reverse',
                        help="Exclude from the gff all the records in the id file.")
    parser.add_argument("--genes", action="store_true",
                        help="""Flag. If set, the program expects as ids
                        only a list of genes, and will exclude/include all the transcripts
                        children of the selected genes.""")
    parser.add_argument('ids', type=argparse.FileType(),
                        help="ID file (format: mrna_id, gene_id - tab separated)")
    parser.add_argument('gff', type=to_gff, help="The GFF file to parse.")
    parser.add_argument('out', nargs='?',
                        type=argparse.FileType('w'), default=sys.stdout,
                        help="Optional output file")
    parser.set_defaults(func=launch)
    return parser
