#!/usr/bin/env python3
# coding: utf-8

"""Simple script to serialise ORF BED files into the SQLite DB"""

import argparse
from mikado_lib.serializers import orf
from mikado_lib import json_utils
from Bio import SeqIO


def main():
    """Main script function"""
    parser = argparse.ArgumentParser("Simple script to serialise ORF BED files into the SQLite DB.")
    parser.add_argument("--fasta", default=None, required=True)
    parser.add_argument("-mo", "--max-objects", dest="max_objects", type=int, default=10 ** 5)
    parser.add_argument("--json-conf", default=None, dest="json_conf", type=json_utils.to_json)
    parser.add_argument("bed12")
    parser.add_argument("db", nargs="?", default=":memory:")
    args = parser.parse_args()
    if args.fasta is not None:
        args.fasta = SeqIO.index(args.fasta, "fasta")

    serializer = orf.OrfSerializer(args.bed12, args.db, fasta_index=args.fasta, maxobjects=args.max_objects,
                                   json_conf=args.json_conf)
    serializer.serialize()


if __name__ == "__main__":
    main()
