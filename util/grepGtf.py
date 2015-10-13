#!/usr/bin/env python3
# coding: utf-8

""" Quick utility to retrieve (or exclude) records from a GTF file."""

import sys
import argparse
from mikado_lib.parsers import GTF


def main():
    """
    Main script function
    """

    parser = argparse.ArgumentParser("Quick utility to retrieve (or exclude) records from a GTF file.")
    parser.add_argument('-v', dest='reverse', action="store_true", default=False,
                        help="Flag. If selected, it excludes the selected ids.")
    parser.add_argument('ids', type=argparse.FileType(), help="The file with the ids.")
    parser.add_argument('gtf', type=argparse.FileType(), help="The GTF file to analyze.")
    parser.add_argument("out", nargs="?", default=sys.stdout,
                        type=argparse.FileType("w"), help="Output file. Default: STDOUT")
    args = parser.parse_args()

    ids = set([f.rstrip() for f in args.ids])

    for record in GTF.GTF(args.gtf):
        if not record:
            continue
        if record.transcript in ids:
            if not args.reverse:
                print(record, file=args.out)
        elif args.reverse:
            print(record, file=args.out)


if __name__ == '__main__':
    main()
