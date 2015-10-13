#!/usr/bin/env python3
# coding: utf-8

"""
Script to add a transcript feature to an input GTF.
"""

import sys
from Mikado.parsers.GTF import GTF
from copy import deepcopy
import operator
import argparse


class Obj(object):
    """ Simple container. Just a basic object."""
    pass


def main():
    """
    Main script function.
    """
    current = Obj()
    current.transcript = None

    rows = []

    parser = argparse.ArgumentParser("Script to add a transcript feature to e.g. Cufflinks GTFs")
    parser.add_argument("gtf", type=argparse.FileType(),
                        help="Input GTF")
    parser.add_argument("out", default=sys.stdout, nargs="?",
                        type=argparse.FileType("w"),
                        help="Output file. Default: stdout.")
    args = parser.parse_args()

    args.gtf.close()

    for record in GTF(args.gtf.name):
        if current.transcript != record.transcript:
            if current.transcript is not None:
                print(current, file=args.out)
                exon_no = 0
                for row in filter(lambda x: x.feature == "exon",
                                  sorted(rows, key=operator.attrgetter("start"))):
                    exon_no += 1
                    row.attributes["exon_number"] = exon_no
                    print(row, file=args.out)
                exon_no = 0
                for row in filter(lambda x: x.feature == "CDS",
                                  sorted(rows, key=operator.attrgetter("start"))):
                    exon_no += 1
                    row.attributes["exon_number"] = exon_no
                    print(row, file=args.out)
            rows = [record]
            current = deepcopy(record)
            current.feature = "transcript"

        else:
            current.end = max(current.end, record.end)
            current.start = min(current.start, record.start)
            rows.append(record)

    print(current, file=args.out)
    for row in rows:
        print(row, file=args.out)

if __name__ == '__main__':
    main()
