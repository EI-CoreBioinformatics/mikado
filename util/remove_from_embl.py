#!/usr/bin/env python3
# coding: utf-8

""" Script to remove sequences specific of a given organism from a SwissProt file. """

import sys
import argparse
from Bio import SeqIO
import os
import re
import gzip


def main():
    """Main script function"""
    def to_file(string):
        """ Function to check the correct parser for the input file.
        :param string: filename
        """
        assert os.path.exists(string)
        if string.endswith(".gz"):
            handle = gzip.open(string)
        else:
            handle = open(string)
        return SeqIO.parse(handle, "swiss")

    parser = argparse.ArgumentParser("Script to remove sequences specific of a given organism from a SwissProt file.")
    parser.add_argument("-o", "--organism", type=str, help="Organism to be excluded", required=True)
    parser.add_argument("--format", type=str, choices=["fasta"], default="fasta",
                        help="Output format. Choices: %(choices)s. Default: %(default)s.")
    parser.add_argument("input", type=to_file)
    parser.add_argument("out", default=sys.stdout, type=argparse.FileType("w"), nargs="?")
    args = parser.parse_args()

    for record in args.input:
        if re.search(args.organism, record.annotations["organism"]) is None:
            SeqIO.write(record, args.out, args.format)
        else:
            continue


if __name__ == "__main__":
    main()
