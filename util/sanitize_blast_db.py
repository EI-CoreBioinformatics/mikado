#!/usr/bin/env python3

import Bio.SeqIO
import sys
import argparse
from collections import Counter
from Mikado.utilities.log_utils import create_default_logger
import re
import gzip


"""Simple utility to discard any non-description information out of the FASTA file, so to reduce
the chances of invalid characters creeping in. It will also check the consistency of
the identifiers and remove duplicated ones."""


def pos(string: str):
    if not string.isdigit():
        raise ValueError("I can only accept numbers")
    string = abs(round(float(string), 0))
    return string


def main():

    logger = create_default_logger("sanitizer")

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("-ml", "--min-length", default=0, type=pos)
    parser.add_argument("-o", "--out", default=sys.stdout, type=argparse.FileType("wt"))
    parser.add_argument("fasta", nargs="+", type=argparse.FileType("rt"))
    args = parser.parse_args()

    found_ids = Counter()

    starter = 96

    for fasta in args.fasta:
        if fasta.name.endswith(".gz"):
            fasta.close()
            fasta = gzip.open(fasta.name, "rt")

        if len(args.fasta) > 1:
            starter += 1
            prefix = "{}_".format(chr(starter))
        else:
            prefix = ""
        for record in Bio.SeqIO.parse(fasta, "fasta"):
            if record.id in found_ids:
                logger.warning("ID found other {} time{} in the input files!".format(
                    found_ids[record.id], "s" if found_ids[record.id] > 1 else ""))
            record.id = "{}{}".format(prefix, re.sub(r"\|", "_", record.id))
            record.description = ""
            if len(record) < args.min_length:
                continue
            Bio.SeqIO.write(record, args.out, "fasta")

    args.out.close()


main()
