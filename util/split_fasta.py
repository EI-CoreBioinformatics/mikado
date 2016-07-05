#!/usr/bin/env python3

import sys
import argparse
import os
import textwrap
import pyfaidx
from math import log, ceil


def positive(string):

    string = int(string)
    if string < 1:
        raise ValueError("Only positive values are acceptable!")
    return string


def main():

    parser = argparse.ArgumentParser("Script to split FASTA sequences in a fixed number of multiple files.")
    parser.add_argument("-m", "--num-files", dest="num_files",
                        default=1000, type=positive,
                        help="Number of files to create. Default: %(default)d")
    parser.add_argument("fasta", type=pyfaidx.Fasta, help="Input FASTA file.")
    parser.add_argument("out", nargs="?", default=None,
                        help="Output prefix. Default: filename+split")
    args = parser.parse_args()

    if args.out is None:
        args.out = os.path.splitext(args.fasta)
    if os.path.dirname(args.out) and not os.path.exists(os.path.dirname(args.out)):
        os.makedirs(os.path.dirname(args.out))
    elif (os.path.dirname(args.out)
          and os.path.exists(os.path.dirname(args.out))
          and not os.path.isdir(os.path.dirname(args.out))):
        raise OSError("Invalid directory selected for output: {}".format(
            os.path.dirname(args.out)))

    max_sequences = max(int(len(args.fasta.keys()) / args.num_files), 1)
    zfiller = max(ceil(log(args.num_files, 10)), 3)

    outfile = None
    counter = 0

    for number, seq in enumerate(args.fasta):
        if number % max_sequences == 0:
            counter += 1
            if outfile:
                outfile.close()
            outfile = open("{}_{}.fasta".format(args.out, str(counter).zfill(zfiller)),
                           "wt")
        print(">{}".format(seq.long_name),
              *textwrap.wrap(str(seq), width=60),
              sep="\n",
              file=outfile)
    if outfile:
        outfile.close()

    while counter < args.num_files:
        counter += 1
        with open("{}_{}.fasta".format(args.out, str(counter).zfill(zfiller)),
                           "wt"):
            pass
    return

if __name__ == "__main__":
    main()
