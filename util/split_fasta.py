#!/usr/bin/env python3

from itertools import chain, repeat
import argparse
import os
import textwrap
import pyfaidx
from math import log, ceil, floor
from itertools import zip_longest

def positive(string):

    string = int(string)
    if string < 1:
        raise ValueError("Only positive values are acceptable!")
    return string


def fixed_grouper(number, iterable, padvalue=None):
    # return [iterable[x:x + n] for x in range(0, len(MyList), n)]
    return zip_longest(*[chain(iterable, repeat(padvalue, number - 1))] * number)


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

    zfiller = max(ceil(log(args.num_files, 10)), 3)

    number = 0

    for number, group in enumerate(fixed_grouper(
            ceil(len(args.fasta.keys())/args.num_files),
            args.fasta.keys())):
        with open("{}_{}.fasta".format(args.out, str(number + 1).zfill(zfiller)),
                           "wt") as outfile:
            for sequence in group:
                if sequence is not None:
                    print(">{}".format(args.fasta[sequence].long_name),
                          *textwrap.wrap(str(args.fasta[sequence]), width=60),
                          sep="\n", file=outfile)

    while number + 1 < args.num_files:
        number += 1
        with open("{}_{}.fasta".format(args.out, str(number + 1).zfill(zfiller)),
                  "wt") as outfile:
            pass

    return

if __name__ == "__main__":
    main()
