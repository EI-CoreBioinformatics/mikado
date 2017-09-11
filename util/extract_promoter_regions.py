#!/usr/bin/env python3

import Mikado
import pyfaidx
import argparse
import sys
import collections
from Mikado.utilities.intervaltree import IntervalTree


__doc__ = """Little script to """


def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("-o", "--out", type=argparse.FileType("wt"), default=sys.stdout)
    parser.add_argument("genome")
    parser.add_argument("gff3")
    parser.add_argument("gene_list")
    args = parser.parse_args()

    with open(args.gff3) as gff3:
        namespace = argparse.Namespace
        namespace.reference = gff3
        namespace.exclude_utr = False
        namespace.protein_coding = False
        null_logger = Mikado.utilities.log_utils.create_null_logger()
        genes, positions = Mikado.scales.compare.load_index(args, null_logger)
        indexer = collections.defaultdict(list).fromkeys(positions)
        for chrom in indexer:
            indexer[chrom] = IntervalTree.from_tuples(positions[chrom].keys())


    with open(args.gene_list) as gene_list:
        gids = [_.rstrip() for _ in gene_list]
        for gid in gids:
            


    return

main.__doc__ = __doc__

if __name__ == "__main__":
    main()

