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
    parser.add_argument("-d", "--distances", nargs="+", type=int, default=[1000, 2000, 5000])
    parser.add_argument("genome")
    parser.add_argument("gff3")
    parser.add_argument("gene_list")
    args = parser.parse_args()

    logger = Mikado.utilities.log_utils.create_default_logger("default")

    genome = pyfaidx.Fasta(args.genome)

    with open(args.gff3) as gff3:
        namespace = argparse.Namespace
        namespace.reference = gff3
        namespace.exclude_utr = False
        namespace.protein_coding = False
        # Use Mikado compare functions to load the index from the GFF3
        # "genes" is a dictionary of Gene objects, having as keys the gene names
        # "positions" is a dictionary of the form: [chrom][(start, end)] = [GID1, GID2, ...]
        genes, positions = Mikado.scales.compare.load_index(args, logger)
        # Create a dictionary of interval trees, one per chromosome
        indexer = collections.defaultdict(list).fromkeys(positions)
        for chrom in indexer:
            indexer[chrom] = IntervalTree.from_tuples(positions[chrom].keys())

    max_distance = max(args.distances)

    with open(args.gene_list) as gene_list:
        gids = [_.rstrip() for _ in gene_list]
        for gid in gids:
            if gid not in genes:
                exc = IndexError("{} not found in the index!".format(gid))
                logger.exception(exc)
                continue
            chrom, start, end, strand = (genes[gid].chrom,
                                         genes[gid].start,
                                         genes[gid].end,
                                         genes[gid].strand)
            if chrom not in genome:
                exc = IndexError("Chromosome {} not found in the genome!".format(chrom))
                logger.exception(exc)
                continue

            # If the gene is on the minus strand, the promoter is further down
            if strand == "-":
                key = (start, min(end + max_distance, len(genome[chrom])))
            else:
                # otherwise it is on the 5' side
                key = (max(0, start - max_distance), end)

            # Find all genes which are near
            neighbours = Mikado.scales.assigner.Assigner.find_neighbours(indexer.get(chrom, IntervalTree()),
                                                                         key, distance=0)

            # Find all the genes which are in the neighbourhood, remove the obvious case ..
            for neighbour in neighbours:

                print(positions[chrom][neighbour[0]])


    return

main.__doc__ = __doc__

if __name__ == "__main__":
    main()

