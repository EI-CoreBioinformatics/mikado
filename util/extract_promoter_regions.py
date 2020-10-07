#!/usr/bin/env python3

import Mikado
import pyfaidx
import argparse
import sys
import collections
from Mikado.utilities.intervaltree import IntervalTree
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gzip


__doc__ = """Little script to extract promoter regions from genes."""


def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("-o", "--out", type=str, default="promoters")
    parser.add_argument("-l", "--log", default=None)
    parser.add_argument("-lv", "--log-level", default="WARN", choices=["DEBUG", "INFO", "WARN", "ERROR", "CRITICAL"],
                        dest="log_level")
    parser.add_argument("-d", "--distances", nargs="+", type=int, default=[1000, 2000, 5000])
    parser.add_argument("-nn", "--no-neighbours", dest="no_neighbours", action="store_true", default=False,
                        help="Ignore the presence of neighbours when extracting genes.")
    parser.add_argument("-eu", "--exclude-utr", dest="exclude_utr", default=False, action="store_true")
    parser.add_argument("-z", "--gzip", default=False, action="store_true",
                        help="Output will be compressed in GZip format.")
    parser.add_argument("genome")
    parser.add_argument("gff3")
    parser.add_argument("gene_list")
    args = parser.parse_args()

    logger = Mikado.utilities.log_utils.create_logger_from_conf({"log_settings": {"log": args.log,
                                                                                  "log_level": args.log_level}},
                                                                "extractor",
                                                                mode="w")

    max_distance = max(args.distances)
    out_files = dict()
    args.distances = sorted([_ for _ in args.distances if _ > 0])
    if not args.distances:
        exc = ValueError("I need at least one positive integer distance!")
        logger.exception(exc)
        sys.exit(1)
    for distance in args.distances:
        if args.gzip is True:
            out_files[distance] = gzip.open("{}-{}bp.fasta.gz".format(os.path.splitext(args.out)[0],
                                                                      distance), "wt")
        else:
            out_files[distance] = open("{}-{}bp.fasta".format(os.path.splitext(args.out)[0],
                                                     distance), "wt")

    logger.info("Starting to load the genome")
    genome = pyfaidx.Fasta(args.genome)
    logger.info("Loaded the genome")

    logger.info("Starting to load the GFF3 index")
    with open(args.gff3) as gff3:
        namespace = argparse.Namespace
        namespace.reference = gff3
        namespace.exclude_utr = args.exclude_utr
        namespace.protein_coding = False
        # Use Mikado compare functions to load the index from the GFF3
        # "genes" is a dictionary of Gene objects, having as keys the gene names
        # "positions" is a dictionary of the form: [chrom][(start, end)] = [GID1, GID2, ...]
        genes, positions = Mikado.scales.compare.load_index(args, logger)
        # Create a dictionary of interval trees, one per chromosome
        indexer = collections.defaultdict(list).fromkeys(positions)
        for chrom in indexer:
            indexer[chrom] = IntervalTree.from_tuples(positions[chrom].keys())
    logger.info("Loaded the index")

    with open(args.gene_list) as gene_list:
        gids = [_.rstrip() for _ in gene_list]
        logger.info("Starting to extract sequences for {} genes".format(len(gids)))
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
            if args.no_neighbours is False:

                neighbours = Mikado.scales.assignment.assigner.Assigner.find_neighbours(indexer.get(chrom, IntervalTree()),
                                                                                        key, distance=0)
                # This is a list of the form [((start, end), distance), ...] where "(start, end)" is a key for the
                # "positions" dictionary, above

                # Find all the genes which are in the neighbourhood, remove the obvious case of the identity ..
                def is_before(gid_coords, key, strand):
                    if strand == "-":
                        return (Mikado.utilities.overlap(gid_coords, key) >= 0) or gid_coords[1] < key[0]
                    else:
                        return (Mikado.utilities.overlap(gid_coords, key) >= 0) or gid_coords[0] > key[1]

                neighbours = [_[0] for _ in neighbours if
                              is_before((start, end), _[0], strand) and gid not in positions[chrom][_[0]]]
            else:
                neighbours = []

            if not neighbours:
                # No neighbours found, we can grab everything
                for distance in args.distances:
                    try:
                        if strand == "-":
                            chunk = (max(0, end), min(end + distance, len(genome[chrom])))
                            seq = genome[chrom][chunk[0]:chunk[1]].reverse.complement.seq
                        else:
                            chunk = (max(0, start - 1 - distance), start - 1)
                            seq = genome[chrom][chunk[0]:chunk[1]].seq

                        seq = SeqRecord(Seq(seq), id="{}-prom-{}".format(gid, distance),
                                        description="{}{}:{}-{}".format(chrom, strand, chunk[0], chunk[1]))
                        print(seq.format("fasta"), file=out_files[distance], end='')
                    except ValueError as err:
                        logger.error("Error extracting the promoter for %s, distance %d. Error:\n%s",
                                     gid, distance, err)
                        continue
            else:
                # We have some neighbours, we have to select the maximum distance we can go to
                logger.warning("{} neighbours found for {}: {}".format(len(neighbours), gid, neighbours))
                if any([Mikado.utilities.overlap((start, end), _) >= 0 for _ in neighbours]):
                    logger.warning("Overlapping genes found for {}. Skipping".format(gid))
                    continue
                for distance in args.distances:
                    try:
                        if strand == "-":
                            max_point = min([_[0] for _ in neighbours])
                            if end + distance > max_point:
                                continue
                            chunk = (max(0, end), min(max_point, min(end + distance, len(genome[chrom]))))
                            seq = genome[chrom][chunk[0]:chunk[1]].reverse.complement.seq
                            description = "{}{}:{}-{}".format(chrom, strand, chunk[1], chunk[0])
                        else:
                            min_point = max([_[1] for _ in neighbours])
                            if start - distance < min_point:
                                continue
                            chunk = (max(0, start - 1 - distance), start - 1)
                            seq = genome[chrom][chunk[0]:chunk[1]].seq
                            description = "{}{}:{}-{}".format(chrom, strand, chunk[0], chunk[1])
                        seq = SeqRecord(Seq(seq), id="{}-prom-{}".format(gid, distance),
                                        description=description)
                        print(seq.format("fasta"), file=out_files[distance], end='')

                    except ValueError as err:
                        logger.error("Error extracting the promoter for %s, distance %d. Error:\n%s",
                                     gid, distance, err)
                        continue


    logger.info("Finished")
    return

main.__doc__ = __doc__

if __name__ == "__main__":
    main()

