#!/usr/bin/env python3

import Mikado
import argparse
# import sys
import collections
from Mikado.utilities.intervaltree import IntervalTree
# import os
import subprocess
import multiprocessing
import tempfile
from scipy.stats.mstats import mquantiles
from numpy import arange


class Geomancer(multiprocessing.Process):

    def __init__(self, gid_pipe: multiprocessing.Pipe,
                 positions: dict,
                 out_file: str,
                 index: dict):

        super(Geomancer, self).__init__()
        self.gid_pipe = gid_pipe
        self.index = index
        self.out = open(out_file, "wt")
        self.positions = positions

    def run(self):

        while True:

            # This part basically replicates the Mikado.scales.assigner.Assigner.find_neighbours function

            (gid, chrom, start, end, strand) = self.gid_pipe.get()
            if gid is None:
                self.gid_pipe.put((None, None, None, None, None))
                self.out.flush()
                self.out.close()
                return

            found = self.index[chrom].search(start, end, max_distance=10**8, num_intervals=1, strict=False)
            if not found:
                line = [gid, "{}:{}..{}".format(chrom, start, end)] + ["NA"] * 6
                self.out.write("\t".join(str(_) for _ in line) + "\n")
                continue

            before, after = [], []

            for key in found:
                key = (key.start, key.end)
                if key[1] < start:
                    before.append(key)
                elif key[0] > end:
                    after.append(key)
                else:
                    found_gids = [_ for _ in self.positions[chrom][key] if _ != gid]
                    if found_gids:
                        before.append(key)
                        after.append(key)
                    else:
                        continue

            before = sorted(before)
            after = sorted(after)

            if len(before) == 0:
                gid_before = "NA"
                position_before = "NA"
                distance_before = "NA"
            else:
                gid_before = ",".join(_ for _ in self.positions[chrom][before[-1]])
                position_before = "{}:{}..{}".format(chrom, *before[-1])
                distance_before = max(0, start - before[-1][1] + 1)
            if len(after) == 0:
                gid_after = "NA"
                distance_after = "NA"
                position_after = "NA"
            else:
                gid_after = ",".join(_ for _ in self.positions[chrom][after[0]])
                distance_after = max(0, after[0][0] - end + 1)
                position_after = "{}:{}..{}".format(chrom, *after[0])

            line = [gid, "{}:{}..{}".format(chrom, start, end),
                    gid_before, distance_before, position_before,
                    gid_after, distance_after, position_after]

            self.out.write("\t".join(str(_) for _ in line) + "\n")
            continue

    def join(self, timeout=None):
        # self.__close_handles()
        # self.terminate()
        super(Geomancer, self).join(timeout=timeout)


def main():

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("-o", "--out", type=str, default="distances.txt")
    parser.add_argument("-l", "--log", default=None)
    parser.add_argument("-lv", "--log-level", default="INFO", choices=["DEBUG", "INFO", "WARN", "ERROR", "CRITICAL"],
                        dest="log_level")
    parser.add_argument("-p", "--procs", type=int, default=1)
    parser.add_argument("gff3")

    args = parser.parse_args()

    logger = Mikado.utilities.log_utils.create_logger_from_conf({"log_settings": {"log": args.log,
                                                                                  "log_level": args.log_level}},
                                                                "geomancer",
                                                                mode="w")
    logger.info("Starting to load the GFF3 index")
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
    logger.info("Loaded the index")

    # Create different processes
    gid_pipe = multiprocessing.JoinableQueue()
    # res_pipe =

    procs = []
    logger.info("Creating children processes")
    for num in range(args.procs):
        temp = tempfile.NamedTemporaryFile(mode="w+t", delete=False)
        temp.close()
        _ = Geomancer(gid_pipe, positions, temp.name, indexer)
        _.start()
        procs.append(_)
    logger.info("Created children processes")

    logger.info("Starting to analyse {} genes".format(len(genes)))
    # Lookup in a small set is MUCH faster than lookup in a small np.array
    steps = set(mquantiles(range(1, len(genes.keys())), prob=arange(0.1, 1.1, 0.1)).astype(int))

    for counter, gid in enumerate(genes.keys(), 1):
        if counter in steps:
            logger.info("Analysed {} genes ({}%)".format(counter, int(counter * 100/len(genes))))
        chrom, start, end, strand = genes[gid].chrom, genes[gid].start, genes[gid].end, genes[gid].strand
        gid_pipe.put((gid, chrom, start, end, strand))

    gid_pipe.put((None, None, None, None, None))
    # [_.join() for _ in procs]

    with open(args.out, "wt") as out:
        out.close()
        for proc in procs:
            subprocess.call("cat {} >> {}".format(proc.out.name, args.out), shell=True)
            # subprocess.call("rm {}".format(proc.out.name), shell=True)

main()
