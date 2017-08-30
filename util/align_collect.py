#!/usr/bin/env python3

import os
import argparse


class Stats:

    file = ""
    seqs = 0
    mapped = 0
    mapped_n_paired = 0
    unmapped = 0
    properly_paired = 0
    paired = 0
    duplicated = 0
    MQ0 = 0

    def __init__(self):
        self.data = []

    def __str__(self):

        return "\t".join([str(_) for _ in [self.file,
                                           self.seqs,
                                           self.mapped,
                                           self.mapped_n_paired,
                                           self.unmapped,
                                           self.properly_paired,
                                           self.paired,
                                           self.duplicated,
                                           self.MQ0
                                           ]])

    @staticmethod
    def header():
        return "\t".join(["File",
                          "sequences",
                          "reads_mapped",
                          "reads_mapped_and_paired",
                          "reads_unmapped",
                          "reads_properly_paired",
                          "reads_paired",
                          "reads_duplicated",
                          "reads_MQ0"])


def main():
    parser = argparse.ArgumentParser("Script to collect info from multiple samtools stats files")
    parser.add_argument("input", nargs='+', help="The list of samtools stats file to process")
    args = parser.parse_args()

    print(Stats.header())

    for stats_file in args.input:
        with open(stats_file) as f:
            s = Stats()
            s.file = os.path.split(stats_file)[1]
            for line in f:
                if line.startswith("SN\t"):
                    parts = line.split("\t")
                    if parts[1] == "sequences:":
                        s.seqs = int(parts[2])
                    elif parts[1] == "reads mapped:":
                        s.mapped = int(parts[2])
                    elif parts[1] == "reads mapped and paired:":
                        s.mapped_n_paired = int(parts[2])
                    elif parts[1] == "reads unmapped:":
                        s.unmapped = int(parts[2])
                    elif parts[1] == "reads properly paired:":
                        s.properly_paired = int(parts[2])
                    elif parts[1] == "reads paired:":
                        s.paired = int(parts[2])
                    elif parts[1] == "reads duplicated:":
                        s.duplicated = int(parts[2])
                    elif parts[1] == "reads MQ0:":
                        s.MQ0 = int(parts[2])
            print(s)


main()
