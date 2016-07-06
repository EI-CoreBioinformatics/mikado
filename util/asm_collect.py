#!/usr/bin/env python3

import os
import argparse


class AsmStats:

    file = ""
    genes = 0
    monoexonic_genes = 0
    transcripts = 0
    transcripts_per_gene = 0.0
    transcript_len_mean = 0.0
    monoexonic_transcripts = 0
    exons = 0
    exons_per_transcript = 0.0
    exon_len_mean = 0.0

    def __init__(self):
        self.data = []

    def __str__(self):

        return "\t".join([str(_) for _ in
                          [self.file, self.genes, self.monoexonic_genes, self.transcripts,
                           self.transcripts_per_gene, self.transcript_len_mean,
                           self.monoexonic_transcripts, self.exons, self.exons_per_transcript,
                           self.exon_len_mean]])

    @staticmethod
    def header():

        return "\t".join(["File",
                          "genes",
                          "monoexonic_genes",
                          "transcripts",
                          "transcripts_per_gene",
                          "transcript_len_mean",
                          "monoexonic_transcripts",
                          "exons",
                          "exons_per_transcript",
                          "exon_len_mean"
                          ])


def main():
    parser = argparse.ArgumentParser("Script to collect info from multiple mikado util stats files")
    parser.add_argument("input", nargs='+', help="The list of mikado util stats file to process")
    args = parser.parse_args()

    print(AsmStats.header())

    for stats_file in args.input:
        with open(stats_file) as f:
            s = AsmStats()
            s.file = os.path.split(stats_file)[1]
            f.readline()
            # f.readline()
            # f.readline()
            # f.readline()
            for line in f:
                parts = line.split("\t")
                col2 = float(parts[1].replace(',','').replace("NA","0"))
                col3 = float(parts[2].replace(',','').replace("NA","0.0"))
                if parts[0] == "Number of genes":
                    s.genes = col2
                elif parts[0] == "Number of monoexonic genes":
                    s.monoexonic_genes = col2
                elif parts[0] == "Transcripts per gene":
                    s.transcripts = col2
                    s.transcripts_per_gene = col3
                elif parts[0] == "CDNA lengths":
                    s.transcript_len_mean = col3
                elif parts[0] == "Monoexonic transcripts":
                    s.monoexonic_transcripts = col2
                elif parts[0] == "Exons per transcript":
                    s.exons = col2
                    s.exons_per_transcript = col3
                elif parts[0] == "Exon lengths":
                    s.exon_len_mean = col3
            print(s)

main()
