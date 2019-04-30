#!/usr/bin/env python3

import argparse
import sys
import pysam
from Mikado.transcripts.transcript import Transcript
from collections import Counter


def to_bam(string):
    return pysam.AlignmentFile(string, mode="rb")


def main():
    parser = argparse.ArgumentParser("Script to convert from BAM to GTF, for PB alignments")
    parser.add_argument("--strict", action="store_true", default=False,
                        help="Switch. If set, this script will never output multiexonic transcripts \
                        without a defined strand.")
    parser.add_argument("--outfmt", choices=["gtf", "bed12"], default="gtf")
    parser.add_argument("bam", type=to_bam, help="Input BAM file")
    parser.add_argument("out", nargs="?", default=sys.stdout, type=argparse.FileType("wt"),
                        help="Optional output file")
    args = parser.parse_args()

    # M 0 alignment match (can be a sequence match or mismatch)
    # I 1 insertion to the reference
    # D 2 deletion from the reference
    # N 3 skipped region from the reference
    # S 4 soft clipping (clipped sequences present in SEQ)
    # H 5 hard clipping (clipped sequences NOT present in SEQ)
    # P 6 padding (silent deletion from padded reference)
    # = 7 sequence match
    # X 8 sequence mismatch

    name_counter = Counter()

    for record in args.bam:
        if record.is_unmapped is True:
            continue
        transcript = Transcript(record, accept_undefined_multi=(not args.strict))
        if name_counter.get(record.query_name):
            name = "{}_{}".format(record.query_name, name_counter.get(record.query_name))
        else:
            name = record.query_name

        if name != transcript.id:
            transcript.alias = transcript.id
            transcript.id = name

        transcript.parent = transcript.attributes["gene_id"] = "{0}.gene".format(name)
        name_counter.update([record.query_name])
        transcript.source = "bam2gtf"
        print(transcript.format(args.outfmt), file=args.out)


main()
