#!/usr/bin/env python3
# coding: utf-8

"""GFF +> BED12 converter for junction files."""

import sys
import argparse
import mikado_lib.parsers
import mikado_lib.loci_objects


def main():
    """Main script function"""

    parser = argparse.ArgumentParser("GFF=>BED12 converter")
    parser.add_argument("gff", type=argparse.FileType())
    parser.add_argument("out", nargs="?", type=argparse.FileType("w"), default=sys.stdout)
    args = parser.parse_args()

    transcript = None

    for row in mikado_lib.parsers.GFF.GFF3(args.gff):
        if row.is_parent and not row.is_transcript:
            continue
        if row.is_exon:
            transcript.add_exon(row)
        elif row.is_transcript:
            if transcript is not None:
                transcript.finalize()
                assert len(transcript.introns) == 1
                introns = sorted(transcript.introns)
                bed12 = mikado_lib.parsers.bed12.BED12()
                bed12.chrom = transcript.chrom
                bed12.strand = transcript.strand
                bed12.score = transcript.score
                bed12.start = transcript.start
                bed12.end = transcript.end
                bed12.thickStart = introns[0][0]
                bed12.thickEnd = introns[0][1]
                bed12.name = transcript.id
                bed12.rgb = "255,0,0"
                bed12.blockCount = 2
                bed12.blockSizes = [(e[1] - e[0]) for e in transcript.exons]
                bed12.blockStarts = []
                bed12.blockStarts = [0, bed12.blockSizes[0] + introns[0][1] - introns[0][0]]
                print(bed12, file=args.out)

            transcript = mikado_lib.loci_objects.transcript.Transcript(row)

    if transcript is not None:
        transcript.finalize()
        assert len(transcript.introns) == 1
        introns = sorted(transcript.introns)
        bed12 = mikado_lib.parsers.bed12.BED12()
        bed12.chrom = transcript.chrom
        bed12.strand = transcript.strand
        bed12.start = transcript.start
        bed12.score = transcript.score
        bed12.end = transcript.end
        bed12.thickStart = introns[0][0]
        bed12.thickEnd = introns[0][1]
        bed12.name = transcript.id
        bed12.rgb = "255,0,0"
        bed12.blockCount = 2
        bed12.blockSizes = [(e[1] - e[0]) for e in transcript.exons]
        bed12.blockStarts = []
        bed12.blockStarts = [0, bed12.blockSizes[0] + introns[0][1] - introns[0][0]]
        print(bed12, file=args.out)


if __name__ == "__main__":
    main()
