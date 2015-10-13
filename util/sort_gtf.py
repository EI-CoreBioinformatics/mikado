#!/usr/bin/env python3
# coding: utf-8

"""Utility to sort GTF files with a transcript attribute"""

import sys
from Mikado.parsers.GTF import GTF
from Mikado.loci_objects.transcript import Transcript
import argparse


def main():
    """
    Main script function
    """

    parser = argparse.ArgumentParser('Utility to sort GTF files with a transcript attribute')
    parser.add_argument("gtf", type=argparse.FileType(),
                        help="Input GTF")
    parser.add_argument("out", default=sys.stdout, nargs="?",
                        type=argparse.FileType("w"),
                        help="Output file. Default: stdout.")
    args = parser.parse_args()

    transcripts = []

    transcript = None

    for record in GTF(args.gtf.name):
        if record.is_transcript is True:
            if transcript is not None:
                assert transcript.exon_num > 0, (transcript.id, record.id)
                transcripts.append(transcript)
            assert record.id is not None, str(record)
            transcript = Transcript(record)
        else:
            transcript.add_exon(record)

    transcripts.append(transcript)
    for tr in sorted(transcripts):
        print(tr.__str__(to_gtf=True), file=args.out)


if __name__ == '__main__':
    main()
