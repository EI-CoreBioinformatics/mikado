#!/usr/bin/env python3

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from loci_objects.GTF import GTF
from loci_objects.transcript import transcript
import argparse

def main():
    
    parser=argparse.ArgumentParser('Utility to sort GTF files with a transcript attribute')
    parser.add_argument("gtf", type=argparse.FileType("r"),
                    help="Input GTF")
    parser.add_argument("out", default=sys.stdout, nargs="?",
                        type=argparse.FileType("w"),
                        help="Output file. Default: stdout.")
    args=parser.parse_args()

    transcripts=[]

    currentTranscript=None

    for record in GTF(args.gtf.name):
        if record.is_transcript is True:
            if currentTranscript is not None:
                assert currentTranscript.exon_num>0, (currentTranscript.id, record.id)
                transcripts.append(currentTranscript)
            assert record.id is not None, str(record)
            currentTranscript=transcript(record)
        else:
            currentTranscript.addExon(record)

    transcripts.append(currentTranscript)
    for tr in sorted(transcripts):
        print(tr.__str__(to_gtf=True), file=args.out)
        
if __name__=='__main__': main()