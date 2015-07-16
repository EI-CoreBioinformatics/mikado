#!/usr/bin/env python3

import sys

import argparse
from mikado_lib.parsers.GTF import GTF
from mikado_lib.loci_objects.transcript import transcript

def main():

    parser=argparse.ArgumentParser("Quick utility to extract features from a GTF with certain coordinates.")
    parser.add_argument("--chrom", required=True)
    parser.add_argument("-as", "--assume-sorted", dest="assume_sorted",
                        default=False, action="store_true")
    parser.add_argument("--start", type=int, default=float("-inf"))
    parser.add_argument("--end", type=int, default=float("inf"))
    parser.add_argument("gtf", type=argparse.FileType("r"))
    parser.add_argument("out", type=argparse.FileType("w"), nargs="?", default=sys.stdout)
    args=parser.parse_args()

    if args.start>=args.end:
        raise ValueError("Start greater than end: {0}\t{1}".format(args.start,args.end))

    currentTranscript=None
    with GTF(args.gtf) as gtf:
        for row in gtf:
            if row.chrom!=args.chrom:
                continue
            else:
                if row.is_transcript is True:
                    if currentTranscript is not None and currentTranscript.start>=args.start and currentTranscript.end<=args.end:
                        print(currentTranscript.__str__(to_gtf=True), file=args.out)
                    if args.assume_sorted is True and row.end>args.end:
                        break
                    currentTranscript=transcript(row)
                else:
                    currentTranscript.addExon(row)

    if currentTranscript is not None and currentTranscript.start>=args.start and currentTranscript.end<=args.end:
        print(currentTranscript.__str__(to_gtf=True), file=args.out)

if __name__=="__main__":
    main()
            
        
                        
                        
                
                
                    
                    
    
