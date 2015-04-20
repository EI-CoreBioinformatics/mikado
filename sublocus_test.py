#!/usr/bin/env python3

from sublocus import sublocus
from transcript import transcript
import sys,argparse
from myRecords.GFF import GFF3
#from myRecords.GTF import GTF

def main():
    
    parser=argparse.ArgumentParser("Quick test utility for the superlocus classes.")
    parser.add_argument("-s","--scores", required=True, type=argparse.FileType("r")  )
    parser.add_argument("gff", type=argparse.FileType("r"))
    parser.add_argument("out", type=argparse.FileType("w"), default=sys.stdout, nargs="?")
    args=parser.parse_args()

    currentTranscript=None
    currentSub=None
    rower=GFF3(args.gff)
    
    scores=dict()
    for line in args.scores:
        tid,score=line.rstrip().split()
        scores[tid]=float(score)
    
    for row in rower:
        
        if row.header is True: continue
        
        elif row.feature=="sublocus":
            if currentSub is not None:
                currentSub.add_transcript_to_locus(currentTranscript)
                currentSub.load_scores(scores)
                currentSub.define_loci()
                print(currentSub, file=args.out)
                currentTranscript=None
            
            currentSub=sublocus(row)
        elif row.feature=="transcript":
            try:
                currentSub.add_transcript_to_locus(currentTranscript)
            except Exception as err:
                print(currentSub)
                raise err
            currentTranscript=transcript(row)
            
        elif row.feature in ("exon","CDS"):
            currentTranscript.addExon(row)
        
        else:
            continue
        
    if currentTranscript is not None:
        if currentSub is not None:
            currentSub.add_transcript_to_locus(currentTranscript)
            currentSub.load_scores(scores)
            currentSub.define_loci()
            print(currentSub, file=args.out)
            


if __name__=="__main__":
    main()