#!/usr/bin/env python3

from sublocus import sublocus
from transcript import transcript
import sys,argparse
import csv
from myRecords.GFF import GFF3
#from myRecords.GTF import GTF

def main():
    
    parser=argparse.ArgumentParser("Quick test utility for the superlocus classes.")
    parser.add_argument("-s","--scores", default=None, required=False, type=argparse.FileType("r")  )
    parser.add_argument("-m", "--metrics_file", type=argparse.FileType('w'),
                        help="File to write the transcript metrics out to.")
    parser.add_argument("gff", type=argparse.FileType("r"))
    parser.add_argument("out", type=argparse.FileType("w"), default=sys.stdout, nargs="?")
    args=parser.parse_args()

    currentTranscript=None
    currentSub=None
    rower=GFF3(args.gff)
    
    
    metrics = ["id","parent"] 
    
    
    if args.scores is not None:
        scores=dict()
        for line in args.scores:
            tid,score=line.rstrip().split()
            scores[tid]=float(score)
    
    for row in rower:
        
        if row.header is True: continue
        
        elif row.feature=="sublocus":
            if currentSub is not None:
                currentSub.add_transcript_to_locus(currentTranscript)
                #print(currentSub, file=args.out)
                if args.scores is None:
                    currentSub.calculate_scores()
                else:
                    currentSub.load_scores(scores)
                if args.metrics_file is not None:
                    if type(args.metrics_file)!=csv.DictWriter: # Initialize
                        metrics.extend(sorted(currentSub.available_metrics))
                        args.metrics_file = csv.DictWriter(args.metrics_file, metrics, delimiter="\t"  )
                        args.metrics_file.writeheader()
                    currentSub.print_metrics(args.metrics_file)
                currentSub.define_monosubloci()
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
            
        elif row.feature in ("exon","CDS") or "utr" in row.feature.lower():
            currentTranscript.addExon(row)
        
        else:
            continue
        
    if currentTranscript is not None:
        if currentSub is not None:
            currentSub.add_transcript_to_locus(currentTranscript)
            if args.scores is not None:
                currentSub.load_scores(scores)
            else:
                currentSub.calculate_scores()
            if args.metrics_file is not None:
                if type(args.metrics_file)!=csv.DictWriter: # Initialize
                    metrics.extend(sorted(currentSub.available_metrics))
                    args.metrics_file = csv.DictWriter(args.metrics_file, metrics, delimiter="\t"  )
                    args.metrics_file.writeheader()
                currentSub.print_metrics(args.metrics_file)
                
            currentSub.define_monosubloci()
            print(currentSub, file=args.out)
            
if __name__=="__main__":
    main()