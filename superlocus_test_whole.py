#!/usr/bin/env python3

from superlocus import superlocus
from transcript import transcript
import sys,argparse
from myRecords.GFF import GFF3
from myRecords.GTF import GTF

def locus_printer( slocus, args ):
    stranded_loci = sorted(list(slocus.split_strands()))
    for stranded_locus in stranded_loci:
        stranded_locus.define_subloci()
        print(stranded_locus, file=args.sub_out)
        stranded_locus.define_monosubloci()
        print(stranded_locus, file=args.mono_out)
        stranded_locus.define_loci()
        print(stranded_locus, file=args.locus_out)
    

def main():
    
    parser=argparse.ArgumentParser("Quick test utility for the superlocus classes.")
    parser.add_argument("--sub_out", type=argparse.FileType("w"), required=True)
    parser.add_argument("--mono_out", type=argparse.FileType("w"), required=True)
    parser.add_argument("--locus_out", type=argparse.FileType("w"), required=True)
    parser.add_argument("gff", type=argparse.FileType("r"))
    
    args=parser.parse_args()

    currentLocus=None
    currentTranscript=None
    currentChrom=None
    
    if args.gff.name[-3:]=="gtf":
        print("Parsing a GTF", file=sys.stderr)
        rower=GTF(args.gff)
    else: rower=GFF3(args.gff)
    
    for row in rower:
        
        if row.header is True: continue
        if row.chrom!=currentChrom:
            if currentChrom is not None:
                if currentTranscript is None:
                    pass
                elif superlocus.in_locus(currentLocus, currentTranscript):
                    currentLocus.add_transcript_to_locus(currentTranscript)
                else:
                    locus_printer(currentLocus, args)
                    currentLocus=superlocus(currentTranscript, stranded=False)
                locus_printer(currentLocus, args)
            #print("Changed chrom", file=sys.stderr)
            currentChrom=row.chrom
            currentTranscript=None
            currentLocus=None
            
        if row.feature=="transcript" or "RNA" in row.feature.upper():
            if currentLocus is not None:
                if currentTranscript is None:
                    pass
                elif superlocus.in_locus(currentLocus, currentTranscript):
                    currentLocus.add_transcript_to_locus(currentTranscript)
                else:
                    locus_printer(currentLocus, args)
                    currentLocus=superlocus(currentTranscript, stranded=False)
            elif currentLocus is None:
                if currentTranscript is not None:
                    currentLocus=superlocus(currentTranscript, stranded=False)
            currentTranscript=transcript(row)
#            print("New transcript: {0}".format(currentTranscript.id))
        elif row.feature in ("exon", "CDS") or "UTR" in row.feature.upper():
            currentTranscript.addExon(row)
        else:
            continue
        
    if currentLocus is not None:
        if currentTranscript is None:
            pass
        elif superlocus.in_locus(currentLocus, currentTranscript):
            currentLocus.add_transcript_to_locus(currentTranscript)
        else:
            locus_printer(currentLocus, args)
            currentLocus=superlocus(currentTranscript, stranded=False)
        locus_printer(currentLocus, args)
        
if __name__=="__main__":
    main()