#!/usr/bin/env python3

from loci_objects.superlocus import superlocus
from loci_objects.transcript import transcript
import sys,argparse
from loci_objects.GFF import GFF3
from loci_objects.GTF import GTF

def add_transcripts_old(transcript_instance, sloci):
    if transcript_instance is None:
        return sloci
    transcript_instance.finalize()
    foundLocus=False
    if len(sloci)==0:
        slocus=superlocus(transcript_instance)
        sloci.append(slocus)
    else:
        for sindex in reversed(range(len(sloci))):
            slocus=sloci[sindex]
            if superlocus.in_locus(slocus, transcript_instance) is True:
                sloci[sindex].add_transcript_to_locus(transcript_instance)
                foundLocus=True
                break
        if foundLocus is False:
            slocus=superlocus(transcript_instance)
            sloci.append(slocus)
    return sloci 

def locus_printer( slocus, out=sys.stdout ):
    stranded_loci = sorted(list(slocus.split_strands()))
    for stranded_locus in stranded_loci:
        print(stranded_locus, file=out)
    

def main():
    
    parser=argparse.ArgumentParser("Quick test utility for the superlocus classes.")
    parser.add_argument("gff", type=argparse.FileType("r"))
    parser.add_argument("out", type=argparse.FileType("w"), default=sys.stdout, nargs="?")
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
                    locus_printer(currentLocus, out=args.out)
                    currentLocus=superlocus(currentTranscript, stranded=False)
                locus_printer(currentLocus, out=args.out)
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
                    locus_printer(currentLocus, out=args.out)
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
            locus_printer(currentLocus, out=args.out)
            currentLocus=superlocus(currentTranscript, stranded=False)
        locus_printer(currentLocus, out=args.out)
        
if __name__=="__main__":
    main()