#!/usr/bin/env python3

import sys,argparse
from collections import namedtuple
from loci_objects.abstractlocus import abstractlocus
from loci_objects.GFF import GFF3
from loci_objects.transcript import transcript

def main():
    
    parser=argparse.ArgumentParser("Simple script to verify the loading of the CDS BEDs")
    parser.add_argument("--bed", type=argparse.FileType("r"), required=True)
    parser.add_argument("gff", type=argparse.FileType("r"))
    parser.add_argument("out", nargs="?", default=sys.stdout, type=argparse.FileType("w"))
    args=parser.parse_args()
    
    bedTup = namedtuple("bed", ["chrom","start","end", "name", "score", "strand",
                                "cdsStart", "cdsEnd", "rgb", "blockCount", "blockSizes", "blockStarts"]
                        )
    
    bed_entries = dict()
    for line in args.bed:
        if line[0]=="#": continue # Skip comments
        bl = line.rstrip().split("\t")[:12]
        for index in range(len(bl)):
            try:
                bl[index]=int(bl[index])
                if index in (1,6,11):
                    bl[index]=bl[index]+1
            except:
                pass
        bed_line = bedTup(*bl) # convert the line on the fly
        assert bed_line.strand in ("+","-"), (bed_line)
        if bed_line.chrom not in bed_entries:
            bed_entries[bed_line.chrom] = [bed_line]
        else:
            to_append = True
            for entry in bed_entries[bed_line.chrom]:
                if entry.strand!=bed_line.strand:
                    continue
                overl = abstractlocus.overlap(
                                              (entry.cdsStart,entry.cdsEnd),
                                              (bed_line.cdsStart, bed_line.cdsEnd)
                                              )
                if overl==entry.cdsEnd-entry.cdsStart: # earlier overlap is completely contained
                    bed_entries[bed_line.chrom].remove(entry)
                elif overl==bed_line.cdsEnd-bed_line.cdsStart:
                    to_append = False
                    break
            if to_append is True:
                bed_entries[bed_line.chrom].append(bed_line)
                
    
    currentTranscript = None
    bed_finals=dict.fromkeys(bed_entries.keys())
    for tid in bed_finals:
        bed_finals[tid]=[]
        for t in bed_entries[tid]:
            bed_finals[tid].append(( t.cdsStart, t.cdsEnd, t.strand, ) )
        
            
    for line in GFF3(args.gff):
        if line.header is True:
            continue
        elif line.feature=="transcript" or "RNA" in line.feature.upper():
            if currentTranscript is not None:
                currentTranscript.load_cds(bed_finals)
#                 print(currentTranscript, file=args.out)
#                 print("###", file=args.out)
                print(currentTranscript.id, currentTranscript.max_internal_cds_length,
                      currentTranscript.cdna_length, currentTranscript.cds_length,
                      currentTranscript.internal_cds_num, currentTranscript.internal_cds_lengths )
            
            currentTranscript=transcript(line)
            
        elif line.feature=="exon":
            currentTranscript.addExon(line)
            
    if currentTranscript is not None:
        currentTranscript.load_cds(bed_finals)
        print(currentTranscript.id, currentTranscript.max_internal_cds_length,
              currentTranscript.cdna_length, currentTranscript.cds_length,
              currentTranscript.internal_cds_num, currentTranscript.internal_cds_lengths )
        
if __name__=="__main__":
    main()