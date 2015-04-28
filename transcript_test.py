#!/usr/bin/env python3

import sys,argparse
from collections import namedtuple
from loci_objects.abstractlocus import abstractlocus
from loci_objects.GFF import GFF3
from loci_objects.transcript import transcript
from loci_objects.bed12 import BED12
from Bio import SeqIO

def main():
    
    parser=argparse.ArgumentParser("Simple script to verify the loading of the CDS BEDs")
    parser.add_argument("--bed", type=argparse.FileType("r"), required=True)
    parser.add_argument("--transcript_fasta", default=None)
    parser.add_argument("gff", type=argparse.FileType("r"), nargs="?")
    parser.add_argument("out", nargs="?", default=sys.stdout, type=argparse.FileType("w"))
    args=parser.parse_args()
    
#     bedTup = namedtuple("bed", ["chrom","start","end", "name", "score", "strand",
#                                 "cdsStart", "cdsEnd", "rgb", "blockCount", "blockSizes", "blockStarts"]
#                         )
#     
    if args.transcript_fasta is not None:
        args.transcript_fasta = SeqIO.index(args.transcript_fasta, "fasta")
    
    bed_entries = dict()
    for line in BED12(args.bed, fasta_index=args.transcript_fasta):
        if line.header is True: continue
        if line.chrom not in bed_entries:
            bed_entries[line.chrom] = [line]
        else:
            to_append = True
            for entry in bed_entries[line.chrom]:
                if entry.strand!=line.strand:
                    continue
                overl = abstractlocus.overlap(
                                              (entry.cdsStart,entry.cdsEnd),
                                              (line.cdsStart, line.cdsEnd)
                                              )
                if overl==entry.cdsEnd-entry.cdsStart: # earlier overlap is completely contained
                    bed_entries[line.chrom].remove(entry)
                elif overl==line.cdsEnd-line.cdsStart:
                    to_append = False
                    break
            if to_append is True:
                bed_entries[line.chrom].append(line)
                
    
    currentTranscript = None
    bed_finals=dict.fromkeys(bed_entries.keys())
    for tid in bed_entries:
        for b in bed_entries[tid]:
            print(str(b), b.has_start, b.has_stop, sep="\t")
            
#     for line in GFF3(args.gff):
#         if line.header is True:
#             continue
#         elif line.feature=="transcript" or "RNA" in line.feature.upper():
#             if currentTranscript is not None:
#                 currentTranscript.load_cds(bed_finals)
# #                 print(currentTranscript, file=args.out)
# #                 print("###", file=args.out)
#                 print(currentTranscript.id, currentTranscript.max_internal_cds_length,
#                       currentTranscript.cdna_length, currentTranscript.cds_length,
#                       currentTranscript.internal_cds_num, currentTranscript.internal_cds_lengths )
#             
#             currentTranscript=transcript(line)
#             
#         elif line.feature=="exon":
#             currentTranscript.addExon(line)
#             
#     if currentTranscript is not None:
#         currentTranscript.load_cds(bed_finals)
#         print(currentTranscript.id, currentTranscript.max_internal_cds_length,
#               currentTranscript.cdna_length, currentTranscript.cds_length,
#               currentTranscript.internal_cds_num, currentTranscript.internal_cds_lengths )
        
if __name__=="__main__":
    main()