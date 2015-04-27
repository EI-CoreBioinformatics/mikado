#!/usr/bin/env python3

from loci_objects.superlocus import superlocus
from loci_objects.transcript import transcript
import sys,argparse
from loci_objects.GFF import GFF3
from loci_objects.GTF import GTF
from loci_objects.bed12 import bed12
from loci_objects.abstractlocus import abstractlocus

def locus_printer( slocus, args, cds_dict=None ):
#     if slocus is None:
#         return
    slocus.load_cds(cds_dict)
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
    parser.add_argument("--cds", type=argparse.FileType("r"), default=None)
    parser.add_argument("gff", type=argparse.FileType("r"))
    
    args=parser.parse_args()

    currentLocus=None
    currentTranscript=None
    currentChrom=None
    
    if args.gff.name[-3:]=="gtf":
        print("Parsing a GTF", file=sys.stderr)
        rower=GTF(args.gff)
    else: rower=GFF3(args.gff)
    
    cds_dict=None
    
    if args.cds is not None:
        cds_dict = dict()
        
        for line in args.cds:
            line=bed12(line)
            if line.header is True: continue
            if line.chrom not in cds_dict:
                cds_dict[line.chrom]=[]
            to_append = True
            indices_to_remove = []
            for index in range(len(cds_dict[line.chrom])):
                entry=cds_dict[line.chrom][index]
                overl=abstractlocus.overlap( (entry.cdsStart,entry.cdsEnd), (line.cdsStart,line.cdsEnd) )
                if overl==entry.cds_len:
                    indices_to_remove.append(index)
                elif overl==line.cds_len:
                    to_append=False
                    break
            if to_append is True:
                for index in indices_to_remove:
                    del cds_dict[line.chrom][index]
                cds_dict[line.chrom].append(line)
    
    for row in rower:
        
        if row.header is True: continue
        if row.chrom!=currentChrom:
            if currentChrom is not None:
                if currentTranscript is None:
                    pass
                elif superlocus.in_locus(currentLocus, currentTranscript):
                    currentLocus.add_transcript_to_locus(currentTranscript)
                else:
                    locus_printer(currentLocus, args, cds_dict=cds_dict)
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
                    locus_printer(currentLocus, args, cds_dict=cds_dict)
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
            locus_printer(currentLocus, args, cds_dict=cds_dict)
            currentLocus=superlocus(currentTranscript, stranded=False)
        locus_printer(currentLocus, args, cds_dict=cds_dict)
        
if __name__=="__main__":
    main()