#!/usr/bin/env python3

import sys,argparse,os
import json

from loci_objects.superlocus import superlocus
from loci_objects.transcript import transcript
from loci_objects.GFF import GFF3
from loci_objects.GTF import GTF
from loci_objects.bed12 import BED12
from loci_objects.abstractlocus import abstractlocus

def locus_printer( slocus, args, cds_dict=None ):
#     if slocus is None:
#         return
    slocus.load_cds(cds_dict)
    stranded_loci = sorted(list(slocus.split_strands()))
    
    if not os.path.exists(args.sub_metrics ):
        sub_metrics = open(args.sub_metrics,"w")
    else:
        sub_metrics = open(args.sub_metrics,"a")
        
    if not os.path.exists(args.locus_metrics ):
        locus_metrics = open(args.locus_metrics,"w")
    else:
        locus_metrics = open(args.locus_metrics,"a")
    
    for stranded_locus in stranded_loci:
        stranded_locus.define_subloci()
        stranded_locus.print_subloci_metrics(sub_metrics)
        print(stranded_locus, file=args.sub_out)
        stranded_locus.define_monosubloci()
        print(stranded_locus, file=args.mono_out)
        stranded_locus.calculate_mono_metrics()
        stranded_locus.print_monoholder_metrics(locus_metrics)
        stranded_locus.define_loci()
        print(stranded_locus, file=args.locus_out)
    
    locus_metrics.close()
    sub_metrics.close()

def main():
    
    def to_json(string):
        with open(string) as json_file:
            json_dict = json.load(json_file)
        return json_dict
    
    def to_index(string):
        try:
            from Bio import SeqIO
        except ImportError as err:
            print("Error importing the Bio module, no indexing performed:\n{0}",format(err) )
            return None
        return SeqIO.index(string, "fasta")
    
    parser=argparse.ArgumentParser("Quick test utility for the superlocus classes.")
    parser.add_argument("--json_conf", type=argparse.FileType("r"), required=True, help="JSON configuration file for scoring transcripts.")
    parser.add_argument("--sub_out", type=argparse.FileType("w"), required=True)
    parser.add_argument("--mono_out", type=argparse.FileType("w"), required=True)
    parser.add_argument("--locus_out", type=argparse.FileType("w"), required=True)
    parser.add_argument("--cds", type=argparse.FileType("r"), default=None)
    parser.add_argument("--transcript_fasta", type=to_index, default=None)
    parser.add_argument("gff", type=argparse.FileType("r"))
    
    args=parser.parse_args()

    args.json_conf = to_json(args.json_conf.name)

    currentLocus=None
    currentTranscript=None
    currentChrom=None

    args.sub_metrics = "{0}.metrics".format( args.sub_out.name) 
    args.locus_metrics = "{0}.metrics".format( args.locus_out.name)
    if os.path.exists(args.sub_metrics):
        os.remove(args.sub_metrics) 
    if os.path.exists(args.locus_metrics):
        os.remove(args.locus_metrics)
    
    if args.gff.name[-3:]=="gtf":
        rower=GTF(args.gff)
    else: rower=GFF3(args.gff)
    
    cds_dict=None
    
    if args.cds is not None:
        cds_dict = dict()
        
        for line in BED12(args.cds, fasta_index=args.transcript_fasta):
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
                    currentLocus=superlocus(currentTranscript, stranded=False, json_dict = args.json_conf)
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
                    currentLocus=superlocus(currentTranscript, stranded=False, json_dict = args.json_conf)
            elif currentLocus is None:
                if currentTranscript is not None:
                    currentLocus=superlocus(currentTranscript, stranded=False, json_dict = args.json_conf)
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
            currentLocus=superlocus(currentTranscript, stranded=False, json_dict = args.json_conf)
        locus_printer(currentLocus, args, cds_dict=cds_dict)
        
if __name__=="__main__":
    main()