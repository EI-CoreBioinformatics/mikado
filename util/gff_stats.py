#!/usr/bin/env python3


import sys, argparse
from shanghai_lib.parsers import GFF,GTF
from shanghai_lib.loci_objects import transcript
from numpy import mean,median,array
from scipy.stats.mstats import mquantiles
import io,os,re
import scipy
import numpy
import shanghai_lib


numpy.seterr(all="ignore") #Suppress warnings
scipy.seterr(all="ignore") #Suppress warnings

class transcriptComputer(transcript.transcript):
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.exon_lengths = []
        self.cds_exon_lengthts = []
        self.utr_exon_lengths = []
        
        self.intron_lengths = []
        self.cds_intron_lengths = []
        self.utr_intron_lengths = []
    
    def finalize(self):
        super().finalize()
        self.exon_lengths = [e[1]-e[0]+1 for e in self.exons]
        self.cds_exon_lengths = [c[1]-c[0]+1 for c in self.selected_cds]
        self.utr_exon_lengths = [u[1]-u[0]+1 for u in self.three_utr+self.five_utr]

        self.intron_lengths = [i[1]-i[0]+1 for i in self.introns]
        self.cds_intron_lengths = [i[1]-i[0] for i in self.selected_cds_introns]
        self.utr_intron_lengths = [i[1]-i[0] for i in filter(lambda j: j not in self.selected_cds_introns, self.introns)]


class Calculator:
    
    def __init__(self, gff, is_gtf=False, only_coding=False):
        
        '''Constructor function'''
        
        self.gff = gff
        self.is_gtf = is_gtf
        self.only_coding = only_coding
        
        
    def __call__(self):
        
        self.genes = dict()
        self.transcript_to_genes = dict()

        currentGene = None
        for record in self.gff:
            if record.header is True:
                continue
            if record.is_parent is True:
                self.genes[record.id] = dict(("gene", record), ("transcripts", dict()))
                to_remove = set()
                for tid in self.genes[currentGene]["transcripts"]:
                    try:
                        self.genes[currentGene]["transcripts"][tid].finalize()
                    except shanghai_lib.exceptions.InvalidTranscript:
                        to_remove.add(tid)
                    if self.only_coding is True and self.genes[currentGene]["transcripts"][tid].selected_cds_length==0:
                        to_remove.add(tid)
                for tid in to_remove: del self.genes[currentGene]["transcripts"][tid]
                if len(self.genes[currentGene]["transcripts"])==0:
                    del self.genes[currentGene]
                    
                currentGene = record.id
            elif record.is_transcript is True:
                assert record.parent[0]==currentGene
                self.genes[record.parent[0]]["transcripts"][record.id] = transcriptComputer(record)
            else:
                for parent in self.parent:
                    self.genes[currentGene]["transcripts"][parent].addExon(record)
        
        
                
            


        
def main():

    global mrna
    def to_gff(string):
        return GFF.GFF3(open(string))

    parser=argparse.ArgumentParser()
    parser.add_argument('--print_exon_num', action="store_true", default=False, help="Flag. If set, print for each mRNA the number of exons and CDSs.")
    parser.add_argument('gff', type=to_gff, help="GFF file to parse.")
    parser.add_argument('out', type=argparse.FileType('w'), default=sys.stdout, nargs='?')
    args=parser.parse_args()

    gff_mrnas, genes, intergenic = dict(), dict(), dict()

    if args.print_exon_num:
        for record in gff_mrnas:
            print(gff_mrnas[record].id, len(gff_mrnas[record].exon_lengths), len(gff_mrnas[record].cds_lengths), sep="\t", file=args.out)
        return

    print("Number of genes", len(genes), sep=":\t", file=args.out)
    print("Number of transcripts", len(gff_mrnas), sep=":\t", file=args.out)
    mrna_count=[genes[gene]['count'] for gene in genes]
    print("Mean number of mRNAs per gene", mean(mrna_count), sep=":\t", file=args.out)
    print("Mean intergenic distance", mean(intergenic), sep=":\t", file=args.out)
    print("Median intergenic distance", median(intergenic), sep=":\t", file=args.out)

    print("Mean mRNA length", mean([gff_mrnas[record].mrna_length for record in gff_mrnas]), sep=":\t", file=args.out)
    print("Mean CDS length", mean([gff_mrnas[record].cds_length for record in gff_mrnas]), sep=":\t", file=args.out)
    exon_lengths = []
    cds_lengths = []
    intron_lengths = []
    cds_intron_lengths = []
    exon_num=[]
    cds_num=[]
    initial_cds=[]
    final_cds=[]

    for parent in gff_mrnas:
        record=gff_mrnas[parent]
        exon_lengths+=record.exon_lengths
        cds_lengths+=record.cds_lengths
        intron_lengths+=record.intron_lengths
        cds_intron_lengths+=record.cds_intron_lengths
        exon_num.append(len(record.exon_lengths))
        if len(record.cds_lengths)>0:
            cds_num.append(len(record.cds_lengths))
            if record.strand=="-":
                initial_cds.append(record.cds_lengths[-1])
                final_cds.append(record.cds_lengths[0])
            else:
                initial_cds.append(record.cds_lengths[0])
                final_cds.append(record.cds_lengths[-1])
                
    print("Mean Exon length", mean(exon_lengths), sep=":\t", file=args.out)
    print("Median Exon length", median(exon_lengths), sep=":\t", file=args.out)
    print("95th percentile of Exon length", *mquantiles(exon_lengths, prob=[0.9]), sep=":\t", file=args.out)
    print("Mean CDS Exon length", mean(cds_lengths), sep=":\t", file=args.out)
    print("Median CDS length", median(cds_lengths), sep=":\t", file=args.out)
    print("95th percentile of CDS length", *mquantiles(cds_lengths, prob=[0.9]), sep=":\t", file=args.out)

    #Initial and Final CDS exon lengths
    
    print("Mean CDS initial length", mean(initial_cds), sep=":\t", file=args.out)
    print("Median CDS initial length", median(initial_cds), sep=":\t", file=args.out)
    print("95th percentile of initial CDS length", *mquantiles(initial_cds, prob=[0.95]), sep=":\t", file=args.out)
    
    print("Mean CDS final length", mean(final_cds), sep=":\t", file=args.out)
    print("Median CDS final length", median(final_cds), sep=":\t", file=args.out)
    print("95th percentile of final CDS length", *mquantiles(final_cds, prob=[0.95]), sep=":\t", file=args.out)

    ##
    if intron_lengths!=[]:
        mean_intron=mean(intron_lengths)
        max_intron=max(intron_lengths)
        q90_intron=mquantiles(intron_lengths, prob=[0.9])
        q95_intron=mquantiles(intron_lengths, prob=[0.95])
        mean_cds_intron=mean(cds_intron_lengths)
    else:
        mean_intron="NA"
        max_intron="NA"
        q90_intron=["NA"]
        q95_intron=["NA"]
        mean_cds_intron="NA"
        
    print("Mean intron length", mean_intron, sep=":\t", file=args.out)
    print("Max intron length", max_intron, sep=":\t", file=args.out)
    print("90th percentile of intron length", *q90_intron, sep=":\t", file=args.out)
    print("95th percentile of intron length", *q95_intron, sep=":\t", file=args.out)
    print("Mean CDS intron length", mean_cds_intron, sep=":\t", file=args.out)
    print("Mean number of exons", mean(exon_num), sep=":\t", file=args.out)
    print("Median number of exons", median(exon_num), sep=":\t", file=args.out)
    print("Mean number of CDS exons", mean(cds_num), sep=":\t", file=args.out)
    print("Median number of CDS exons", median(cds_num), sep=":\t", file=args.out)
    monoexonic = [record for record in gff_mrnas if len(gff_mrnas[record].exon_lengths)==1]
    monocds = [record for record in gff_mrnas if len(gff_mrnas[record].cds_lengths)==1]
    print("Monoexonic", len(monoexonic), sep=":\t", file=args.out)
    print("% monoexonic", "{0}%".format(
            round(100*len(monoexonic)/len(gff_mrnas),4)), sep=":\t", file=args.out)
    monoexonic_lengths=[]
    for parent in monoexonic:
        record=gff_mrnas[parent]
        monoexonic_lengths+=record.exon_lengths
    monocds_lengths=[]
    for parent in monocds:
        record=gff_mrnas[parent]
        monocds_lengths+=record.cds_lengths

    if monoexonic_lengths!=[]:
        mean_mono=mean(monoexonic_lengths)
    else:
        mean_mono="--"
    print("Monoexonic mean length", mean_mono, sep=":\t", file=args.out)
    print("10th percentile of monoexonic length", *mquantiles(monoexonic_lengths, prob=[0.1]), sep=":\t", file=args.out)
    print("95th percentile of monoexonic length", *mquantiles(monoexonic_lengths, prob=[0.95]), sep=":\t", file=args.out)

    print("MonoCDS", len(monocds), sep=":\t", file=args.out)
    print("% monoCDS", "{0}%".format(round(100*len(monocds)/len(gff_mrnas),4)), sep=":\t", file=args.out)
    if len(monocds_lengths)>0:
        mean_monocds=mean(monocds_lengths)
    else:
        mean_monocds="--"
    print("MonoCDS mean length", mean_monocds, sep=":\t", file=args.out)
    print("10th percentile of monoCDS length", *mquantiles(monocds_lengths, prob=[0.1]), sep=":\t", file=args.out)
    print("95th percentile of monoCDS length", *mquantiles(monocds_lengths, prob=[0.95]), sep=":\t", file=args.out)

if __name__=='__main__': main()
