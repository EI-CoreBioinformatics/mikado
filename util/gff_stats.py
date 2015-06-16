#!/usr/bin/env python3


import sys, argparse, os.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from shanghai_lib.parsers import GFF,GTF
from shanghai_lib.loci_objects import transcript
from scipy.stats.mstats import mquantiles
# import io,os,re
import scipy
import numpy
import shanghai_lib
from collections import namedtuple, Counter
from array import array as cArray

numpy.seterr(all="ignore") #Suppress warnings
numpy.warnings.filterwarnings("ignore")
scipy.seterr(all="ignore") #Suppress warnings

class gene_object:
    
    def __init__(self, record, only_coding = False):
        
        self.transcripts = dict()
        
        for field in ['chrom', 'strand', 'start', 'end', 'id']:
            setattr(self, field, getattr(record, field))
        self.only_coding = only_coding
        
    def add_transcript(self, tcomputer):
        self.transcripts[tcomputer.id]=tcomputer
        self.start=min(self.start, tcomputer.start)
        self.end=max(self.end, tcomputer.end)

    def finalize(self):
        if len(self.transcripts)==0: return
        
        to_remove = set()
        for tid in self.transcripts:
            try:
                self.transcripts[tid]=self.transcripts[tid].as_tuple()
                if self.only_coding is True and self.transcripts[tid].selected_cds_length == 0:
                    to_remove.add(tid)            
            except shanghai_lib.exceptions.InvalidTranscript:
                to_remove.add(tid)
        if len(to_remove)==len(self.transcripts):
            self.transcripts=dict()
        else:
            for x in to_remove: del self.transcripts[x]


    def __lt__(self, other):
        if self.chrom!=other.chrom:
            return self.chrom<other.chrom
        else:
            if self.start!=other.start:
                return self.start<other.start
            else:
                return self.end<other.end
            
    def __eq__(self, other):
        if type(self)!=type(other): return False
        if self.chrom==other.chrom and self.strand == other.strand and self.start == other.start and self.end==other.end:
            return True
        return False
    
    @property
    def introns(self):
        return set(self.transcripts[tid].introns for tid in self.transcripts)
    
    @property
    def exons(self):
        return set(self.transcripts[tid].exons for tid in self.transcripts)
    
    @property
    def has_monoexonic(self):
        return any(len(self.transcripts[tid].introns)==0 for tid in self.transcripts.keys())
    
    @property
    def num_transcripts(self):
        return len(self.transcripts)

class transcriptComputer(transcript.transcript):
    
    data_fields = ["parent", 'chrom',
                   'start', 'end',
                   'introns', 'exons',
                   'exon_lengths', 'intron_lengths',
                   'cdna_length', 'selected_cds_length',
                   'cds_intron_lengths', 'cds_exon_lengths'   ]
    data_tuple = namedtuple("transcript_data", data_fields, verbose=False)
    
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
        self.utr_exon_lengths = [u[2]-u[1]+1 for u in self.three_utr+self.five_utr]

        self.intron_lengths = [i[1]-i[0]+1 for i in self.introns]
        self.cds_intron_lengths = [i[1]-i[0] for i in self.selected_cds_introns]
        self.utr_intron_lengths = [i[1]-i[0] for i in filter(lambda j: j not in self.selected_cds_introns, self.introns)]

    def as_tuple(self):
        '''Method to build a namedtuple containing only the basic information for stat building.
        
        We want to analyze the following:
        - cDNA length
        - CDS length
        - Exons (number and length)
        - CDS Exons (number and length)
        - Introns (number and length)
        - CDS Introns (number and length)
        '''
        
        self.finalize()
        constructor = dict()
        for field in self.data_fields:
            constructor[field]=getattr(self, field)
                    
        return self.data_tuple(**constructor)
        

class Calculator:
    
    def __init__(self, parsed_args):
        
        '''Constructor function'''
        
        self.gff = parsed_args.gff
        if type(self.gff) is GFF.GFF3:
            self.is_gff = True
        else:
            self.is_gff = False
        self.only_coding = parsed_args.only_coding
        self.out = parsed_args.out
        self.genes = dict()
        
    def parse_input(self):
                
        currentGene = None
        for record in self.gff:
            if record.header is True:
                continue
#             if record.feature == "superlocus":
#                 continue

            if record.feature == "locus" or (record.is_parent is True and record.is_transcript is False):
                if currentGene is not None:
                    self.genes[currentGene].finalize()
                    if self.genes[currentGene].num_transcripts == 0:
                        del self.genes[currentGene]
                self.genes[record.id] = gene_object(record, only_coding=self.only_coding)    
                currentGene = record.id
            elif record.is_transcript is True:
                if self.is_gff is False:
                    new_record = record.copy()
                    new_record.feature = "gene"
                    if new_record.gene not in self.genes:
                        if currentGene is not None:
                            self.genes[currentGene].finalize()
                            if self.genes[currentGene].num_transcripts == 0:
                                del self.genes[currentGene]
                        
                        self.genes[new_record.gene] = gene_object(new_record, only_coding=self.only_coding)
                        currentGene = new_record.gene
                    else:
                        assert record.parent[0]==currentGene, (currentGene,record.parent)
                else:
                    assert record.parent[0]==currentGene, (currentGene,record.parent)
                    
                self.genes[record.parent[0]].transcripts[record.id] = transcriptComputer(record)
            else:
                for parent in record.parent:
                    if parent not in self.genes[currentGene].transcripts: continue #Hack for the AT gff .. stupid "protein" features ...
                    self.genes[currentGene].transcripts[parent].addExon(record)

        self.genes[currentGene].finalize()

        
    def __call__(self):

        #self.transcript_to_genes = dict()


        self.parse_input()
        ordered_genes = sorted(self.genes.values())
        distances = []
        for index,gene in enumerate(ordered_genes[:-1]):
            next_gene = ordered_genes[index+1]
            if next_gene.chrom!=gene.chrom: continue
            distances.append(next_gene.start+1-gene.end)
            
        self.distances = numpy.array(distances)

        self.writer()


    def get_stats(self, row, array):
        
        row["Average"] = "{0:,.2f}".format( round(numpy.mean(array),2)) #Decimal to second digit precision
        
        counter_object = Counter(array)
        moder = [x for x in counter_object if counter_object[x]==counter_object.most_common(1)[0][1]]
        row["Mode"] = ";".join( str(x) for x in moder )
        quantiles = mquantiles(array, prob=[0,0.05, 0.10, 0.25, 0.5, 0.75,0.9, 0.95, 1]  )
        for key,val in zip(['Min', '5%', '10%', '25%', 'Median', '75%', '90%', '95%', 'Max'   ], quantiles):
            try:
                row[key] = "{0:,.0f}".format(val) #No decimal
            except:
                row[key] = val
        return row

    def writer(self):
        '''Method which creates the final output'''

        import csv
        fieldnames = ['Stat', 'Total', 'Average', 'Mode', 'Min', '5%', '10%', '25%', 'Median', '75%', '90%', '95%', 'Max'   ]
        rower = csv.DictWriter(self.out,
                               fieldnames,
                               delimiter="\t"  )

        rower.writeheader()
        row=dict()
        for key in fieldnames:
            row[key] = "NA"
        row["Stat"] = 'Number of genes'
        row['Total'] = len(self.genes)
        rower.writerow(row)
        row["Stat"] = 'Number of transcripts'
        row['Total'] = sum(self.genes[x].num_transcripts for x in self.genes  )   
        rower.writerow(row)
        
        row["Stat"] = 'Transcript per gene'
        t_per_g = numpy.array(list( self.genes[x].num_transcripts for x in self.genes  ))
        row = self.get_stats(row, t_per_g)
        rower.writerow(row)

        exons = cArray('i') #Done
        exons_coding = cArray('i')
        exon_num = cArray('i') #Done
        exon_num_coding = cArray('i')
        introns = cArray('i') #Done
        introns_coding = cArray('i')
        cds_introns = cArray('i') #Done
        cds_exons = cArray('i') # Done
        cds_exon_num = cArray('i') # Done
        cds_exon_num_coding = cArray('i')
        cdna_lengths = cArray('i') # Done
        cdna_lengths_coding = cArray('i')
        cds_lengths = cArray('i') # Done
        monoexonic_lengths = cArray('i')
        multiexonic_lengths = cArray('i')
        monocds_lengths = cArray('i')
        
        for gene in self.genes:
            for tid in self.genes[gene].transcripts:
                exons.extend(  self.genes[gene].transcripts[tid].exon_lengths  )
                exon_number =  len(self.genes[gene].transcripts[tid].exon_lengths)
                exon_num.append(exon_number)
                if exon_number==1:
                    monoexonic_lengths.append(self.genes[gene].transcripts[tid].cdna_length)
                else:
                    multiexonic_lengths.append(self.genes[gene].transcripts[tid].cdna_length)
                introns.extend( self.genes[gene].transcripts[tid].intron_lengths  )
                cds_introns.extend( self.genes[gene].transcripts[tid].cds_intron_lengths  )
                cds_exons.extend( self.genes[gene].transcripts[tid].cds_exon_lengths  )
                cds_num=len(self.genes[gene].transcripts[tid].cds_exon_lengths)
                if cds_num==1 or (exon_num==1 and self.genes[gene].transcripts[tid].selected_cds_length>0):
                    monocds_lengths.append(self.genes[gene].transcripts[tid].selected_cds_length  )
                cds_exon_num.append(cds_num)
                cdna_lengths.append( self.genes[gene].transcripts[tid].cdna_length  )
                cds_lengths.append( self.genes[gene].transcripts[tid].selected_cds_length  )
                if self.only_coding is False and self.genes[gene].transcripts[tid].selected_cds_length>0:
                    cdna_lengths_coding.append(self.genes[gene].transcripts[tid].cdna_length)
                    exons_coding.extend(self.genes[gene].transcripts[tid].exon_lengths)
                    exon_num_coding.append(len(self.genes[gene].transcripts[tid].exon_lengths))
                    cds_exon_num_coding.append(len(self.genes[gene].transcripts[tid].cds_exon_lengths))
                    introns_coding.extend(self.genes[gene].transcripts[tid].intron_lengths)

        row["Stat"] = 'CDNA lengths'
        row["Total"] = 'NA'
        ar = cdna_lengths
        row = self.get_stats(row, ar)
        rower.writerow(row)

        row["Stat"] = 'CDS lengths'
        row["Total"] = 'NA'
        ar = cds_lengths
        row = self.get_stats(row, ar)
        rower.writerow(row)

        if self.only_coding is False:
            row["Stat"] = "CDS lengths (mRNAs)"
            row["Total"] = 'NA'
            ar =  cArray( 'i', list(filter(lambda x: x>0, cds_lengths) ) )
            row = self.get_stats(row, ar)
            rower.writerow(row)

        row["Stat"] = 'Monoexonic transcripts'
        row['Total'] = len(monoexonic_lengths)
        ar = monoexonic_lengths
        row = self.get_stats(row, ar)
        rower.writerow(row)

        row["Stat"] = 'MonoCDS transcripts'
        row['Total'] = len(monocds_lengths)
        ar = monoexonic_lengths
        row = self.get_stats(row, ar)
        rower.writerow(row)
        
        row["Stat"] = 'Exons per transcript'
        ar = exon_num
        row = self.get_stats(row, ar)
        row["Total"]=len(exons)
        rower.writerow(row)
                
        if self.only_coding is False:
            row["Stat"] = 'Exons per transcript (mRNAs)'
            ar = exon_num_coding
            row = self.get_stats(row, ar)
            row["Total"]=len(exons)
            rower.writerow(row)
                
        row["Stat"] = 'Exon lengths'
        ar=exons
        row=self.get_stats(row, ar)
        row["Total"]="NA"
        rower.writerow(row)

        if self.only_coding is False:
            row["Stat"] = 'Exon lengths (mRNAs)'
            ar=exons_coding
            row=self.get_stats(row, ar)
            row["Total"]="NA"
            rower.writerow(row)
        
        row["Stat"] = "Intron lengths"
        ar = introns
        row = self.get_stats(row, ar)
        row["Total"]="NA"
        rower.writerow(row)

        if self.only_coding is False:
            row["Stat"] = "Intron lengths (mRNAs)"
            ar = introns_coding
            row = self.get_stats(row, ar)
            row["Total"]="NA"
            rower.writerow(row)

        row["Stat"] = "CDS exons per transcript"
        ar = cds_exon_num
        row = self.get_stats(row, ar)
        row["Total"] = len(cds_exons)
        rower.writerow(row)
        
        if self.only_coding is False:
            row["Stat"] = "CDS exons per transcript (mRNAs)"
            ar = cds_exon_num_coding
            row = self.get_stats(row, ar)
            row["Total"] = len(cds_exons)
            rower.writerow(row)
        
        row["Stat"] = "CDS exon lengths"
        ar = cds_exons
        row = self.get_stats(row, ar)
        row["Total"] = "NA"
        rower.writerow(row)        
        
        row["Stat"] = "CDS Intron lengths"
        ar = cds_introns
        row = self.get_stats(row, ar)
        row["Total"]="NA"
        rower.writerow(row)

        

        
def main():

    def to_gff(string):
        if string.endswith("gtf"):
            return GTF.GTF(open(string))
        elif string.endswith('gff') or string.endswith('gff3'):
            return GFF.GFF3(open(string))
        else:
            raise ValueError('Unrecognized format')

    parser=argparse.ArgumentParser()
    parser.add_argument('--only-coding', dest="only_coding", action="store_true", default=False )
    parser.add_argument('gff', type=to_gff, help="GFF file to parse.")
    parser.add_argument('out', type=argparse.FileType('w'), default=sys.stdout, nargs='?')
    args=parser.parse_args()

    calculator = Calculator(args)
    calculator()


if __name__=='__main__': main()
