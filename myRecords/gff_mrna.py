#!/usr/bin/env python3

import sys, argparse,copy, operator,re
from myRecords import GFF
from numpy import mean,median
from scipy.stats.mstats import mquantiles

class mrna:
    
    def __init__(self, genome=None, stop_out=False):
        '''Optional parameter: genome (a FASTA SeqIO.dict object).'''

        self.mrna=None
        self.fasta=genome
        self.id =None
        self.chrom=""
        self.start=None
        self.end=None
        self.strand=""
        self.exons=[]
        self.cds=[]
        self.utrs=[]
        self.cds_introns=[]
        self.introns=[]
        self.monoexonic=False
        self.cds_monoexonic=False
        self.gene=None
        self.added_codons=False
        self.stop_out=stop_out


    def append(self, feature):
        '''Add a gffLine to this mrna record. The feature type is used to determine which action should be taken.
        If the feature is "mRNA" or "transcript", this is considered the "mother" feature and will be used to set the following attributes:
        - id
        - mrna (a copy of the original gffLine object)
        - chrom
        - start
        - end
        - strand
        - gene
        If the feature is "exon", it will be added to the list of exons (attribute self.exons)
        If the feature is "CDS", it will be added to the list of CDSs (attribute self.cds)'''

        if feature.feature=="CDS":
            self.cds.append(feature)
        elif re.search("UTR", feature.feature):
            self.utrs.append(feature)
        elif feature.feature=="exon":
            self.exons.append(feature)
        elif feature.feature=="transcript" or "rna" in feature.feature.lower():
            self.id=feature.id
            self.mrna=copy.deepcopy(feature)
            self.chrom=self.mrna.chrom
            self.start=self.mrna.start
            self.end=self.mrna.end
            self.strand=self.mrna.strand
            if "Parent" in feature.attributes:
                self.gene = feature.attributes['Parent']
            else:
                #Added for compatibility with non-gene records.
                self.gene=None

    def calc_introns(self):

        '''This method calculates the introns of the class, as gffLine objects.
        All the introns are stored in a specified list (self.introns)'''

        assert len(self.exons)>0, str(self.mrna)
        if len(self.exons)==1: return #No introns
        self.exons = sorted(self.exons, key=operator.attrgetter("start","end"))
        #Create introns
        for x in range(len(self.exons)-1):
            first=self.exons[x]
            second=self.exons[x+1]
            intron=copy.deepcopy(first)
            intron.feature="intron"
            intron.start=first.end+1
            intron.end=second.start-1
            intron.id = "{0}:intron_{1}".format(
                re.sub(".?exon", "", intron.id),
                x)
            self.introns.append(intron)
        
    def calc_cds_introns(self):
        '''This method calculates the introns of the class, as gffLine objects, but only for the CDS object.
        All the introns are stored in a specified list (self.cds_introns)'''

        if len(self.cds)==1: return #No introns
        self.cds = sorted(self.cds, key=operator.attrgetter("start","end"))
        #Create introns
        for x in range(len(self.cds)-1):
            first=self.cds[x]
            second=self.cds[x+1]
            intron=copy.deepcopy(first)
            intron.feature="intron"
            intron.start=first.end+1
            intron.end=second.start-1
            intron.id = "{0}:cds_intron_{1}".format(
                re.sub(".?cds", "", intron.id),
                x)
            self.cds_introns.append(intron)

    def calc_stats(self):

        '''This method calculates the lengths of introns, exons, mRNA, etc. It is used most predominantly in the gff_stats program.'''

        #Calculate introns
        self.calc_introns()
        self.calc_cds_introns()
        self.exon_lengths = [exon.end-exon.start for exon in self.exons]
        self.cds_lengths = [cds.end-cds.start for cds in self.cds]
        self.intron_lengths = [intron.end-intron.start for intron in self.introns]
        self.cds_intron_lengths = [intron.end - intron.start for intron in self.cds_introns]
        self.mrna_length = sum(self.exon_lengths)
        self.cds_length = sum(self.cds_lengths)
        if len(self.exons)==1: self.monoexonic=True
        if len(self.cds)==1: self.monocds=True

        #Delete raw data
        del self.exons
        del self.cds
        del self.introns
        del self.mrna

    def add_start_stop(self, genome=None):

        '''This method calculates the position of the start and stop codons of the mRNA and uses the genome file to verify
        that the positions indicated really contain a start or a stop codon.'''

        print("Starting to add codons", file=sys.stderr)

        if self.added_codons==True:
            print("Nothing to do", file=sys.stderr)
            return
        if self.fasta is None and genome is not None:
            
            self.fasta=genome

        if self.fasta is not None:
            self.sequence=self.fasta[self.chrom]
            print(self.sequence.id, file=sys.stderr)

        self.cds = sorted(self.cds, key=operator.attrgetter("start","end"))
        self.exons = sorted(self.exons, key=operator.attrgetter("start","end"))

        stop_exon = None
        start_exon = None

        if self.strand=="-":
            start_exon=copy.deepcopy(self.cds[-1])
            stop_exon=copy.deepcopy(self.cds[0])

            start_exon.feature="start_codon"
            start_exon.phase=0
            start_exon.start=start_exon.end-2
            start_sequence=str(self.sequence[start_exon.start-1:start_exon.end].seq.reverse_complement())
            if start_sequence!="ATG":
                start_exon=None
            
            stop_exon.feature="stop_codon"
            stop_exon.phase=0
            stop_exon.end=stop_exon.start+2
            stop_sequence=str(self.sequence[stop_exon.start-1:stop_exon.end].seq.reverse_complement())
            if stop_sequence not in ("TAG","TAA","TGA"):
                stop_exon=None
            elif self.stop_out:
                self.cds[0].start+=3

                        
        elif self.strand=="+":
            start_exon=copy.deepcopy(self.cds[0])
            stop_exon=copy.deepcopy(self.cds[-1])

            start_exon.feature="start_codon"
            start_exon.phase=0
            start_exon.end=start_exon.start+2
            try:
                start_sequence=str(self.sequence[start_exon.start-1:start_exon.end].seq)
            except KeyError:
                raise KeyError(str(start_exon),list(self.fasta.keys()))
            except ValueError:
                raise ValueError(str(start_exon), list(self.fasta.id))
                
            if start_sequence!="ATG":
                start_exon=None

            self.cds[-1].end-=3

            stop_exon.feature="stop_codon"
            stop_exon.phase=0
            stop_exon.start=stop_exon.end-2
            self.cds=[start_exon]+self.cds+[stop_exon]
            stop_sequence=str(self.sequence[stop_exon.start-1:stop_exon.end].seq)
            if stop_sequence not in ("TAG","TAA","TGA"):
                stop_exon=None
            elif self.stop_out:
                self.cds[-1].end-=3
            
        if stop_exon!=None:
            stop_exon.id = "{0}:stop_codon".format( re.sub(".?(cds|CDS:)", "", stop_exon.id) )
            self.cds.append(stop_exon)
        if start_exon is not None:
            start_exon.id = "{0}:start_codon".format( re.sub("(.?cds|CDS:)", "", start_exon.id) )
            self.cds.append(start_exon)
            

        self.added_codons=True
        self.cds=[x for x in self.cds if x!=None]
        

    def __str__(self):
        '''This method will print all the records inside the object, starting from the mRNA and continuing with CDS and exons.
        CDS, exons, start and stop codons will be interleaved according to their coordinates on the genome.'''
        if self.fasta is not None:
            self.add_start_stop()
        if "mrna" in self.__dict__:
            lines=[str(self.mrna)]
            exons=sorted([x for x in self.cds+self.exons+self.utrs if x is not None], key=operator.attrgetter("start","end"))
            for exon in exons:
                lines.append(str(exon))
        else:
             line="\t".join([str(s) for s in [self.id, self.start, self.end]])
             lines=[line]
             
        return "\n".join(lines)
