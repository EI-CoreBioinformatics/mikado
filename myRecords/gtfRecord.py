#!/usr/bin/env python3

import sys,re
from myRecords import GTF
from Bio import SeqIO, Seq, SeqRecord
import logging


default_logger=logging.getLogger()
stream_handler=logging.StreamHandler()
default_logger.addHandler(stream_handler)

##In the maker-derived GTF the stop codon is EXTERNAL to any other annotation, be it CDS or UTR.
##Start codons are instead considered part of the CDS.


class gtfRecord(object):
    def __init__(self, *args, logger=default_logger):
        self.name=None
        if len(args)>0: self.name=args[0]
        self.start=None
        self.strict=True
        self.stop=None
        self.exons=[]
        self.strand=None
        self.chrom=None
        self.introns=[]
        self.cds=None
        self.donors, self.acceptors = {}, {}
        self.length=0
        self.region=None
        self.false_starts=[]
        self.false_stops=[]
        self.false_donors=[]
        self.false_acceptors=[]
        self.failed=False #Flag. Set to True if something goes amiss.
        self.non_canonical=False #Flag. Set to True if some splicing junction is non-canonical
        self.logger=logger

    def add(self, gtfLine):
        '''This method adds the input line to the GTF model.'''
        if not self.name: self.name=gtfLine.transcript
        if not self.strand: self.strand=gtfLine.strand
        if not self.chrom: self.chrom=gtfLine.chrom

        if gtfLine:
            if not self.name==gtfLine.transcript:
                raise AttributeError("This record has TID {0} and you tried to add a record for TID {1}".format(self.name, gtfLine.transcript))
            if gtfLine.feature=="start_codon":
                self.start=gtfLine
            elif gtfLine.feature=="stop_codon":
                self.stop=gtfLine
            elif gtfLine.feature=="CDS": #This is a brutish hack. I will have to add exon features to the GTF
                self.exons.append(gtfLine)

    def calculateStart(self):
        length=self.length
        self.logger.debug("Calculating start signal FASTA.")
        if not self.start:
            self.logger.error("Start information absent!")
            self.failed=True
            return

        '''This function will print the Start codon padded by the length selected when calling the class.'''

        if self.strand=="+":
            seq=self.fasta[self.start.start-length-1:self.start.start+length-1]
        else:
            seq=self.fasta[self.start.start-length+2:self.start.start+length+2]
            seq=seq.reverse_complement()
        seq.id=self.name
        seq.description="Start"
        self.start_fasta=seq

    def calculateStop(self):
        '''This function will print the Start codon padded by the length selected when calling the class.'''
        length=self.length
        if not self.stop:
            self.logger.error("Stop information absent!")
            self.failed=True

        if self.strand=="+":
            seq=self.fasta[self.stop.start-(length-1):self.stop.start+length+1]
        else:
            seq=self.fasta[self.stop.start-length:self.stop.start+length]
            seq=seq.reverse_complement()
        seq.id=self.name
        seq.description="Stop"
        self.stop_fasta=seq

    def calculateCDS(self):
        '''This method derives and returns the CDS sequence for the model.'''
        self.logger.debug("Calculating the CDS features.")
        self.exons=sorted(set(self.exons), key=lambda x: x.start)
        cds=""
        if self.strand=="+":
            for exon in self.exons: cds+=str(self.fasta[exon.start-1:exon.stop].seq)
        else:
            for exon in reversed(self.exons):
                cds_seq=str(self.fasta[exon.start-1:exon.stop].reverse_complement().seq)
                # assert isinstance(cds_seq, str)
                cds+=cds_seq
            
        self.cds=SeqRecord.SeqRecord(Seq.Seq(cds), id=self.name, description=self.strand) 


    def calculateSplice(self):

        '''This function calculates three things:
              - intron: a **list** of fasta records
              - acceptors: a **dictionary** of fasta records, indexed on the position.
              - donors: a **dictionary** of fasta records, indexed on the position.'''

        length=self.length
        self.calculateCDS()
        self.logger.debug("Starting to calculate the splice junctions.")
        if len(self.exons)==1:
            self.logger.info("Gene is monoexonic: no splice junctions")
            return 

        introns, donors, acceptors = [], {} , {}

        # print(len(self.exons), file=sys.stderr)
        for pos in range(len(self.exons)-1):
            intron=self.fasta[self.exons[pos].stop:self.exons[pos+1].start-1]
            if self.strand=="+":

                donor=self.exons[pos]
                acc=self.exons[pos+1]
                donor_start=donor.stop+1
                donor_end=donor_start+2
                acc_start=acc.start-2
                acc_stop=acc.start
                donor_fasta=self.fasta[donor_start-(length+1):donor_start+length-1]
                acc_fasta=self.fasta[acc_start-length-1:acc_start+length-1]

            else:
                intron=intron.reverse_complement()
                acc=self.exons[pos]
                donor=self.exons[pos+1]
                donor_start=donor.start-2
                acc_start=acc.stop+1

                donor_fasta=self.fasta[donor_start-length+1:donor_start+length+1]
                acc_fasta=self.fasta[acc_start-length+1:acc_start+length+1]
                donor_fasta=donor_fasta.reverse_complement()
                acc_fasta=acc_fasta.reverse_complement()

            acc_fasta.id=self.name+"_"+str(acc_start)
            donor_fasta.id=self.name+"_"+str(donor_start)
            intron.id=self.name+"."+str(pos)
            intron.description="Intron"
            acc_fasta.description="Acceptor"
            donor_fasta.description="Donor"
            #Only yield canonical sites.
            if not self.strict: donor_signal=("GT","GC")
            else: donor_signal=("GT",)

            donor_observed=str(donor_fasta[length:length+2].seq).upper()
            acceptor_observed=str(acc_fasta[length:length+2].seq)
            if (donor_observed in donor_signal and acceptor_observed=="AG"):
                donors[donor_start]=donor_fasta
                acceptors[acc_start]=acc_fasta
            else:
                self.logger.warn("Exons {0}-{1} have a non-canonical splice site: {donor}, {acceptor}".format(pos, pos+1,
                                                                                                         donor=donor_observed,
                                                                                                         acceptor=acceptor_observed))
                self.non_canonical=True
                if self.strict:
                    self.failed=True
                    return
                
            introns.append(intron)
            #Only if ALL splices
            self.introns=introns
            self.donors=donors
            self.acceptors=acceptors

    def calculateRegion(self):
        '''This function produces the genomic region of the gene, with the specified flank.'''
        if not self.start or not self.stop:
            self.logger.error("No Start or Stop information provided.")
            return
        start,stop=self.start.start, self.stop.stop
        start,stop=sorted([start,stop])
        start=max(start-self.flank,0)
        stop=min(stop+self.flank, len(self.fasta))
        seq=self.fasta[start:stop]

        ## If the genic region is on the - strand, we have to take the RC of the sequence.
        if self.strand=="-": seq=seq.reverse_complement()

        seq.id=self.name
        seq.description=''
        self.region=seq

    def calculateFalseSignals(self):
        '''This function has the purpose of finding all "false" signals inside the sequence. It returns:
           - false_starts
           - false_stops
           - false_donors
           - false_acceptors'''
        
        self.logger.info("Calculate all the false splicing signals.")
        if not self.region:
            self.calculateRegion()
        
        self.false_starts=[]
        self.false_stops=[]
        self.false_donors=[]
        self.false_acceptors=[]

        ##Donors and Acceptors
        self.logger.debug("Calculating splice false signals.")
        for signal, real, falsies in zip( ["GT","AG"], [self.donors, self.acceptors], [self.false_donors, self.false_acceptors]):
            pattern="(?=([ACGTNacgtn]{{{0}}}{1}[ACGTNacgtn]{{{2}}}))".format(self.length, signal, self.length-2 )
            pattern=re.compile(pattern) #Search all occurrences of the signal. Automatic generation!
            for hit in pattern.finditer(str(self.region.seq).upper()):
                    if hit.start()+self.length in real: continue #Exclude true hits
                    falsies.append(hit.groups()[0])

        #Starts and Stops
        self.logger.debug("Calculating Start and Stop false signals.")
        for signal, real, falsies in zip( ["ATG", "(TGA|TAG|TAA)"], [self.start, self.stop], [self.false_starts, self.false_stops]):
            pattern=re.compile("(?=([ACGTNacgtn]{{{0}}}{1}[ACGTNacgtn]{{{2}}}))".format(self.length, signal, self.length-3 )) #Search all occurrences of the signal. Automatic generation!
            for hit in pattern.finditer(str(self.region.seq).upper()):
                if hit.start()+self.length==real.start: continue
                falsies.append(hit.groups()[0])


    def empty(self):
        '''Mock method to set to None all the derived features - start, stop, splices, etc.'''
        self.logger.warn("Removing all features and going back to a blank state")
        self.logger.warn
        self.start=None
        self.stop=None
        self.introns=None
        self.donors=None
        self.acceptors=None
                           
        
    def calculateSignals(self, fasta, **kwargs):
        '''Method to calculate all the relevant signals. If one of them fails, it indicates which one and resets everything through the empty() method.'''

        if 'strict' in kwargs: self.strict=kwargs['strict']
        if 'length' in kwargs: self.length=int(kwargs['length']/2)
        elif not self.length: self.length=int(60/2)
        if 'flank' in kwargs:
            self.flank=kwargs['flank']
        else:
            self.flank=0
        self.fasta=fasta #Set this as a class parameter
        #Do I really need this?

        order=["calculateStart", "calculateStop",
               "calculateSplice",
               "calculateRegion",
               "calculateFalseSignals"]

        for method in order:
            getattr(self, method)() #Execute method :-)
            if self.failed:
                self.empty
                return
            
    def printSignals(self, fasta, **kwargs):
        '''This master function retrieves and yields the signal sequences
        so that they can be printed by programs which access the method.'''

        self.strict=False
        if 'length' in kwargs: self.length=kwargs['length']/2
        if 'strict' in kwargs: self.strict=True
        elif not self.length: self.length=60/2
        if 'flank' not in kwargs: kwargs['flank']=None
        self.calculateSignals(fasta)#, flank=kwargs['flank']) 
        self.calculateRegion(fasta)
        yield self.start_fasta.format('fasta')
        for donor,intron,acceptor in zip(self.donors, self.introns, self.acceptors):
            yield donor.format('fasta')
            yield intron.format('fasta')
            yield acceptor.format('fasta')
        yield self.stop_fasta.format('fasta')

    def printMatch(self):
        '''Function to transform the GTF record into a match/match_part-like GFF record.'''

        match_line="\t".join([ str(attribute) for attribute in
                               [self.chrom, self.source, "match_part", 
