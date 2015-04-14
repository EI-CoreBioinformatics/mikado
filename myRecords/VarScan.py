#!/usr/bin/env python3
#coding: utf_8

from myRecords import HeaderError,tabParser
import re

class varLine(tabParser):
    '''This class allows easy parsing of a Varscan line.
    Options:
        - type: either "SNP" or "Indel".

    Important fields:
       - chrom
       - pos
       - refBase
       - cons - the base consensus for this position.
       - reads1 - reads supporting the WT allele.
       - reads2 - reads supporting the mutation.
       - varFreq - Variation frequency.
       - qual1 - quality of mapping for reads supporting the WT allele.
       - qual2 - quality of mapping for reads supporting the mutant allele.
       - pValue - derived from a Fisher test on the nmber of reads.
       - var (mutazione osservata) - nucelotides only.
       - clas enum('SNP','deletion','insertion')
    '''
    def __init__(self,line):
        super(varLine, self).__init__(line)
        assert isinstance(line,str)
	if self._fields[0]=='Chrom': raise HeaderError
        self.chrom=self._fields[0]
        self.pos=int(self._fields[1])
        self.refBase=self._fields[2]
        assert len(self.refBase)==1
        self.cons=self._fields[3]
        #Understand the type of mutation we have here
        #Swear Swarm...VarScan does not always prepend the right + or - sign to the mutation in the var field.
        #This part has been put here, although it logically should go after, to minimize the resources used before the assesment is made.
        self.clas='SNP'
        self.legit=True
        if self.cons.find('+')!=-1: self.clas='insertion'
#        if re.findall(r"\+",self.cons): self.clas='insertion'
        if self.cons.find('-')!=-1: self.clas='deletion'
#        elif re.findall(r"\-",self.cons): self.clas='deletion'
        self.var=self._fields[-1]
        if self.clas=='SNP': assert len(self.var)==1
        if self.clas=='insertion' and self.var.find('+')==-1: self.legit=False; return #In these three cases, I have something strange.
        elif self.clas=='deletion' and self.var.find('-')==-1: self.legit=False; return
        elif self.clas=='SNP' and ( self.var.find('+')!=-1 or  self.var.find('-')!=-1 ): self.legit=False; return

        self.reads1=int(self._fields[4])
        self.reads2=int(self._fields[5])
        self.varFreq=float(re.sub(r",",'.',self._fields[6].rstrip('%')))
        self.strands1=int(self._fields[7])
        assert self.strands1 in (0,1,2)
        self.strands2=int(self._fields[8])
        assert self.strands2 in (0,1,2)
        self.qual1,self.qual2=tuple(int(i) for i in self._fields[9:11])
        self.pValue=float(self._fields[11])
        self.mapQual1,self.mapQual2,self.readsPlus1,\
        self.readsMinus1,self.readsPlus2,self.readsMinus2=tuple(int(i) for i in self._fields[12:18])

    def __str__(self): return '\t'.join(self._fields)

    def toGff(self):
        '''This easy function converts the Varscan line to a Gff line.'''
        if self.strands2==2: self._Strand='.'
        else:
                if self.readsPlus2>self.readsMinus2: self._Strand='+'
                else: self._Strand='-'
        gff=self.chrom,'Varscan',self.type,str(self.pos),str(self.pos+1),\
                str(round(self.pValue,2)),self._Strand,'.','.'
        return '\t'.join(gff)
