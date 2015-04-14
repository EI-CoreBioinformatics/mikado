#!/usr/bin/env python3
#coding: utf_8


import re
import sys,subprocess
from math import log
import io

'''Deprecated. Use the pysam module instead.'''


def parseSAM(handle):
    for line in handle:
        yield AlignedRead(line)

def parseBam(handle):
    handle.close()
    bam=subprocess.Popen(['samtools','view',handle.name],stdout=subprocess.PIPE)
    for line in bam.stdout: yield AlignedRead(line)

class AlignedRead(object):
    '''
    This in-house class defines a SAM/BAM read line. Arguments of the object (following the convention established at the SAM format specification):
         qname: query name
         flag: bitwise flag
         rname: reference sequence name
         pos: 1-based leftmost mapping position
         mapq: mapping quality
         cigar: cigar term.
         rnext: reference name of the mate/next fragment.
         tlen: template length. Negative if the read is reverse.
         seq: fragment sequence.
         qual: fragment Quality. It is by default in the Sanger FASTQ format (ASCII+33)
    
    Following, for reference, the flags as defined by the SAM format:
          
               0x1       template having multiple fragments in sequencing  
               0x2       each fragment properly aligned according to the aligner 
               0x4       fragment unmapped 
               0x8       next fragment in the template unmapped 
               0x10      SEQ being reverse complemented 
               0x20      SEQ of the next fragment in the template being reversed 
               0x40      the first fragment in the template 
               0x80      the last fragment in the template 
               0x100     secondary alignment 
               0x200     not passing quality controls 
               0x400     PCR or optical duplicate

    Arguments derived from the bitwise flag:
        Each argument is a boolean. If set to None, it indicates that no information can be reliably inferred from the read.

        paired: The read is in a pair.
        proper_pair: Both fragments are mapped.
        is_mapped: the read is properly mapped.
        mate_unmapped: the mate of the read is not properly mapped.
        reverse: the read is aligned on the reverse strand.
        mate_reverse: the mate of the read is aligned on the reverse strand.
        first: the read is the first in a pair.
        second: the read is the second in a pair.
        secondary: this is not the primary alignment.
        qc_fail: the read failed the quality control test.
        duplicate: the read is a PCR duplicate.
'''


    def __init__(self,line):
        self.line=line
        self._fields=line.rstrip().split('\t')


        self.qname=self._fields[0]
        self.flag=int(self._fields[1])#(bin(self._fields[1]))
        self.rname=self._fields[2];
        self.pos,self.mapq=tuple(int(i) for i in self._fields[3:5]);
        if self.mapq==255: self.mapq=None;
        self.cigar,self.rnext=self._fields[5:7];
        if self.rnext=='=': self.rnext=self.rname
        self.pnext,self.tlen=tuple(int(i) for i in self._fields[7:9])
        self.seq,self.qual=self._fields[9:11]
        self._tags=self._fields[11:]

        #Blocco di informazioni che possono essere ricavate dal flag

        self.paired=False
        self.proper_pair=None #Li inizializzo come None
        self.mate_unmapped=None
        self.mate_reverse=None
        self.first=None
        self.second=None
        self.qc_fail=False #Quality control
        self.duplicate=False

        if self.flag & 0x1:
           self.paired=True

           if self.flag & 0x2: self.proper_pair=True
           else: self.proper_pair=False

           if not(self.flag & 0x8): self.mate_unmapped=False
           else:
               self.mate_unmapped=True
               self.rnext=None

           if self.flag & 0x20: self.mate_reverse=True
           else: self.mate_reverse=False

           if self.flag & 0x40: self.first=True
           else: self.first=False

           if self.flag & 0x80: self.second=True
           else: self.second=False

        self.is_mapped=True
        self.reverse=False
        self.secondary=False
        if self.flag & 0x4:
           self.is_mapped=False
           self.proper_pair=None
           self.reverse=None
           self.secondary=None
        else:
           if self.flag & 0x10: self.reverse=True
           if self.flag & 0x100: self.secondary=True

        if self.flag & 0x200: self.qc_fail=True
        if self.flag & 0x400: self.duplicate=True


        #Fine blocco flags
        if self.pnext==0: self.rnext=None; self.mate_reverse=None #Secondo l'indicazione del SAM format
        self.tags=dict()

        for tag in self._tags:
                        sTag=tag.split(':')
                        if sTag[1]=='i':
                                self.tags[sTag[0]]=int(sTag[2])
                        elif sTag[1]=='f':
                                self.tags[sTag[0]]=float(sTag[2])
                        elif sTag[1]=='H':
                                self.tags[sTag[0]]=int(sTag[2],16)
                        elif sTag[1] in ('A','Z'):
                                if sTag[1]=='A': assert len(sTag[2])==1
                                self.tags[sTag[0]]=sTag[2]

    def lenAlig(self):
        '''This function serves the only purpose of calculating the number of bases actually mapped.'''
        if not self.is_mapped: return 0
        d={}
        n=''
        for i in self.cigar:
            try: int(i)
            except:
                if i not in d: d[i]=0
                d[i]+=int(n)
                n=''
            else: n+=i
        lun=0
        cigars=('M','I')
        for l in cigars:
            if l in d: lun+=d[l]
        rlun=0
        rcigars=('S','=','X')
        for l in rcigars:
            if l in d: rlun+=d[l]
        if (rlun+lun)!=len(self.seq):
            print(d,len(self.seq),file=sys.stderr)
            raise AssertionError
        return lun

    def __str__(self):
        '''Return the original SAM line.'''
        return self.line.rstrip()

    def __reverse(self,string): #Reverse a string.
         s=str('')
         assert type(string)==type(str())
         for i in range(len(string)-1,-1,-1):
             try: s+=str(string[i])
             except: raise
         return s

    def __revFa(self,string): #Reverse-complement a Fasta sequence.
         string=self.__reverse(string.upper())
         s=''
         for i in string:
            assert i in ('A','C','G','T','N')
            if i=='A': i='T'
            elif i=='T': i='A'
            elif i=='C': i='G'
            elif i=='G': i='C'
            s+=i
         return s

    def __len__(self):
        '''Returns the length of the aligned sequence.'''
        return len(self.seq)

    def toFastq(self):
         '''This method returns the Fastq-Sanger representation of the read.'''
         suff='/1\n'
         if self.second==True: suff='/2\n'
         if self.reverse:
             seq=self.__revFa(self.seq)
             qual=self.__reverse(self.qual)
         else:
            seq=self.seq;
            qual=self.qual;
         assert len(seq)==len(qual);
         fastq='@'+self.qname+suff + \
             seq + '\n' + \
             '+'+ '\n'+\
             qual
         return fastq

    def toFasta(self):
         '''This method returns the Fasta representation of the read.'''
         suff='/1\n'
         if self.second==True: suff='/2\n'
         if self.reverse: seq=self.__revFa(self.seq)
         else: seq=self.seq
         fasta='>'+self.qname+suff + \
               seq
         return fasta
