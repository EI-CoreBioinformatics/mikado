#!/usr/bin/env python3
#coding: utf_8

from myRecords import Parser,HeaderError,tabParser,SizeError
import re

class cuffFile(Parser):
    def __init__(self,handle):
        super(cuffFile, self).__init__(handle)
        self.transcript=True

    def __next__(self):
        while True:
            line=self._handle.readline()
            if line=='': raise StopIteration
            if self.transcript:
                try: return cuffTrans(line)
                except SizeError:
                    self.transcript=False
                    try: return cuffGene(line)
                    except: raise

                except HeaderError: continue

            else:
                try: return cuffGene(line)
                except: raise

class cuffLine(tabParser):
    def __init__(self,line):
        super(cuffLine, self).__init__(line)
        if self._fields[0] in ('gene_id','trans_id'): raise HeaderError
        self.id=self._fields[0]
        self.bundle_id=int(self._fields[1])
        self.chrom=self._fields[2]
        self.start,self.stop=tuple(int(i) for i in self._fields[3:5])
        self.fpkm=float(self._fields[5])
        self.status=self._fields[-1]

class cuffGene(cuffLine):
    def __init__(self,line):
        super(cuffGene, self).__init__(line)
        self.type='gene'
        self.fpkm_lo,self.fpkm_hi=tuple(float(i) for i in self._fields[6:8])

    def __str__(self):
        return '\t'.join([str(i) for i in [self.id, self.bundle_id,
                                           self.chrom, self.start, self.stop, self.fpkm,
                                           self.fpkm_lo, self.fpkm_hi, self.status]])

class cuffTrans(cuffLine):
    def __init__(self,line):
        super(cuffTrans, self).__init__(line)
        if len(self._fields)<12: raise SizeError
        self.type='transcript'
        self.fpkm_lo,self.fpkm_hi,self.coverage=tuple(float(i) for i in self._fields[8:11])
        self.length=int(self._fields[11])
        self.eff_len=float(self._fields[12])

    def __str__(self):
        return '\t'.join([str(i) for i in [self.id, self.bundle_id,
                                           self.chrom, self.start, self.stop, self.fpkm,
                                           self.fpkm_lo, self.fpkm_hi, self.coverage,
                                           self.length, self.eff_len, self.status]])

class Loci(Parser):
    def __init__(self,handle): super(Loci, self).__init__(handle)

    def __next__(self):
        try: return LociLine(self._handle.readline())
        except: raise

class LociLine(tabParser):
    def __init__(self,line):
        super(LociLine, self).__init__(line)
        self.gene=self._fields[0]
        try: self.chrom,self.strand,self._start_stop=re.split(r'[\[, \]]',self._fields[1])
        except ValueError: raise ValueError(self._fields[1])
        self.start,self.stop=[int(i) for i in self._start_stop.split('-')]
        self.name=self._fields[2].split('|')[0]
        if self.name=='-': self.name=None
