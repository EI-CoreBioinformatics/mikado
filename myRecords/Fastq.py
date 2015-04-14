#coding: utf_8


import re
import sys,subprocess
from math import log
import io

def parseFastq(handle,enc):
    '''A fast (I hope) parsing function for the Fastq record.
    Type supported:
    illumina
    illumina-old
    sanger'''
    enc=enc.lower()
    assert enc in ('illumina','illumina_old','sanger')
    while True:
        try:line1=handle.readline()
        except:
            try: handle=open(handle); line1=handle.readline()
            except: raise
        if line1=='': raise StopIteration
        while line1[0]!='@': line1=handle.readline()
        line2=handle.readline()
        line3=handle.readline()
        line4=handle.readline()
        yield fastqRecord(line1+line2+line3+line4,encoding=enc)

class fastqRecord(object):
    '''This dummy class provides a fast access to a Fastq record.
    The arguments of each istance are:
        - id; the name of the record.
        - seq; the read sequence.
        - qual; the read qualities.'''
    def __init__(self,block,encoding='illumina'):
        self._block=io.StringIO()
        self._enc=encoding
        assert self._enc in ('illumina','illumina_old','sanger'), \
        'The quality must be chosen among "illumina", "illumina_old", "sanger".'
        self._block.write(block); self._block.seek(0)
        try:
            self.id=self._block.readline().lstrip('@').rstrip()
            self.seq=self._block.readline().rstrip()
            for letter in self.seq: assert letter.upper() in ('A','C','G','T','N')
            self._block.readline();
            self.qual=self._block.readline().rstrip()
            assert len(self.qual)==len(self.seq)
        except: raise ValueError
        self._block.write(block); self._block.seek(0)

    def __str__(self):
        return '@'+self.id+'\n'+ \
            self.seq + '\n' + \
            '+' + '\n'+ \
            self.qual

    def __len__(self): return len(self.seq)

    def toSam(self):
        flag=12

        pass #Remember to implement this!

    def qualList(self):
        '''This method provides a list with the qualities.'''
        self.quals=[]
        types={'illumina': (64,40),\
               'illumina_old': (33,40),\
               'sanger': (64,62)}

        for letter in self.qual:
            q=ord(letter)-types[self._enc][0]
            assert q<= types[self._enc][1], 'Encountered too large a value in {0}.'.format(self.id)
            assert q>0, 'Encountered too little a value in {0}'.format(self.id)
            self.quals.append(q)


    def qualZip(self):
        '''This method returns a list of tuples, of the form (nucleotide, numerical quality.)'''#ascii quality, numerical quality.)'''
        dummy=[]
        for letter in self.seq: dummy.append(letter)
        dummy2=[]
        for letter in self.qual: dummy2.append(letter)
        self.qualList()
        return list(zip(dummy,self.quals))
