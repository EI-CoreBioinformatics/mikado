#!/usr/bin/env python3
#coding: utf_8


import re
import sys,subprocess
from math import log
import io

class Blat(object):
    '''This class functions as interface between the main program and the blatLine object. It implements the classic iterator functions (__iter__ and next).'''
    def __init__(self,handle,unique=False):
        '''The class must be instantiated with a file handle. The kwarg "unique" is used to indicate whether the program should keep only the first hit 
        for each query, or not. Default: False.'''
        if isinstance(handle,str): self.__handle=open(handle)
        elif isinstance(handle, file): self.__handle=handle
        self.__curr=''
        self.__unique=unique

    def __next__(self):
        while True:
            line=self.__handle.readline()
            self.__bestScore=float('-inf')
            if line=='':
                raise StopIteration
            try:
                self.record=blatLine(line)
                if self.__unique:
                    if self.__curr!=self.record.qName:
                        self.__curr=self.record.qName
                        break
                    else:
                        continue
                else:
                    break
            except:
                continue
        return self.record
    def __iter__(self):
        return self

class blatLine(object):
    '''This class defines a typical BLAT psl line. It returns the various fields and other derived arguments (e.g. identity,coverage,offset).
       The functions here contained apart from __init__ have been translated from C, the original ones can be found at the official UCSC page: http://genome.ucsc.edu/FAQ/FAQblat#blat4
    
    Each record contains the following attributes:
    line: the original BLAT line.
    matches - number of bases which match. (int)
    mismatches - number of aligned bases with mismatch. (int)
    rep - number of aligned bases in repetitive regions. (int)
    ns - number of Ns in the alignment. (int)
    qNumInserts- number of gaps opened in the query sequence. (int)
    qBaseInserts - number of bases in opened gaps in the query sequence. (int)
    tNumInserts- number of gaps opened in the target sequence. (int)
    tBaseInserts - number of bases in opened gaps in the target sequence. (int)
    strand - strandness. Must be either + or -. (string)
    qName - query name. (string)
    qSize - query size. (int)
    qStart - start of the alignment on the query sequence. (int)
    qEnd - end of the alignment on the query sequence. (int)
    tName - target name. (string)
        tSize - target size. (int)
        tStart - start of the alignment on the target sequence. (int)
        tEnd - end of the alignment on the target sequence. (int)
    blocks - number of blocks the alignment is broken up into. For mRNA, equivalent to number of exons. (int).
    blockSizes - list of the sizes of the various blocks. List of ints.
    qStarts - list of the positions at which the various blocks start on the query sequence. List of ints.
    qEnd - list of the positions at which the various blocks end on the query sequence. List of ints.
    tStarts - list of the positions at which the various blocks start on the target sequence. List of ints.
    tEnds - list of the positions at which the various blocks end on the target sequence. List of ints.
    identity - identity percentage. (Float)
    coverage - coverage percentage (Float)
    alignment - length of the alignment (minimum between the query and the target alignment). (int)
    offset - difference between the aligned region and the query size (qSize-Alig). (int)
    Score - the alignment score.
    '''
    def __init__(self,line):
        self.line=line.rstrip()
        fields=self.line.split('\t')
        self.matches,self.mismatches,self.rep,self.ns,self.qNumInserts,self.qBaseInserts,self.tNumInserts,self.tBaseInserts=(int(i) for i in fields[0:8])
        self.strand,self.qName=fields[8:10]
        assert self.strand in ('+','-')
        self.qSize,self.qStart,self.qEnd=[int(i) for i in fields[10:13]]
        self.tName=fields[13]
        self.tSize,self.tStart,self.tEnd=[int(i) for i in fields[14:17]]
        self.blocks=int(fields[17])
        self.blockSizes=[int(i) for i in fields[18].split(',')[:-1]]
        self.qStarts=[int(i) for i in fields[19].split(',')[:-1]]
        self.qEnds=[x+y for x,y in zip(self.blockSizes,self.qStarts)]
        self.tStarts=[int(i) for i in fields[20].split(',')[:-1]]
        self.tEnds=[x+y for x,y in zip(self.blockSizes,self.tStarts)]
        self.identity=round(100-self.calcMilliBad()*0.1)
        self.coverage=round((self.alignment)/self.qSize*100)
        self.offset=self.qSize-self.alignment
        self.Score=self.pslScore()

    def isProt(self):
        '''This function returns a boolean; it determines whether the input was protein vs. DNA or not.'''
        if self.strand=='+':
            return self.tEnd==self.tStarts[-1]+3*self.blockSizes[-1]
        else:
            return self.tStart==self.tSize-self.tStarts[-1]+3*self.blockSizes[-1]


    def calcMilliBad(self):
        '''This function calculates the identity score.'''
        if self.isProt():
            sizeMul=3

            self.insertFactor=self.tNumInserts
        else:
            sizeMul=1
            self.insertFactor=self.qNumInserts

        self.qAlig=self.qEnd-self.qStart
        self.tAlig=self.qEnd-self.qStart
        self.alignment=min(self.qAlig,self.tAlig)
        if self.alignment==0:
            return 0

        sizeDif=self.qAlig-self.tAlig
        if sizeDif<0:
            if sizeMul==1:
                sizeDif=0
            elif sizeMul==3:
                sizeDif=-sizeDif

        total=sizeMul*(self.matches + self.rep + self.mismatches)
        if total!=0:
            milliBad=(1000*(self.mismatches*sizeMul+self.insertFactor+round(3*log(1+sizeDif))))/total

        return milliBad

    def pslScore(self):
        '''This function calculates the base score.'''

        if self.isProt(): sizeMul=3
        else: sizeMul=1

        score=sizeMul*(self.matches+(self.rep>>1)) - sizeMul*self.mismatches- self.qNumInserts-self.tNumInserts #Il >> è un bitwise shift.

        return score
