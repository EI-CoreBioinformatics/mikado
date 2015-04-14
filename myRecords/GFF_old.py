#!/usr/bin/env python3
#coding: utf_8


class gffLine(object):

    def __init__(self,line):
        
        self._line=line
        self._fields=line.rstrip().split('\t')
        self.chrom,self.source,self.feature=self._fields[0:3]
        self.start,self.end=tuple(int(i) for i in self._fields[3:5])

        if self._fields[5]=='.': self.score=None
        else: self.score=float(self._fields[5])

        self.strand=self._fields[6]
        if self.strand=='.': self.strand=None
        assert self.strand in (None,'+','-','?')

        if self._fields[7]=='.': self.phase=None
        else: 
            try: 
                self.phase=int(self._fields[7]); assert self.phase in (0,1,2)
            except: raise

        self._Attr=self._fields[8]
        self.attributes={}

        for item in self._Attr.rstrip().split(';'):
            item=item.split('=')
            self.attributes[item[0]]=item[1]
        assert ('ID' in self.attributes) or ('Parent' in self.attributes)
        tags=['Name','ID','Parent']
        if 'ID' in self.attributes: self.id=self.attributes['ID']
        else: self.id=self.attributes['Parent']

        for tag in tags:
            if tag in self.attributes: self.__dict__[tag.lower()]=self.attributes[tag]

    def __str__(self): return '\t'.join(self._fields)


class GFF3(object):
    def __init__(self,handle):
        if isinstance(handle,file):
            self._handle=handle

        else:
            assert isinstance(handle,str)
            try: self._handle=open(handle)
            except: raise ValueError('File not found: {0}'.format(handle))

        self.header=''

    def __iter__(self): return self

    def __next__(self):
        line=self._handle.readline()
        if line=='': raise StopIteration
        while line[0]=='#':
            self.header+=line
            line=self._handle.readline()

        return gffLine(line)
