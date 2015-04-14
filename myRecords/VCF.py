#!/usr/bin/env python3


from myRecords import HeaderError
from sqlalchemy import Column, Numeric, String, Integer, create_engine, Index, UniqueConstraint
from sqlalchemy.orm import *
import sys,itertools
from sqlalchemy.ext.declarative import declarative_base

vcfBase=declarative_base()

class VCFparser(object):
	def __init__(self,handle):
		if not isinstance(handle,file):
			try: handle=open(handle)
			except: raise TypeError
		self._handle=handle
		self.header=''
		pos=self._handle.tell(); line=self._handle.readline()
		while line[0]=='#': self.header+=line; pos=self._handle.tell(); line=self._handle.readline(); #Iterate until we get out of the comment section
		self._handle.seek(pos)

	def __iter__(self): return self

	def __next__(self): 
		line=self._handle.readline()
		if line=='': raise StopIteration
		return VCF(line)
		
		

class VCF(vcfBase):

	__tablename__ = 'dbsnp'
	#unid = Column(Integer, autoincrement=True, primary_key=True)
	chrom = Column(String(10), autoincrement=False, primary_key=True)
	pos = Column(Integer, autoincrement=False, primary_key=True)
	rs = Column(String(20), default=None, primary_key=True)
	ref = Column(String(50), default=None, nullable=False, primary_key=True)
	alt = Column(String(50), default=None, nullable=False, primary_key=True)
	
	__table_args__ = (
		UniqueConstraint('pos','alt','rs', name='uix1'),
	{} )

	def __init__(self,line):
		try: 
			if line[0]=="#": raise HeaderError
		except: self.line=line.rstrip(); return None
		self.line=line.rstrip()
		self._columns=line.rstrip().split('\t')
		self.chrom=self._columns[0]

		try: self.pos=int(self._columns[1])
		except: print(self.line); raise

		if self._columns[2]=='.': self.rs=None
		else: self.rs=self._columns[2]

		self.ref,self.alt=self._columns[3:5]

		try: self.qual=float(self._columns[5])
		except: self.qual=None

		if self._columns[6]=='.': self.Filter=None
		else: self.Filter=self._columns[6]
		self.info=dict()
		for x in self._columns[7].split(';'):
			x=x.split('=')
			if len(x)==1: self.info[x[0]]=True
			else: 
				try: self.info[x[0]]=int(x[1])
				except: 
					try: self.info[x[0]]=float(x[1])
					except: self.info[x[0]]=x[1]


		self._genotype=dict()
		#for x in itertools.izip(self._columns[8].split(':'),self._columns[9].split(':')):
			#GT:AD:DP:GQ:PL  0/1:4,2:6:61.95:62,0,142
		#	self._genotype[x[0]]=x[1]

		self.hom_ref=False
		self.het=False
		self.hom_alt=False
		self.hom_alt_probability=0
		self.het_probability=0
		self.hom_ref_probability=0
		self.depth=0
		self.genq=0
		self.refCount=0
		self.altCount=0

		#if self._genotype['GT']=='0/0': self.hom_ref=True
		#elif self._genotype['GT']=='0/1': self.het=True
		#elif self._genotype['GT']=='1/1': self.hom_alt=True

		if 'DP' in self._genotype:
			self._genotype['DP']=int(self._genotype['DP'])
			self.depth=self._genotype['DP']

		if 'GQ' in self._genotype:
			self._genotype['GQ']=float(self._genotype['GQ'])
			self.genq=self._genotype['GQ']

		if 'AD' in self._genotype:
			self._genotype['AD']=tuple(int(i) for i in self._genotype['AD'].split(','))
			self.refCount=self._genotype['AD'][0]
			self.altCount=self._genotype['AD'][1]

		if 'PL' in self._genotype:
			self._genotype['PL']=tuple(float(i) for i in self._genotype['PL'].split(','))
			self.hom_ref_probability=10**(-(self._genotype['PL'][0]/10))
			self.het_probability=10**(-(self._genotype['PL'][1]/10))
			self.hom_alt_probability=10**(-(self._genotype['PL'][2]/10))

	def __str__(self): return self.line


if __name__=='__main__':
	vcf=open('s_1_1_chr21.sort.chr21.realigned.bam.snp.vcf')
	vcfPar=VCFparser(vcf)

	v=next(vcfPar)
	print(v)

	print(v.hom_ref_probability,v.het_probability,v.hom_alt_probability,v.refCount,v.altCount,v.depth,sep='\t')
	print(v._genotype)
	print(v.info)
#	print(v.__str__().split('\t')[-1])
