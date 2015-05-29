import random
try:
    import Bio.File
except:
    pass
from loci_objects import Parser
from sqlalchemy import Column,Index,String,Integer,CHAR,Boolean #Basic constructs for the tables                                                                                                               
from sqlalchemy.ext.declarative import declarative_base #Necessary to create the metadata                                                                                                                    
 
 
Base=declarative_base()

class bed12:
    
    __tablename__="bed12"
     
    id=Column(Integer, primary_key=True)
    chrom=Column(String)
    start=Column(Integer)
    end=Column(Integer)
    name=Column(String)
    strand=Column(CHAR)
    thickStart=Column(Integer)
    thickEnd=Column(Integer)
    score=Column(Integer)
    hasStart=Column(Boolean)
    hasStop=Column(Boolean)
    __table_args__ = ( Index("index", "chrom", "start","end", "thickStart", "thickEnd"  ) , Index("nameindex", "name", "thickStart", "thickEnd"))
    
    
    def __init__(self, line, fasta_index = None, transcriptomic=True):
        if type(line) is str:
            if line[0]=="#":
                self.header=True
                return
            
            line=line.rstrip().split("\t")
        elif type(line) not in (list, tuple):
            raise TypeError("I need an ordered array, not {0}".format(type(line)))
        if len(line)!=12:
            self.header=True
            return
            #raise ValueError("Erroneous number of fields detected")
        
        self.transcriptomic = transcriptomic
        self.header=False
        self.chrom, self.start, self.end, \
            self.name, self.score, self.strand, \
            self.thickStart, self.thickEnd, self.rgb, \
            self.blockCount, self.blockSizes, self.blockStarts = line
            
        self.start=int(self.start)+1
        self.end = int(self.end)
        self.score = float(self.score)
        self.thickStart = int(self.thickStart)+1
        self.thickEnd = int(self.thickEnd)
        self.blockCount = int(self.blockCount)
        self.blockSizes = [int(x) for x in self.blockSizes.split(",")]
        self.blockStarts = [int(x) for x in self.blockStarts.split(",")]
        self.has_start_codon = False
        self.has_stop_codon = False
        
        if fasta_index is not None:
            assert self.id in fasta_index
            start_codon = fasta_index[self.id][self.thickStart-1:self.thickStart+2] # I have translated into 1-offset
            stop_codon = fasta_index[self.id][self.thickEnd-3:self.thickEnd]
            if self.strand=="-":
                start_codon=start_codon.reverse_complement()
                stop_codon=stop_codon.reverse_complement()
                start_codon,stop_codon=stop_codon,start_codon
            if str(start_codon.seq)=="ATG":
                self.has_start_codon=True
            if str(stop_codon.seq) in ("TAA", "TGA", "TAG"):
                self.has_stop_codon=True
        
        assert self.blockCount==len(self.blockStarts)==len(self.blockSizes)
        
        
    def __str__(self):
        
        line = [self.chrom, self.start-1, self.end, self.name, self.score]
        if self.strand is None:
            line.append(".")
        else:
            line.append(self.strand)
        line.extend( [self.score, self.thickStart-1, self.thickEnd, self.blockCount] )
        line.append( ",".join([str(x) for x in self.blockSizes]  ) )
        line.append( ",".join([str(x) for x in self.blockStarts]  ) )
        return "\t".join([str(x) for x in line])
        
    def __eq__(self,other):
        for key in ["chrom","strand","start","end","thickStart","thickEnd","blockCount","blockSizes","blockStarts"]:
            if getattr(self,key)!=getattr(other,key): return False
        return True
    
    def __hash__(self):
        return super().__hash__()
    
    
    @property
    def strand(self):
        return self.__strand
    
    @strand.setter
    def strand(self,strand):
        if strand in (".","?"):
            self.__strand = None
        elif strand in ("+", "-"):
            self.__strand = strand
        else:
            raise ValueError("Erroneous strand provided: {0}".format(self.strand))
        
    @property
    def cds_len(self):
        return self.thickEnd-self.thickStart+1
    
    @property
    def has_start_codon(self):
        return self.__has_start
    
    @has_start_codon.setter
    def has_start_codon(self,value):
        if value not in (None,True,False):
            raise ValueError()
        self.__has_start = value 

    @property
    def has_stop_codon(self):
        return self.__has_stop
    
    @has_stop_codon.setter
    def has_stop_codon(self,value):
        if value not in (None,True,False):
            raise ValueError()
        self.__has_stop = value 

    @property
    def full_orf(self):
        return self.has_stop_codon and self.has_start_codon

    @property
    def id(self):
        if self.transcriptomic:
            return self.chrom
        else:
            return self.name

class BED12(Parser):

    '''Parser class for a BED12 file. It accepts optionally a fasta index which is used to determine whether an ORF has start/stop codons.'''
    
    def __init__(self,handle, fasta_index=None):
        
        super().__init__(handle)

        if type(fasta_index) is dict:
                #check that this is a bona fide dictionary ...
            assert type( fasta_index[random.sample(fasta_index.keys(),1)] ) is Bio.SeqRecord.SeqRecord
        elif fasta_index is not None:
            assert type(fasta_index) in (Bio.File._IndexedSeqFileDict, Bio.File._SQLiteManySeqFilesDict), (type(fasta_index))

        self.fasta_index = fasta_index

        self.header=False

    def __iter__(self): return self

    def __next__(self):
        line=self._handle.readline()
        if line=='': raise StopIteration
        return bed12(line, fasta_index=self.fasta_index)
    