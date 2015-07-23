import random
import os
from Bio import SeqIO
try:
    import Bio.File
except:
    pass
from mikado_lib.parsers import Parser

'''Generic module for parsing bed12Parser files.'''
 
class BED12:
    
    def __init__(self, *args:str, fasta_index = None, transcriptomic=False):
        if len(args)==0:
            self.header=True
            return 
        self.invalid = False
        line=args[0]
        if type(line) in (str,None):
            if line is None or line[0]=="#":
                self.header=True
                return
            
            line=line.rstrip().split("\t")
        elif line is None:
            self.header=True
            return
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
        self.has_start_codon = None
        self.has_stop_codon = None
        
        if fasta_index is not None:
            assert self.id in fasta_index
            if len(fasta_index[self.id])<self.length:
                self.invalid = True
            
            start_codon = fasta_index[self.id][self.thickStart-1:self.thickStart+2] # I have translated into 1-offset
            stop_codon = fasta_index[self.id][self.thickEnd-3:self.thickEnd]
            if self.strand=="-":
                start_codon=start_codon.reverse_complement()
                stop_codon=stop_codon.reverse_complement()
                start_codon,stop_codon=stop_codon,start_codon
            if str(start_codon.seq)=="ATG":
                self.has_start_codon=True
            else:
                self.has_start_codon=False
            if str(stop_codon.seq) in ("TAA", "TGA", "TAG"):
                self.has_stop_codon=True
            else:
                self.has_stop_codon=False
                
            #This is a bug in TD by which sometimes truncates ORFs even when there would be a read through
            if self.has_stop_codon is False and (self.strand!="-" and self.thickEnd<self.length-2) or (self.strand=="-" and self.thickStart>2):
                sequence = fasta_index[self.id]
                if self.strand!="-":
                    for num in range(self.thickEnd+3, self.end, 3  ):
                        codon= sequence[num-3:num]
                        if str(codon) in  ("TAA", "TGA", "TAG"):
                            self.has_stop_codon=True
                            break
                    self.thickEnd = num
                else:
                    
                    for num in reversed(range(self.start, self.thickStart,3)):
                        codon = sequence[num-3:num]
                        if str(codon) in  ("TTA", "TCA", "CTA"): #Reversed version, save on reversal, should save time
                            self.has_stop_codon=True
                            break
                    self.thickStart=num
        
        assert self.blockCount==len(self.blockStarts)==len(self.blockSizes)
        
        
    def __str__(self):
        
        line = [self.chrom, self.start-1, self.end, self.name, self.score]
        if self.strand is None:
            line.append(".")
        else:
            line.append(self.strand)
        line.extend( [self.thickStart-1, self.thickEnd, self.rgb, self.blockCount] )
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
    def strand(self,strand:str):
        if strand in (".","?"):
            self.__strand = None
        elif strand in ("+", "-"):
            self.__strand = strand
        else:
            raise ValueError("Erroneous strand provided: {0}".format(self.strand))
        
    @property
    def cds_len(self):
        return self.thickEnd-self.thickStart+1
    
    def length(self):
        return self.end-self.start
    
    @property
    def has_start_codon(self):
        return self.__has_start
    
    @has_start_codon.setter
    def has_start_codon(self,value:bool):
        if value not in (None,True,False):
            raise ValueError()
        self.__has_start = value 

    @property
    def has_stop_codon(self):
        return self.__has_stop
    
    @has_stop_codon.setter
    def has_stop_codon(self,value:bool):
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
        
    @property
    def invalid(self):
        return self.__invalid
    
    @invalid.setter
    def invalid(self, value):
        if type(value) is not bool:
            raise ValueError("Invalid value for invalid (should be boolean): {0}".format(value))
        self.__invalid = value
        

class bed12Parser(Parser):

    '''Parser class for a bed12Parser file. It accepts optionally a fasta index which is used to determine whether an ORF has start/stop codons.'''
    
    def __init__(self,handle, fasta_index=None, transcriptomic=False):
        
        super().__init__(handle)
        self.transcriptomic=transcriptomic

        if type(fasta_index) is dict:
                #check that this is a bona fide dictionary ...
            assert type( fasta_index[random.sample(fasta_index.keys(),1)] ) is Bio.SeqRecord.SeqRecord
        elif fasta_index is not None:
            if type(fasta_index) is str:
                assert os.path.exists(fasta_index)
                fasta_index = SeqIO.index(fasta_index, "fasta")
            else:
                assert type(fasta_index) in (Bio.File._IndexedSeqFileDict, Bio.File._SQLiteManySeqFilesDict), (type(fasta_index))

        self.fasta_index = fasta_index

        self.header=False

    def __iter__(self): return self

    def __next__(self):
        line=self._handle.readline()
        if line=='': raise StopIteration
        return BED12(line, fasta_index=self.fasta_index, transcriptomic=self.transcriptomic)
    