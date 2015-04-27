from builtins import str

class bed12:
    
    def __init__(self, line):
        if type(line) is str:
            if line[0]=="#":
                self.header=True
                return
            
            line=line.rstrip().split("\t")
        elif type(line) not in (list, tuple):
            raise TypeError("I need an ordered array, not {0}".format(type(line)))
        if len(line)!=12:
            raise ValueError("Erroneous number of fields detected")
        
        self.header=False
        self.chrom, self.start, self.end, \
            self.name, self.score, self.strand, \
            self.cdsStart, self.cdsEnd, self.rgb, \
            self.blockCount, self.blockSizes, self.blockStarts = line
            
        self.start=int(self.start)+1
        self.end = int(self.end)
        self.score = float(self.score)
        self.cdsStart = int(self.cdsStart)+1
        self.cdsEnd = int(self.cdsEnd)
        self.blockCount = int(self.blockCount)
        self.blockSizes = [int(x) for x in self.blockSizes.split(",")]
        self.blockStarts = [int(x) for x in self.blockStarts.split(",")]
        assert self.blockCount==len(self.blockStarts)==len(self.blockSizes)
        
        
    def __str__(self):
        
        line = [self.chrom, self.start-1, self.end, self.name, self.score]
        if self.strand is None:
            line.append(".")
        else:
            line.append(self.strand)
        line.extend( [self.score, self.cdsStart-1, self.cdsEnd, self.blockCount] )
        line.append( ",".join([str(x) for x in self.blockSizes]  ) )
        line.append( ",".join([str(x) for x in self.blockStarts]  ) )
        return "\t".join(line)
        
        
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
        return self.cdsEnd-self.cdsStart+1