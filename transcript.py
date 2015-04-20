import operator
import re

class transcript:
    
    def __init__(self, gffLine):
        
        '''Initialise the transcript object, using a mRNA/transcript line.
        Note: I am assuming that the input line is an object from my own "GFF" class.
        The transcript instance must be initialised by a "(m|r|lnc|whatever)RNA" or "transcript" gffLine.'''
        
#         if gffLine.feature!="transcript" or "RNA" not in gffLine.feature.upper():
#             raise AttributeError("Wrong feature line provided:\n{0}\n{1}".format( gffLine.feature, gffLine ))
        
        self.chrom = gffLine.chrom
        self.feature="transcript"
        self.id = gffLine.attributes["ID"]
        self.start=gffLine.start
        self.strand = gffLine.strand
        self.end=gffLine.end
        self.exons = []
        self.junctions = []
        self.splices = []
        self.monoexonic = False
        self.finalized = False # Flag. We do not want to repeat the finalising more than once.
        self.parent = gffLine.parent
        self.attributes = gffLine.attributes
        
    def __str__(self):
        '''Each transcript will be printed out in the GFF style.
        This is pretty rudimentary, as the class does not hold any information on the original source, feature, score, etc.'''
        
        self.finalize() #Necessary to sort the exons
        attr_field = "ID={0}".format(self.id)
        if self.parent is not None:
            attr_field = "{0};Parent={1}".format(attr_field, self.parent)
        if self.strand is None:
            strand="."
        else:
            strand=self.strand
            
        for attribute in self.attributes:
            if attribute in ("Parent","ID"): continue
            value=self.attributes[attribute]
            #ttribute=attribute.lower()
            attribute=re.sub(";",":", attribute.lower())
            attr_field="{0};{1}={2}".format(attr_field,attribute, value)
        
        parent_line = [self.chrom, "locus_pipeline", "transcript", self.start, self.end, ".", strand, ".",  attr_field ]
        
        parent_line ="\t".join( str(s) for s in parent_line )
        
        exon_lines = []
        for index in range(len(self.exons)):
            exon=self.exons[index]
            exon_line = [self.chrom, "locus_pipeline", "exon", exon[0], exon[1],
                         ".", strand, ".",
                         "Parent={0};ID={0}.{1}".format(self.id, index) ]
            exon_lines.append("\t".join(str(s) for s in exon_line))
        
        lines=[parent_line]
        lines.extend(exon_lines) 
        return "\n".join(lines)
        
#     def __eq__(self, other):
#         if self.id==other.id and self.start==other.start and self.end==other.end and self.chrom==other.chrom: return True
#         else: return False
        
    def addExon(self, gffLine):
        '''This function will append an exon/CDS feature to the object.'''
        if gffLine.feature not in ("exon", "CDS"):
            raise AttributeError()
        
        if gffLine.attributes["Parent"]!=self.id:
            raise AssertionError("""Mismatch between transcript and exon:\n
            {0}\n
            {1}
            """.format(self.id, gffLine))
        start,end=sorted([gffLine.start, gffLine.end])
        self.exons.append((start, end) )
        
    def finalize(self):
        '''Function to calculate the internal introns from the exons.
        In the first step, it will sort the exons by their internal coordinates.'''
        
        # We do not want to repeat this step multiple times
        if self.finalized is True:
            return
        self.exons = sorted(self.exons, key=operator.itemgetter(0,1) ) # Sort the exons by start then stop
        if self.exons[0][0]!=self.start or self.exons[-1][1]!=self.end:
            raise ValueError("The transcript {id} has coordinates {tstart}:{tend}, but its first and last exons define it up until {estart}:{eend}!".format(
                                                                                                                                                            tstart=self.start,
                                                                                                                                                            tend=self.end,
                                                                                                                                                            id=self.id,
                                                                                                                                                            eend=self.exons[-1][1],
                                                                                                                                                            estart=self.exons[0][0],
                                                                                                                                                            ))
        
        if len(self.exons)==1:
            self.monoexonic = True
            return # There is no sense in performing any operation on single exon transcripts
        
        
        
        
        for index in range(len(self.exons)-1):
            exonA, exonB = self.exons[index:index+2]
            if exonA[1]>=exonB[0]:
                raise ValueError("Overlapping exons found!")
            self.junctions.append( (exonA[1]+1, exonB[0]-1) ) #Append the splice junction
            self.splices.extend( [exonA[1]+1, exonB[0]-1] ) # Append the splice locations
            
        self.junctions = set(self.junctions)
        self.splices = set(self.splices)
        self.finalized = True
        return