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
        self.exons, self.cds, self.utr = [], [], []
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
        
        cds_begin = False
        
        cds_count=0
        exon_count=0
        utr_count=0
        
        for segment in self.segments:
            if cds_begin is False and segment[0]=="CDS": cds_begin = True
            if segment[0]=="UTR":
                utr_count+=1
                index=utr_count
                if cds_begin is True:
                    if self.strand=="-": feature="five_prime_utr"
                    else: feature="three_prime_utr"
                else:
                    if self.strand=="-": feature="three_prime_utr"
                    else: feature="five_prime_utr"
            else:
                if segment[0]=="CDS":
                    cds_count+=1
                    index=cds_count
                else:
                    exon_count+=1
                    index=exon_count
                feature=segment[0]
            
            exon_line = [self.chrom, "locus_pipeline", feature, segment[1], segment[2],
                         ".", strand, ".",
                         "Parent={0};ID={0}.{1}".format(self.id, index) ]
            exon_lines.append("\t".join(str(s) for s in exon_line))
        
        
        lines=[parent_line]
        lines.extend(exon_lines) 
        return "\n".join(lines)
        
#     def __eq__(self, other):
#         if self.id==other.id and self.start==other.start and self.end==other.end and self.chrom==other.chrom: return True
#         else: return False
        
    @property
    def cds_length(self):
        return sum([ c[1]-c[0]+1 for c in self.cds ])
        
    @property
    def length(self):
        return sum([ e[1]-e[0]+1 for e in self.exons ])
    
    @property
    def utr_length(self):
        return sum([ e[1]-e[0]+1 for e in self.utr ])
        
    def addExon(self, gffLine):
        '''This function will append an exon/CDS feature to the object.'''
#         if gffLine.feature not in ("exon", "CDS") or "UTR" not in gffLine.feature:
#             raise AttributeError(gffLine.feature)
    
        if self.finalized is True:
            raise RuntimeError("You cannot add exons to a finalized transcript!")
        
        if gffLine.attributes["Parent"]!=self.id:
            raise AssertionError("""Mismatch between transcript and exon:\n
            {0}\n
            {1}
            """.format(self.id, gffLine))
        if gffLine.feature=="CDS":
            store=self.cds
        elif "UTR" in gffLine.feature:
            store=self.utr
        elif gffLine.feature=="exon":
            store=self.exons
            
        start,end=sorted([gffLine.start, gffLine.end])
        store.append((start, end) )
        
    def finalize(self):
        '''Function to calculate the internal introns from the exons.
        In the first step, it will sort the exons by their internal coordinates.'''
        
        # We do not want to repeat this step multiple times
        if self.finalized is True:
            return

        if self.utr!=[] and self.cds==[]:
            raise ValueError("Transcript {tid} has defined UTRs but no CDS feature!".format(tid=self.id))

        assert self.length == self.utr_length + self.cds_length, (self.length, self.utr_length, self.cds_length )

        self.exons = sorted(self.exons, key=operator.itemgetter(0,1) ) # Sort the exons by start then stop
        
        if self.exons[0][0]!=self.start or self.exons[-1][1]!=self.end:
            raise ValueError("The transcript {id} has coordinates {tstart}:{tend}, but its first and last exons define it up until {estart}:{eend}!".format(
                                                                                                                                                            tstart=self.start,
                                                                                                                                                            tend=self.end,
                                                                                                                                                            id=self.id,
                                                                                                                                                            eend=self.exons[-1][1],
                                                                                                                                                            estart=self.exons[0][0],
                                                                                                                                                            ))
        self.cds = sorted(self.cds, key=operator.itemgetter(0,1))
        self.utr = sorted(self.utr, key=operator.itemgetter(0,1))
        self.segments = [ ("exon",e[0],e[1]) for e in self.exons] + \
                    [("CDS", c[0],c[1]) for c in self.cds ] + \
                    [ ("UTR", u[0], u[1]) for u in self.utr ]
        self.segments = sorted(self.segments, key=operator.itemgetter(1,2) )
        
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