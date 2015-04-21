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

    def addExon(self, gffLine):
        '''This function will append an exon/CDS feature to the object.'''

        if self.finalized is True:
            raise RuntimeError("You cannot add exons to a finalized transcript!")
        
        if gffLine.parent!=self.id:
            raise AssertionError("""Mismatch between transcript and exon:\n
            {0}\n
            {1}
            """.format(self.id, gffLine))
        if gffLine.feature=="CDS":
            store=self.cds
        elif "utr" in gffLine.feature or "UTR" in gffLine.feature:
            store=self.utr
        elif gffLine.feature=="exon":
            store=self.exons
        else:
            raise AttributeError("Unknown feature: {0}".format(gffLine.feature))
            
        start,end=sorted([gffLine.start, gffLine.end])
        store.append((start, end) )

    def finalize(self):
        '''Function to calculate the internal introns from the exons.
        In the first step, it will sort the exons by their internal coordinates.'''
        
        # We do not want to repeat this step multiple times
        if self.finalized is True:
            return

        self.metrics=dict() # create the store for the metrics
        if self.utr!=[] and self.cds==[]:
            raise ValueError("Transcript {tid} has defined UTRs but no CDS feature!".format(tid=self.id))

        assert self.cds_length==self.utr_length==0 or  self.cdna_length == self.utr_length + self.cds_length, (self.id, self.cdna_length, self.utr_length, self.cds_length,
                                                                                                               self.utr, self.cds, self.exons )

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
            self.finalized = True
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


    ####################Class properties##################################

    @property
    def cds_length(self):
        '''This property return the length of the CDS part of the transcript.'''
        return sum([ c[1]-c[0]+1 for c in self.cds ])
        
    @property
    def cdna_length(self):
        '''This property returns the length of the transcript.'''
        return sum([ e[1]-e[0]+1 for e in self.exons ])
    
    @property
    def utr_length(self):
        '''This property return the length of the UTR part of the transcript.'''
        return sum([ e[1]-e[0]+1 for e in self.utr ])
    
    @property
    def internal_cds_num(self):
        '''This property returns the number of CDSs inside a transcript.'''
        
        return len(self.internal_cds)

    @property
    def internal_cds(self):
        '''This property calculates the CDSs inside a transcript.
        The property is a list of tuples, each of which
        describes a CDS segment and holds as tuples (start, end) the 
        CDS segments inside the single CDS.'''
        
        self.finalize()
        self.__internal_cds = []
        if len(self.utr)==0:
            self.__internal_cds = [ tuple(self.cds) ] 
        elif len(self.cds)>0:
            current_cds=[]
            in_utr=True
            #The sense of this cycle is to look for multiple CDSs. I exploit the fact that 
            #the transcript is a directed segment .. so by parsing linearly I can detect easily
            #instances where I have a UTR placed between two CDSs
            for segment in filter(lambda x: x[0]!="exon", self.segments):
                if segment[0]=="CDS":
                    if in_utr is True:
                        in_utr=False
                    current_cds.append(tuple([segment[1],segment[2]]))
                elif segment[0]=="UTR":
                    if in_utr is False:
                        if len(current_cds)>0:
                            self.__internal_cds.append(tuple(current_cds) )
                        current_cds=[]
                    in_utr=True
            if len(current_cds)>0:  
                self.__internal_cds.append(tuple(current_cds))
        assert sum(len(x) for x in self.__internal_cds) == len(self.cds)
        return self.__internal_cds
    
#     @internal_cds.setter
#     def internal_cds(self, *args):
#         if len(args)==0:
#             args.append(None)
#         assert len(args)==1 and type(args[0]) in (None,list)
#         self.__internal_cds = args[0]
    
    #internal_cds=property(get_internal_cds, set_internal_cds) 
        
    @property
    def max_internal_cds_length(self):
        '''This property calculates the length of the greatest CDS inside the cDNA.'''
        #_ = self.max_internal_cds #Calculate on the fly  
        if len(self.cds)==0:
            self.__max_internal_cds_length=0
        else:
            self.__max_internal_cds_length=sum(x[1]-x[0]+1 for x in self.max_internal_cds)
        return self.__max_internal_cds_length

#     @max_internal_cds_length.setter
#     def max_internal_cds_length(self, *args):
#         if len(*args)==0:
#             self.__max_internal_cds_index = None
#             args=[0]
#         assert len(args)==1 and type(args[0]) is int
#         self.__max_internal_cds_length = args[0]

    @property
    def max_internal_cds(self):
        '''This property will return the tuple of tuples of the CDS with the greatest length
        inside the transcript. To avoid memory wasting, the tuple is accessed in real-time using 
        a token (__max_internal_cds_index) which holds the position in the __internal_cds list of the longest CDS.'''
        if len(self.cds)==0: # Non-sense to calculate the maximum CDS for transcripts without it
            self.__max_internal_cds_length=0
            return tuple([])
        else:
            greatest=0
            self.__max_internal_cds_index=-1
            for index in range(len(self.__internal_cds)):
                cds=self.__internal_cds[index]
                length = sum( c[1]-c[0]+1  for c in cds  )
                if length>greatest:
                    #self.__max_internal_cds_length=length
                    self.__max_internal_cds_index=index
            if self.__max_internal_cds_index==-1: raise ValueError("""Index not modified for transcript {0}!
            Monoexonic: {1}; CDS length: {2};
            CDS: {3};
            greatest: {4};
            internal: {5};
            __internal: {6}""".format(self.id,
                               self.monoexonic,
                               self.cds_length,
                               self.cds,
                               greatest,
                               self.internal_cds,
                               self.__internal_cds))
        return self.__internal_cds[self.__max_internal_cds_index]

    @max_internal_cds.setter
    def max_internal_cds(self,*args):
        if len(args)==0:
            args.append(None)
        assert len(args)==1
        self.__max_internal_cds=args[0]
        