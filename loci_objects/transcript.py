import operator
from loci_objects.abstractlocus import abstractlocus # Needed for the BronKerbosch algorithm ...
from loci_objects.GTF import gtfLine
from loci_objects.GFF import gffLine


class transcript:
    
    '''This class defines a transcript, down to its exon/CDS/UTR components. It is instantiated by a transcript
    GT/GFF3 line.
    Key attributes:
    - chrom    The chromosome
    - source
    - feature            mRNA if at least one CDS is defined
    - start
    - end
    - score
    - strand            one of +,-,None
    - phase            set to None
    - id            the ID of the transcripts (or tid)
    - parent        
    - attributes    a dictionary with additional informations from the GFFline
    
    After all exons have been loaded into the instance (see "addExon"), the class must be finalized with the appropriate method.

    CDS locations can be uploaded from the external, using a dictionary of indexed BED12 entries. 
    '''
    
    
    ######### Class special methods ####################
    
    def __init__(self, gffLine, source=None):
        
        '''Initialise the transcript object, using a mRNA/transcript line.
        Note: I am assuming that the input line is an object from my own "GFF" class.
        The transcript instance must be initialised by a "(m|r|lnc|whatever)RNA" or "transcript" gffLine.'''
        
        self.chrom = gffLine.chrom
        assert "transcript"==gffLine.feature or "RNA" in gffLine.feature.upper()
        self.feature="transcript"
        self.id = gffLine.id
        self.name = gffLine.name
        if source is None:
            self.source=gffLine.source
        else:
            self.source=source
        self.start=gffLine.start
        self.strand = gffLine.strand
        self.end=gffLine.end
        self.score = gffLine.score
        self.exons, self.combined_cds, self.combined_utr = [], [], []
        self.junctions = []
        self.splices = []
        self.monoexonic = False
        self.finalized = False # Flag. We do not want to repeat the finalising more than once.
        self.parent = gffLine.parent
        self.attributes = gffLine.attributes
        self.best_internal_orf_index = None
        self.has_start, self.has_stop = False,False
        self.non_overlapping_cds = None
        
    def __str__(self, to_gtf=False):
        '''Each transcript will be printed out in the GFF style.
        This is pretty rudimentary, as the class does not hold any information on the original source, feature, score, etc.'''
        
        self.finalize() #Necessary to sort the exons
        lines = []
        transcript_counter = 0
#         assert self.best_internal_orf_index > -1
        
        for index in range(len(self.internal_cds)):
            
            if self.number_internal_orfs>1:
                if index==self.best_internal_orf_index: self.attributes["maximal"]=True
                else: self.attributes["maximal"]=False
            cds_run = self.internal_cds[index]
            
#             if len(list(filter(lambda x: x[0]=="UTR", cds_run)  )  )>0:
#                 assert len(list(filter(lambda x: x[0]=="CDS", cds_run)  )  )>0
            if self.number_internal_orfs>1:
                transcript_counter+=1
                tid = "{0}.orf{1}".format(self.id, transcript_counter)
            else:
                tid = self.id
            
            if self.strand is None:
                strand="."
            else:   
                strand=self.strand
            
#             for attribute in self.attributes:
#                 if attribute in ("Parent","ID"): continue
#                 value=self.attributes[attribute]
#                 #ttribute=attribute.lower()
#                 attribute=re.sub(";",":", attribute.lower())
#                 attr_field="{0};{1}={2}".format(attr_field,attribute, value)
            
            
            if self.score is None or self.score==float("-inf"):
                score="."
            else:
                score=round(self.score,2)
                
            if to_gtf is True:
                parent_line = gtfLine('')
            else:
                parent_line=gffLine('')
                
            parent_line.chrom=self.chrom
            parent_line.source=self.source
            parent_line.feature=self.feature
            parent_line.start,parent_line.end=self.start,self.end
            parent_line.score=score
            parent_line.strand=strand
            parent_line.phase='.'
            parent_line.attributes=self.attributes
            
            parent_line.parent=self.parent
            parent_line.id=tid
        
            exon_lines = []
        
            cds_begin = False
        
            cds_count=0
            exon_count=0
            five_utr_count=0
            three_utr_count=0

            for segment in cds_run:
                if cds_begin is False and segment[0]=="CDS": cds_begin = True
                if segment[0]=="UTR":
                    if cds_begin is True:
                        if to_gtf is True:
                            if self.strand=="-": feature="5UTR"
                            else: feature="3UTR"
                        else:
                            if self.strand=="-": feature="five_prime_UTR"
                            else: feature="three_prime_UTR"
                    else:
                        if to_gtf is True:
                            if self.strand=="-": feature="3UTR"
                            else: feature="5UTR"
                        else:
                            if self.strand=="-": feature="three_prime_UTR"
                            else: feature="five_prime_UTR"
                    if "five" in feature or "5" in feature:
                        five_utr_count+=1
                        index=five_utr_count
                    else:
                        three_utr_count+=1
                        index=three_utr_count
                else:
                    if segment[0]=="CDS":
                        cds_count+=1
                        index=cds_count
                    else:
                        exon_count+=1
                        index=exon_count
                    feature=segment[0]
                if to_gtf is True:
                    exon_line=gtfLine('')
                else:
                    exon_line=gffLine('')
                    
                exon_line.chrom=self.chrom
                exon_line.source=self.source
                exon_line.feature=feature
                exon_line.start,exon_line.end=segment[1],segment[2]
                exon_line.strand=strand
                exon_line.phase=None
                exon_line.score = None
                if to_gtf is True:
                    exon_line.gene=self.parent
                    exon_line.transcript=tid
                else:
                    exon_line.id="{0}.{1}{2}".format(tid, feature,index)
                    exon_line.parent=tid
                    
                exon_lines.append(str(exon_line))
        
        
            lines.append(str(parent_line))
            lines.extend(exon_lines) 
        
        return "\n".join(lines)
    
    def __eq__(self, other):
        '''Two transcripts are considered identical if they have the same
        start, end, chromosome, strand and internal exons.
        IDs are not important for this comparison; two transcripts coming from different
        methods and having different IDs can still be identical.'''
        
        if not type(self)==type(other): return False
        self.finalize()
        other.finalize()
           
        if self.strand == other.strand and self.chrom == other.chrom and \
            self.start==other.start and self.end == other.end and \
            self.exons == other.exons:
            return True
          
        return False
    
    def __hash__(self):
        '''Returns the hash of the object (call to super().__hash__()).
        Necessary to be able to add these objects to hashes like sets.'''

        return super().__hash__()
    
    def __len__(self):
        '''Returns the length occupied by the unspliced transcript on the genome.'''
        return self.end-self.start+1

     
    def __lt__(self, other):
        '''A transcript is lesser than another if it is on a lexicographic inferior chromosome,
        or if it begins before the other, or (in the case where they begin at the same location)
        it ends earlier than the other.'''
        if self.chrom!=other.chrom:
            return self.chrom<other.chrom
        if self==other:
            return False
        if self.start<other.start:
            return True
        elif self.start==other.start and self.end<other.end:
            return True
        return False
     
    def __gt__(self, other):
        return not self<other
     
    def __le__(self, other):
        return (self==other) or (self<other)
     
    def __ge__(self, other):
        return (self==other) or (self>other)          
    
    ######### Class instance methods ####################


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
            store=self.combined_cds
        elif "combined_utr" in gffLine.feature or "UTR" in gffLine.feature:
            store=self.combined_utr
        elif gffLine.feature=="exon":
            store=self.exons
        elif gffLine.feature=="start_codon":
            self.has_start = True
            return
        elif gffLine.feature=="stop_codon":
            self.has_stop = True
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

        if len(self.exons)>1 and self.strand is None:
            raise AttributeError("Multiexonic transcripts must have a defined strand! Error for {0}".format(self.id))

        if self.combined_utr!=[] and self.combined_cds==[]:
            raise ValueError("Transcript {tid} has defined UTRs but no CDS feature!".format(tid=self.id))

        assert self.combined_cds_length==self.combined_utr_length==0 or  self.cdna_length == self.combined_utr_length + self.combined_cds_length, (self.id, self.cdna_length, self.combined_utr_length, self.combined_cds_length,
                                                                                                               self.combined_utr, self.combined_cds, self.exons )

        self.exons = sorted(self.exons, key=operator.itemgetter(0,1) ) # Sort the exons by start then stop
#         assert len(self.exons)>0
        try:
            if self.exons[0][0]<self.start or self.exons[-1][1]>self.end:
                raise ValueError("The transcript {id} has coordinates {tstart}:{tend}, but its first and last exons define it up until {estart}:{eend}!".format(
                                                                                                                                                            tstart=self.start,
                                                                                                                                                            tend=self.end,
                                                                                                                                                            id=self.id,
                                                                                                                                                            eend=self.exons[-1][1],
                                                                                                                                                            estart=self.exons[0][0],
                                                                                                                                                            ))
        except IndexError as err:
            raise IndexError(err, self.id, str(self.exons))
            
        
        self.combined_cds = sorted(self.combined_cds, key=operator.itemgetter(0,1))
        self.combined_utr = sorted(self.combined_utr, key=operator.itemgetter(0,1))
        if len(self.combined_utr)>0 and self.combined_utr[0][0]<self.combined_cds[0][0]:
            if self.strand=="+": self.has_start=True
            elif self.strand=="-": self.has_stop=True
        if len(self.combined_utr)>0 and self.combined_utr[-1][1]>self.combined_cds[-1][1]:
            if self.strand=="+": self.has_stop=True
            elif self.strand=="-": self.has_start=True
        
        self.internal_cds, self.junctions, self.splices = [], [], []
        self.segments = [ ("exon",e[0],e[1]) for e in self.exons] + \
                    [("CDS", c[0],c[1]) for c in self.combined_cds ] + \
                    [ ("UTR", u[0], u[1]) for u in self.combined_utr ]
        self.segments =  sorted(self.segments, key=operator.itemgetter(1,2,0) )
                
        self.internal_cds.append(self.segments)
        
        if len(self.exons)==1:
            #self.finalized = True
            self.monoexonic = True
            #return # There is no sense in performing any operation on single exon transcripts
        else:
            for index in range(len(self.exons)-1):
                exonA, exonB = self.exons[index:index+2]
                if exonA[1]>=exonB[0]:
                    raise ValueError("Overlapping exons found!")
                self.junctions.append( (exonA[1]+1, exonB[0]-1) ) #Append the splice junction
                self.splices.extend( [exonA[1]+1, exonB[0]-1] ) # Append the splice locations
            
        self.junctions = set(self.junctions)
        self.splices = set(self.splices)
        _ = self.best_internal_orf
#         assert self.best_internal_orf_index > -1
        if len(self.combined_cds)>0:
            self.feature="mRNA"
        
        self.finalized = True
        return


    def load_cds(self, cds_dict, trust_strand=False):
        
        '''Arguments:
        - cds_dict        a dictionary (indexed on the TIDs) that holds BED12 information on the transcript
        - trust_strand    for monoexonic transcripts, whether to trust the strand information or to overwrite it with the one provided by TD.
        
        
        This function is used to load the various CDSs from an external dictionary, loaded from a BED file.
        It replicates what is done internally by the "cdna_alignment_orf_to_genome_orf.pl" utility in the
        TransDecoder suite.
        The method expects as argument a dictionary containing BED entries, and indexed by the transcript name.
        The indexed name *must* equal the "id" property, otherwise the method returns immediately. 
        If no entry is found for the transcript, the method exits immediately. Otherwise, any CDS information present in
        the original GFF/GTF file is completely overwritten.
        Briefly, it follows this logic:
        - Finalise the transcript
        - Retrieve from the dictionary (input) the CDS object
        - Sort in decreasing order the CDSs on the basis of:
            - Presence of start/stop codon
            - CDS length (useful for monoexonic transcripts where we might want to set the strand)
        - For each CDS:
            - If the ORF is on the + strand:
                - all good
            - If the ORF is on the - strand:
                - if the transcript is monoexonic: invert its genomic strand
                - if the transcript is multiexonic: skip
            - Start looking at the exons
        '''
        self.finalize()

        if self.id not in cds_dict:
            return        

        self.combined_utr = []
        self.combined_cds = []
        self.internal_cds = []
        
        self.finalized = False
        
        #Ordering the CDSs by: presence of start/stop codons, combined_cds length
        original_strand = self.strand
        new_strand = None
        
        best_cds=True # Token to be set to False after the first CDS is exhausted 
        for cds_run in sorted(cds_dict[self.id], reverse=True, key=operator.attrgetter("cds_len") ):
            
            cds_start, cds_end, strand = cds_run.cdsStart, cds_run.cdsEnd, cds_run.strand
            if not (cds_start>=1 and cds_end<=self.cdna_length):
                continue
            
            if best_cds:
                self.has_start, self.has_stop = cds_run.has_start, cds_run.has_stop
            best_cds=False # Exhaust the token

            # assert cds_start>=1 and cds_end<=self.cdna_length, ( self.id, self.cdna_length, (cds_start,cds_end) )
            
            if self.strand is None:
                self.strand=new_strand=strand
            elif strand == "-":
                #Basic case
                if self.monoexonic is False:
                    continue
                if trust_strand is True:
                    continue
                #Case 1 we trust the strand
                if new_strand is not None:
                    if original_strand is None:
                    #We already assigned to + strand
                        if self.strand=="+":
                            continue
                        else:
                            pass
                    elif original_strand=="+":
                        if self.strand=="+":
                            continue
                        else:
                            pass
                    elif original_strand=="-":
                        if self.strand=="+":
                            pass
                        else:
                            continue
                else:
                    if self.strand=="+":
                        self.strand=new_strand="-"
                    elif self.strand=="-":
                        self.strand=new_strand="+"
                
            cds_exons = []
            current_start, current_end = 0,0
            if self.strand == "+":
                for exon in sorted(self.exons, key=operator.itemgetter(0,1)):
                    cds_exons.append(("exon", exon[0], exon[1] ) )
                    current_start+=1
                    current_end+=exon[1]-exon[0]+1
                    #Whole UTR
                    if current_end<cds_start or current_start>cds_end:
                        cds_exons.append( ("UTR", exon[0], exon[1])  )
                    else:
                        c_start = exon[0] + max(0, cds_start-current_start )
                        if c_start > exon[0]:
                            u_end = c_start-1
                            cds_exons.append( ("UTR", exon[0], u_end) )
                        c_end = exon[1] - max(0, current_end - cds_end )
#                         assert c_end>=exon[0] 
                        if c_start<c_end:
                            cds_exons.append(("CDS", c_start, c_end))
                        if c_end < exon[1]:
                            cds_exons.append( ("UTR", c_end+1, exon[1]  ) )
                    current_start=current_end
                            
            elif self.strand=="-":
                for exon in sorted(self.exons, key=operator.itemgetter(0,1), reverse=True):
                    cds_exons.append(("exon", exon[0], exon[1] ) )
                    current_start+=1
                    current_end+=exon[1]-exon[0]+1
                    if current_end<cds_start or current_start>cds_end:
                        cds_exons.append( ("UTR", exon[0], exon[1] ))
                    else:
                        c_end = exon[1] - max(0,cds_start - current_start ) 
#                         assert c_end>=exon[0]
                        if c_end < exon[1]:
                            cds_exons.append(("UTR", c_end+1, exon[1]))
                        c_start = exon[0] + max(0, current_end - cds_end )
                        cds_exons.append( ("CDS", c_start, c_end) )
                        if c_start>exon[0]:
                            cds_exons.append( ("UTR", exon[0], c_start-1) )
                    current_start=current_end
        
            self.internal_cds.append( sorted(cds_exons, key=operator.itemgetter(1,2)   ) )

        if len(self.internal_cds)==1:
            self.combined_cds = sorted(
                              [(a[1],a[2]) for a in filter(lambda x: x[0]=="CDS", self.internal_cds[0])],
                              key=operator.itemgetter(0,1)
                              
                              )
            self.combined_utr = sorted(
                              [(a[1],a[2]) for a in filter(lambda x: x[0]=="UTR", self.internal_cds[0])],
                              key=operator.itemgetter(0,1)
                              
                              )
            
            
        elif len(self.internal_cds)>1:
            
            cds_spans = []
            candidates = []
            for internal_cds in self.internal_cds:
                candidates.extend([tuple([a[1],a[2]]) for a in filter(lambda tup: tup[0]=="CDS", internal_cds  )])
                              
            candidates=set(candidates)
            original=candidates.copy()
            for mc in self.merge_cliques(list(self.BronKerbosch(set(), candidates, set(), original,
                                                           ))):
                span=tuple([min(t[0] for t in mc),
                            max(t[1] for t in mc)                        
                            ])
                cds_spans.append(span)
                

            self.combined_cds = sorted(cds_spans, key = operator.itemgetter(0,1))
#             cds_spans=sorted(cds_spans, key = operator.itemgetter(0,1))
#             new_cds_spans = []
#             for span in cds_spans:
#                 if new_cds_spans == []:
#                     new_cds_spans.append(span)
#                 elif new_cds_spans[-1][1]<span[0]:
#                     new_cds_spans.append(span)
#                 elif new_cds_spans[-1]==span[0]:
#                     new_cds_spans[-1]=(new_cds_spans[-1][0], span[1])
#             self.combined_cds = new_cds_spans
            
            #This method is probably OBSCENELY inefficient, but I cannot think of a better one for the moment.
            curr_utr_segment = None

            utr_pos = set.difference( 
                                                  set.union(*[ set(range(exon[0],exon[1]+1)) for exon in self.exons]),
                                                  set.union(*[ set(range(cds[0],cds[1]+1)) for cds in self.combined_cds])
                                                  )
            for pos in sorted(list(utr_pos)):
                if curr_utr_segment is None:
                    curr_utr_segment = (pos,pos)
                else:
                    if pos==curr_utr_segment[1]+1:
                        curr_utr_segment = (curr_utr_segment[0],pos)
                    else:
                        self.combined_utr.append(curr_utr_segment)
                        curr_utr_segment = (pos,pos)
                        
            if curr_utr_segment is not None:
                self.combined_utr.append(curr_utr_segment)           
                                   
            assert self.cdna_length == self.combined_cds_length + self.combined_utr_length, (self.cdna_length, self.combined_cds, self.combined_utr)                            
        
        if self.internal_cds == []:
            self.finalize()
        else:
            self.feature="mRNA"
            self.finalized=True
        return
                        

      
    @classmethod
    ####################Class methods#####################################  
    def BronKerbosch(cls, clique, candidates, non_clique, original):
        '''Wrapper for the abstractlocus method. It will pass to the function the class's "is_intersecting" method
        (which would be otherwise be inaccessible from the abstractlocus class method)'''
        for clique in abstractlocus.BronKerbosch(clique, candidates, non_clique, original,
                                                 inters=cls.is_intersecting,
                                                 neighbours = None   ):
            yield clique
    
    @classmethod
    def is_intersecting(cls, first, second):
        '''Implementation of the is_intersecting method.'''
        if first==second or cls.overlap(first,second)<0:
            return False
        return True

    @classmethod
    def overlap(cls, first,second):
        lend = max(first[0], second[0])
        rend = min(first[1], second[1])
        return rend-lend
    
    @classmethod
    def merge_cliques(cls, cliques):
        '''Wrapper for the abstractlocus method.'''
        return abstractlocus.merge_cliques(cliques)

    ####################Class properties##################################

    @property
    def id(self):
        '''ID of the transcript - cannot be an undefined value.'''
        return self.__id
    
    @id.setter
    def id(self, Id):
        if type(Id) is not str:
            raise ValueError("Invalid value for id: {0}, type {1}".format(
                                                                          Id, type(Id)))
        self.__id = Id

    @property
    def tid(self):
        '''ID of the transcript - cannot be an undefined value. Alias of id.'''
        return self.id
    
    @tid.setter
    def tid(self,tid):
        self.id=tid
        
    @property
    def parent(self):
        '''Name of the parent feature of the transcript.'''
        return self.__parent
    
    @parent.setter
    def parent(self,parent):
        if parent is not None and type(parent) is not str:
            raise ValueError("Invalid value for parent: {0}, type {1}".format(
                                                                          parent, type(parent)))
            
        self.__parent=parent
        
    @property
    def score(self):
        '''Numerical value which summarizes the reliability of the transcript.'''
        return self.__score
        
    @score.setter
    def score(self,score):
        if score is not None and type(score) not in (float,int):
            try: score=float(score)
            except:
                raise ValueError("Invalid value for score: {0}, type {1}".format(
                                                                          score, type(score)))
        self.__score=score
                
        

    @property
    def combined_cds_length(self):
        '''This property return the length of the CDS part of the transcript.'''
        return sum([ c[1]-c[0]+1 for c in self.combined_cds ])
    
    @property
    def combined_cds_num(self):
        '''This property returns the number of non-overlapping CDS segments in the transcript.'''
        return len( self.combined_cds )

    @property
    def combined_cds_num_fraction(self):
        '''This property returns the fraction of non-overlapping CDS segments in the transcript
        vs. the total number of exons'''
        return len( self.combined_cds )/len(self.exons)

    @property
    def combined_cds_fraction(self):
        '''This property return the length of the CDS part of the transcript.'''
        return self.combined_cds_length/self.cdna_length
    
    @property
    def combined_utr_length(self):
        '''This property return the length of the UTR part of the transcript.'''
        return sum([ e[1]-e[0]+1 for e in self.combined_utr ])
    
    @property
    def combined_utr_fraction(self):
        '''This property returns the fraction of the cDNA which is not coding according
        to any ORF. Complement of combined_cds_fraction'''
        return 1-self.combined_cds_fraction
        
    @property
    def cdna_length(self):
        '''This property returns the length of the transcript.'''
        return sum([ e[1]-e[0]+1 for e in self.exons ])
    
    @property
    def number_internal_orfs(self):
        '''This property returns the number of CDSs inside a transcript.'''
        return len(self.internal_cds)

    @property
    def cds_length(self):
        '''This property calculates the length of the greatest CDS inside the cDNA.'''
        if len(self.combined_cds)==0:
            self.__max_internal_orf_length=0
        else:
            self.__max_internal_orf_length=sum(x[2]-x[1]+1 for x in filter(lambda x: x[0]=="CDS", self.best_internal_orf))
        return self.__max_internal_orf_length
    
    @property
    def cds_num(self):
        '''This property calculates the number of CDS exons for the greatest ORF'''
        return len(list( filter(lambda exon: exon[0]=="CDS", self.best_internal_orf) ))
    
    @property
    def cds_fraction(self):
        '''This property calculates the fraction of the greatest CDS vs. the cDNA length.'''
        return self.__max_internal_orf_length/self.cdna_length
    
    @property
    def best_cds_exons_num(self):
        '''Returns the number of CDS segments in the greatest ORF (irrespective of the number of exons involved)'''
        return len(list(filter(lambda x: x[0]=="CDS", self.best_internal_orf)))
    
    @property
    def best_cds_exons_fraction(self):
        '''Returns the fraction of CDS segments in the greatest ORF (irrespective of the number of exons involved)'''
        return len(list(filter(lambda x: x[0]=="CDS", self.best_internal_orf)))/len(self.exons)
    

    @property
    def best_cds_number(self):
        '''This property returns the maximum number of CDS segments among the ORFs; this number
        can refer to an ORF *DIFFERENT* from the maximal ORF.'''
        cds_numbers = []
        for cds in self.internal_cds:
            cds_numbers.append(len(list(filter(lambda x: x[0]=="CDS", cds))))
        return max(cds_numbers)
    
    @property
    def best_cds_number_fraction(self):
        '''This property returns the proportion of best possible CDS segments vs. the number of exons.
        See best_cds_number.'''
        return self.best_cds_number/self.exon_num

    @property
    def best_internal_orf(self):
        '''This property will return the tuple of tuples of the CDS with the greatest length
        inside the transcript. To avoid memory wasting, the tuple is accessed in real-time using 
        a token (__max_internal_orf_index) which holds the position in the __internal_cds list of the longest CDS.'''
        if len(self.combined_cds)==0: # Non-sense to calculate the maximum CDS for transcripts without it
            self.__max_internal_orf_length=0
            self.best_internal_orf_index=0
            return tuple([])
        else:
            return self.internal_cds[self.best_internal_orf_index]
    
    @property
    def cds_not_maximal(self):
        '''This property returns the length of the CDS excluded from the longest CDS.'''
        if len(self.internal_cds)<2:
            return 0
        return self.combined_cds_length-self.cds_length
    
    @property
    def cds_not_maximal_fraction(self):
        '''This property returns the fraction of bases not in the greatest ORF compared to
        the total number of CDS bases in the cDNA.'''
        if self.combined_cds_length==0:
            return 0
        else:
            return self.cds_not_maximal/self.combined_cds_length
    
    @property
    def five_utr_length(self):
        '''Returns the length of the 5' UTR of the greatest ORF.'''
        if len(self.combined_cds)==0:
            return 0
        return sum(x[2]-x[1]+1 for x in self.five_utr)
                            
    @property
    def five_utr(self):
        '''Returns the exons in the 5' UTR of the greatest ORF.'''
        if len(self.combined_cds)==0:
            return []
        if self.strand=="+":
            cds_start=min(x[1] for x in filter(lambda exon: exon[0]=="CDS", self.best_internal_orf)  )
            return list(filter( lambda exon: exon[0]=="UTR" and exon[2]<cds_start, self.best_internal_orf  )  )
        elif self.strand=="-":
            cds_start=max(x[2] for x in filter(lambda exon: exon[0]=="CDS", self.best_internal_orf)  )
            return list(filter( lambda exon: exon[0]=="UTR" and exon[1]>cds_start, self.best_internal_orf  )  )

    @property
    def five_utr_num(self):
        '''This property returns the number of 5' UTR segments for the greatest ORF.'''
        return len(self.five_utr)

    @property
    def three_utr_length(self):
        '''Returns the length of the 5' UTR of the greatest ORF.'''
        if len(self.combined_cds)==0:
            return 0
        return sum(x[2]-x[1]+1 for x in self.five_utr)
                            
    @property
    def three_utr(self):
        '''Returns the exons in the 3' UTR of the greatest ORF.'''
        if len(self.combined_cds)==0:
            return []
        if self.strand=="-":
            cds_end=min(x[1] for x in filter(lambda exon: exon[0]=="CDS", self.best_internal_orf)  )
            return list(filter( lambda exon: exon[0]=="UTR" and exon[2]<cds_end, self.best_internal_orf  )  )
        elif self.strand=="+":
            cds_end=max(x[2] for x in filter(lambda exon: exon[0]=="CDS", self.best_internal_orf)  )
            return list(filter( lambda exon: exon[0]=="UTR" and exon[1]>cds_end, self.best_internal_orf  )  )

    @property
    def three_utr_num(self):
        '''This property returns the number of 3' UTR segments (referred to the greatest ORF).'''
        return len(self.three_utr)

    @property
    def utr_num(self):
        '''Returns the number of UTR segments (referred to the greatest ORF).'''
        return len(self.three_utr+self.five_utr)

    @property
    def utr_fraction(self):
        '''This property calculates the length of the UTR of the greatest CDS vs. the cDNA length.'''
        return 1-self.cds_fraction

    @property
    def utr_length(self):
        '''Returns the sum of the 5'+3' UTR lengths'''
        return self.three_utr_length+self.five_utr_length
    
    @property
    def has_start(self):
        '''Boolean. True if the greatest ORF has a start codon.'''
        return self.__has_start
    
    @has_start.setter
    def has_start(self, *args):
        if args[0] not in (None, False,True):
            raise TypeError("Invalid value for has_start: {0}".format(type(args[0])))
        self.__has_start=args[0]
        
    @property
    def has_stop(self):
        '''Boolean. True if the greatest ORF has a stop codon.'''
        return self.__has_stop
    
    @has_stop.setter
    def has_stop(self, *args):
        if args[0] not in (None, False,True):
            raise TypeError("Invalid value for has_stop: {0}".format(type(args[0])))
        self.__has_stop=args[0]

    @property
    def is_complete(self):
        '''Boolean. True if the best ORF has both start and end.'''
        return self.has_start and self.has_stop

    @property
    def best_internal_orf_index(self):
        '''Token which memorises the position in the ORF list of the best ORF.'''
        return self.__max_internal_orf_index
    
    @best_internal_orf_index.setter
    def best_internal_orf_index(self, index):
        if index is None:
            self.__max_internal_orf_index = index
            return
        if type(index) is not int:
            raise TypeError()
        if index<0 or index>=len(self.internal_cds):
            raise IndexError("No ORF corresponding to this index: {0}".format(index))
        self.__max_internal_orf_index = index
        
    @property
    def internal_orf_lengths(self):
        '''This property returns a list of the lengths of the internal ORFs.'''
        lengths = []
        for internal_cds in self.internal_cds:
            lengths.append( sum( x[2]-x[1]+1 for x in filter(lambda c: c[0]=="CDS", internal_cds) ) )
        lengths = sorted(lengths, reverse=True)
        return lengths
        
    @property
    def non_overlapping_cds(self):
        '''This property returns a set containing the set union of all CDS segments inside the internal CDSs.
        In the case of a transcript with no CDS, this is empty.
        In the case where there is only one CDS, this returns the combined_cds holder.
        In the case instead where there are multiple CDSs, the property will calculate the set union of all CDS segments.'''
        if self.__non_overlapping_cds is None: 
            self.finalize()
            self.__non_overlapping_cds=set()
            for internal_cds in self.internal_cds:
                segments = set([(x[1],x[2]) for x in filter(lambda segment: segment[0]=="CDS", internal_cds   )])
                self.__non_overlapping_cds.update(segments)
        return self.__non_overlapping_cds
    
    @non_overlapping_cds.setter
    def non_overlapping_cds(self,arg):
        '''Setter for the non_overlapping_cds property.'''
        self.__non_overlapping_cds = arg
    
    @property
    def exons(self):
        '''This property stores the exons of the transcript as (start,end) tuples.'''
        return self.__exons
    
    @exons.setter
    def exons(self, *args):
        if type(args[0]) not in (set,list):
            raise TypeError(type(args[0]))
        self.__exons=args[0]
    
    @property
    def exon_num(self):
        '''This property returns the number of exons of the transcript.'''
        return len(self.exons)
    
    @property
    def exon_fraction(self):
        '''This property returns the fraction of exons of the transcript which are contained in the sublocus.
        If the transcript is by itself, it returns None.'''
        
        return self.__exon_fraction
    
    @exon_fraction.setter
    def exon_fraction(self, *args):
        if type(args[0]) is not float or (args[0]<=0 or args[0]>1):
            raise TypeError("Invalid value for the fraction: {0}".format(args[0]))
        self.__exon_fraction=args[0]
    
    @property
    def intron_fraction(self):
        '''This property returns the fraction of introns of the transcript vs. the total number of introns in the locus.
        If the transcript is by itself, it returns 1.'''
        return self.__intron_fraction
    
    @intron_fraction.setter
    def intron_fraction(self, *args):
        if type(args[0]) not in (float,int) or (args[0]<0 or args[0]>1):
            raise TypeError("Invalid value for the fraction: {0}".format(args[0]))
        self.__intron_fraction=args[0]
    
    @property
    def retained_intron_num(self):
        '''This property records the number of introns in the transcripts which are marked as being retained.
        See the corresponding method in the sublocus class.'''
        return len(self.retained_introns)
    
    @property
    def retained_fraction(self):
        '''This property returns the fraction of the cDNA which is contained in retained introns.'''
        return self.__retained_fraction        

    @retained_fraction.setter
    def retained_fraction(self, *args):
        if type(args[0]) not in (float,int) or (args[0]<0 or args[0]>1):
            raise TypeError("Invalid value for the fraction: {0}".format(args[0]))
        self.__retained_fraction=args[0]
        
