#!/usr/bin/env python3

import operator # standard python library
import re

class transcript:
    
    def __init__(self, gffLine):
        
        '''Initialise the transcript object, using a mRNA/transcript line.
        Note: I am assuming that the input line is an object from my own "GFF" class.
        The transcript instance must be initialised by a "(m|r|lnc|whatever)RNA" or "transcript" gffLine.'''
        
#         if gffLine.feature!="transcript" or "RNA" not in gffLine.feature.upper():
#             raise AttributeError("Wrong feature line provided:\n{0}\n{1}".format( gffLine.feature, gffLine ))
        
        self.chrom = gffLine.chrom
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
            attribute=attribute.lower()
            attribute=re.sub(";",":", attribute)
            value=re.sub(";",":",attribute)
            attr_field="{0};{1}={2}".format(attribute, value)
        
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
        
        if len(self.exons)==1:
            self.monoexonic = True
            return # There is no sense in performing any operation on single exon transcripts
        
        self.exons = sorted(self.exons, key=operator.itemgetter(0,1) ) # Sort the exons by start then stop
        
        
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
            
class superlocus:
    
    '''The superlocus class is used to define overlapping regions on the genome, and it receives as input
    transcript class instances.'''
    
    def __init__(self, transcript):
        
        '''The superlocus class is instantiated from a transcript class, which it copies in its entirety.
        
        It will therefore have the following attributes:
        - chrom, strand, start, end
        - splices - a *set* which contains the position of each splice site
        - junctions - a *set* which contains the positions of each *splice junction* (registered as 2-tuples)
        - transcripts - a *set* which holds the transcripts added to the superlocus'''
        
        transcript.finalize()
        self.__dict__.update(transcript.__dict__)
        self.splices = set(self.splices)
        self.junctions = set(self.junctions)
        self.transcripts = set()
        self.transcripts.add(transcript)
        return
    
    @staticmethod
    def overlap(a, b, flank=0):
        '''This static method returns the overlap between two intervals. Values<=0 indicate no overlap.
        The optional "flank" argument (default 0) allows to expand a superlocus upstream and downstream.
        As a static method, it can be used also outside of any instance - "superlocus.overlap()" will function.
        The method takes as input either two objects/namedtuples that have the "start", "end" attributes,
        or two 2-tuples. 
        '''
        if hasattr(a, "start") and hasattr(b,"start"):
            right_boundary=min(a.end+flank, b.end+flank)
            left_boundary=max(a.start-flank, b.start-flank)
        elif type(a)==type(b)==tuple and len(a)==len(b)==2:
            right_boundary=min(a[1]+flank, b[1]+flank)
            left_boundary=max(a[0]-flank, b[0]-flank)
        
        return right_boundary + 1 - left_boundary #+1 because otherwise the intervals (15,20),(20,25) have an overlap of 0
    
    @classmethod
    def in_superlocus(cls, superlocus, transcript, flank=0):
        '''Function to determine whether a transcript should be added or not to the superlocus.
        This is a class method, i.e. it can be used also unbound from any specific instance of the class.
        It will be possible therefore to use it to compare any superlocus to any transcript.
        Arguments: 
        - a "superlocus" object
        - a "transcript" object (it must possess the "finalize" method)
        - flank - optional keyword'''
        transcript.finalize()
        if superlocus.chrom == transcript.chrom and \
            superlocus.strand == transcript.strand and \
            cls.overlap( superlocus, transcript, flank=flank  ) > 0:
            return True
        return False 
        
    def add_to_superlocus(self, transcript):
        '''This method checks that a transcript is contained within the superlocus (using the "in_superlocus" class method) and
        upon a successful check extends the superlocus with the new transcript.
        More precisely, it updates the boundaries (start and end), adds the transcript to the internal "transcripts" store,
        and extends the splices and junctions with those found inside the transcript.'''
        transcript.finalize()
        if self.in_superlocus(self, transcript) is True:
            self.start = min(self.start, transcript.start)
            self.end = max(self.end, transcript.end)
            self.transcripts.add(transcript)
            self.splices=set.union(self.splices, transcript.splices)
            self.junctions=set.union(self.junctions, transcript.splices)
        return
    
    @classmethod
    def is_intersecting(cls,transcript, other):
        '''I use this class method to verify that two transcripts are intersecting.
        If both are multiexonic, this amounts to checking whether there is at least one intron in common.
        If both are monoexonic, this amounts to checking whether there is some overlap between them.
        If one is monoexonic and the other is not, this function will return False by definition.        
        '''
        
        if transcript.id==other.id: return False # We do not want intersection with oneself
        monoexonic_check = len( list(filter(lambda x: x.monoexonic is True, [transcript, other]   )  )   )
        
        if monoexonic_check==0: #Both multiexonic
            inter=set.intersection(transcript.junctions, other.junctions )
            #print(transcript.id, other.id, inter)
            if len(inter)>0: return True
        
        elif monoexonic_check==1: #One monoexonic, the other multiexonic: different subloci by definition
            return False
        
        elif monoexonic_check==2:
            if cls.overlap(transcript, other)>0: #A simple overlap analysis will suffice
                return True
        
        return False
    
    @classmethod    
    def BronKerbosch(cls, clique, candidates, non_clique ):
        '''Implementation of the Bron-Kerbosch algorithm with pivot to define the subloci.
        We are using the class method "is_intersecting" to define the neighbours.
        Wiki: http://en.wikipedia.org/wiki/Bron%E2%80%93Kerbosch_algorithm '''
        
        if len(candidates)==0 and len(non_clique)==0:
            
            return clique
       
        #Define the pivot
        pivot = list(set.union(candidates, non_clique ))[0] #Choose a random member of the remaining non-neighbours
        pivot_neighbours = set(filter(
                                lambda candidate: cls.is_intersecting(pivot, candidate),
                                candidates )
                         )
        #Choose the vertices not neighbours of the pivot
        
        for vertex in set.difference(candidates, pivot_neighbours):
            vertex_neighbours = set(filter(lambda candidate: cls.is_intersecting(vertex, candidate), candidates ))
            #Now the recursive part:
            #Arguments:
            # clique u {vertex}
            # candidates intersect {vertex_neighbours}
            # non-clique intersect {vertex_neighbours}
            
            result = cls.BronKerbosch(
                                      set.union(clique, set([vertex])),
                                      set.intersection(candidates, vertex_neighbours),
                                      set.intersection(non_clique, vertex_neighbours)
                                      )
            if result is not None:
                return result
            
            candidates.remove(vertex)
            non_clique.add(vertex)
    
    def define_subloci(self):
        '''This method will define all subloci inside the superlocus.
        The method performs multiple calls to the BronKerbosch class method to define the possible groups of transcripts.        
        '''
        candidates = set(self.transcripts) # This will order the transcripts based on their position
        if len(candidates)==0:
            raise ValueError("This superlocus has no transcripts in it!")
        
        subloci = dict()
        index=0
        
        while len(candidates)>0:
            result=self.BronKerbosch(set(), candidates, set())
            index+=1
            subloci[index]=result
            candidates = set.difference( candidates, result )

        #Now we should define each sublocus and store it in a permanent structure of the class
        self.subloci = []
    
        for sublocus in subloci:
            transcripts = list(subloci[sublocus])
            if len(transcripts)==0:
                continue
            new_superlocus = superlocus(transcripts[0])
            if len(transcripts)>1:
                for ttt in transcripts[1:]:
                    new_superlocus.add_to_superlocus(ttt)
                    
            self.subloci.append((new_superlocus, transcripts[0].monoexonic))
    
    def __str__(self):
        
        '''Before printing, the class calls the define_subloci method. It will then print:
        # a "superlocus" line
        # for each "sublocus":
        ## a "sublocus" line
        ## all the transcripts inside the sublocus (see the transcript class)'''

        if self.strand is not None:
            strand=self.strand
        else:
            strand="."
        
        self.define_subloci()
        superlocus_id = "superlocus:{0}{3}:{1}-{2}".format(self.chrom, self.start, self.end,strand)

        superlocus_line = [self.chrom, "locus_pipeline", "superlocus", self.start, self.end, ".", strand, ".", "ID={0}".format(superlocus_id) ]
        superlocus_line = "\t".join(str(s) for s in superlocus_line)
        
        sublocus_lines = []
        counter=0
        for sublocus in self.subloci:
            counter+=1
            sublocus, monoexonic = sublocus
            sublocus_id = "{0}.{1}".format(superlocus_id, counter)
            attr_field = "ID={0};Parent={1};".format(sublocus_id, superlocus_id)
            if monoexonic is True:
                tag="multiexonic=false;"
            else:
                tag="multiexonic=true;"
            attr_field="{0}{1}".format(attr_field, tag)
            if sublocus.strand is None:
                substrand="."
            else:
                substrand=sublocus.strand
            sublocus_line = [ self.chrom, "locus_pipeline", "sublocus", sublocus.start, sublocus.end,
                             ".", substrand, ".", attr_field]
            sublocus_line = "\t".join([str(s) for s in sublocus_line])
            
            sublocus_lines.append(sublocus_line)
            for transcript in sublocus.transcripts:
                transcript.parent=sublocus_id
                sublocus_lines.append(str(transcript).rstrip())
                
        sublocus_lines="\n".join(sublocus_lines)
        lines=[superlocus_line]
        lines.append(sublocus_lines)
        try:
            return "\n".join(lines)
        except TypeError:
            raise TypeError(lines)
            
