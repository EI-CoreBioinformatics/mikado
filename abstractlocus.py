import abc
import random
#import sys

class abstractlocus(metaclass=abc.ABCMeta):
    
    '''This abstract class defines the basic features of any locus-like object.
    It also defines methods/properties that are needed throughout the program,
    e.g. the Bron-Kerbosch algorithm for defining cliques, or the find_retained_introns method.'''
    
    __name__ = "abstractlocus"
    
    ###### Special methods #########
    
    @abc.abstractmethod
    def __init__(self):
        self.transcripts = dict()
        self.junctions, self.exons, self.splices = set(), set(), set()
        self.start, self.end, self.strand = float("Inf"), float("-Inf"), None
        self.stranded=True
        self.initialized = False
        #raise NotImplementedError("This is an abstract class and should not be called directly!")
    
    @abc.abstractmethod
    def __str__(self):
        pass
    
    
    def __eq__(self, other):
        if type(self)!=type(other):
            return False
        for feature in ["chrom", "strand","start","end","exons","junctions","splices","stranded"]:
            if getattr(self,feature )!=getattr(other,feature):
                return False
        return True
    
    def __hash__(self):
        '''This has to be defined, otherwise abstractloci objects won't be hashable
        (and therefore operations like adding to sets will be forbidden)'''
        return super().__hash__()
    
    def __len__(self):
        return self.end-self.start+1
    
    def __lt__(self, other):
        if self.strand!=other.strand or self.chrom!=other.chrom:
            return False
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
    
    ##### Static methods #######
    
    @staticmethod
    def overlap(a, b, flank=0):
        '''This static method returns the overlap between two intervals. Values<=0 indicate no overlap.
        The optional "flank" argument (default 0) allows to expand a superlocus upstream and downstream.
        As a static method, it can be used also outside of any instance - "superlocus.overlap()" will function.
        Input: two 2-tuples of integers.
        '''
        
        left_boundary=max(a[0]-flank, b[0]-flank)
        right_boundary=min(a[1]+flank, b[1]+flank)
        
        return right_boundary - left_boundary 
    
    ##### Class methods ########

    @classmethod
    def in_locus(cls, locus_instance, transcript, flank=0):
        '''Function to determine whether a transcript should be added or not to the locus_instance.
        This is a class method, i.e. it can be used also unbound from any specific instance of the class.
        It will be possible therefore to use it to compare any locus_instance to any transcript.
        Arguments: 
        - a "locus_instance" object
        - a "transcript" object (it must possess the "finalize" method)
        - flank - optional keyword'''
        transcript.finalize()
        #We want to check for the strand only if we are considering the strand
        if locus_instance.chrom == transcript.chrom and \
            (locus_instance.stranded is False or locus_instance.strand == transcript.strand) and \
            cls.overlap( (locus_instance.start,locus_instance.end), (transcript.start,transcript.end), flank=flank  ) > 0:
                return True
        return False 

    @classmethod    
    def BronKerbosch(cls, clique, candidates, non_clique, original ):
        '''Implementation of the Bron-Kerbosch algorithm with pivot to define the subloci.
        We are using the class method "is_intersecting" to define the neighbours.
        Wiki: http://en.wikipedia.org/wiki/Bron%E2%80%93Kerbosch_algorithm '''

        pool=set.union(candidates,non_clique)

        if not any((candidates, non_clique)) or len( pool )==0:
            yield clique
            return

        pivot = random.sample( pool, 1)[0]
        pivot_neighbours = cls.neighbours(pivot, original, )
        excluded = set.difference( candidates, pivot_neighbours)

        for vertex in excluded:
            vertex_neighbours = cls.neighbours(vertex, original )
            clique_vertex = set.union(clique, set([vertex]))
            for result in cls.BronKerbosch(
                    clique_vertex,
                    set.intersection(candidates, vertex_neighbours),
                    set.intersection(non_clique, vertex_neighbours),
                    original):
                yield result
            candidates.remove(vertex)
            non_clique.add(vertex)

    @classmethod
    def neighbours( cls, vertex, graph):
        '''Function to define the vertices which are near a given vertex in the graph.'''
        return set(filter(lambda x: cls.is_intersecting(vertex, x), graph))
        
    @classmethod
    def merge_cliques(cls, cliques):
        '''This class method will merge together intersecting cliques found by the Bron-Kerbosch algorithm.
        It is therefore used to e.g. create the subloci.
        It is a somewhat naive implementation; it might be made better by looking for a more specific algorithm.
        Usually the method should be called as follows:
            - cliques = self.BronKerbosch( set(), candidates, set(), copy(candidates))
            - merged_cliques = self.merge_cliques(cliques) 
        '''
        merged_cliques = set()
        
        while len(cliques)>0:
            node=random.sample(cliques,1)[0]
            cliques.remove(node)
            new_node = set(node)
            intersecting=set()
            for merged_clique in merged_cliques:
                mc = set(merged_clique)
                if set.intersection(mc, new_node ) != set():
                    new_node=set.union(new_node, mc)
                    intersecting.add(merged_clique)
            for s in intersecting: merged_cliques.remove(s)
            merged_cliques.add(tuple(new_node))

        return merged_cliques

    @classmethod
    def choose_best(cls, metrics):
        '''Given a dictionary of metrics, this function will select the best according to the "score" item inside the dictionary itself.
        Form of the dictionary:
        
        dict( (<transcript_id>, dict( (key,val),.. (<"score">, val), ...)))
        
        It returns two items:
            transcript_id, score
        
        '''
        
        best_score,best_tid=float("-Inf"),[]
        for tid in metrics:
            score=metrics[tid]["score"]
            if score>best_score:
                best_score,best_tid=score,[tid]
            elif score==best_score:
                best_tid.append(tid)
        if len(best_tid)!=1:
            if len(best_tid)==0:
                raise ValueError("Odd. I have not been able to find the transcript with the best score: {0}".format(best_score))
            else:
#                 print("WARNING: multiple transcripts with the same score ({0}). I will choose one randomly.".format(", ".join(best_tid)
#                                                                                                                     ) , file=sys.stderr)
                best_tid=random.sample(best_tid, 1) #this returns a list
        best_tid=best_tid[0]
        return best_tid, best_score


    ####### Class instance methods  #######


    def add_transcript_to_locus(self, transcript_instance, check_in_locus = True):
        '''This method checks that a transcript is contained within the superlocus (using the "in_superlocus" class method) and
        upon a successful check extends the superlocus with the new transcript.
        More precisely, it updates the boundaries (start and end), adds the transcript to the internal "transcripts" store,
        and extends the splices and junctions with those found inside the transcript.'''
        transcript_instance.finalize()
        if self.initialized is True:
            if check_in_locus is False:
                pass
            elif not self.in_locus(self, transcript_instance):
                raise AssertionError("""Trying to merge a locus with an incompatible transcript!
                Locus: {lchrom}:{lstart}-{lend} {lstrand}
                Transcript: {tchrom}:{tstart}-{tend} {tstrand} {tid}
                """.format(
                           lchrom = self.chrom, lstart=self.start, lend = self.end, lstrand = self.strand,
                           tchrom = transcript_instance.chrom,
                           tstart =  transcript_instance.start,
                           tend =  transcript_instance.end,
                           tstrand = transcript_instance.strand,
                           tid = transcript_instance.id
                           
                           ))
        else:
            self.strand = transcript_instance.strand
            self.chrom = transcript_instance.chrom
#         if self.in_locus(self, transcript) is True:
#             if transcript.id in self.transcripts:
#                 raise KeyError("Trying to add transcript {0} to the monosublocus, but a different transcript with the same name is already present!".format(transcript.id))
        self.start = min(self.start, transcript_instance.start)
        self.end = max(self.end, transcript_instance.end)
        self.transcripts[transcript_instance.id]=transcript_instance
        self.splices.update(transcript_instance.splices)
        self.junctions.update(transcript_instance.junctions)
        self.exons.update(set(transcript_instance.exons))
        for transcript_id in self.transcripts:
            self.transcripts[transcript_id].parent=self.id

        if self.initialized is False:
            self.initialized = True
        
        return

    def find_retained_introns(self, transcript_instance):
         
        '''This method checks the number of exons that are possibly retained introns for a given transcript.
        To perform this operation, it checks for each non-CDS exon whether it exists a sublocus intron that
        is *completely* contained within a transcript exon.
        CDS exons are ignored because their retention might be perfectly valid.
        The results are stored inside the transcript instance, in the "retained_introns" tuple.'''
         
        transcript_instance.retained_introns=[]
        for exon in filter(lambda e: e not in transcript_instance.cds, transcript_instance.exons):
            #Check that the overlap is at least as long as the minimum between the exon and the intron.
            if any(filter(
                          lambda junction: self.overlap(exon,junction)>=junction[1]-junction[0],
                          self.junctions                          
                          )) is True:
                    transcript_instance.retained_introns.append(exon)
        transcript_instance.retained_introns=tuple(transcript_instance.retained_introns)

    @classmethod
    @abc.abstractmethod
    def is_intersecting(self):
        '''This class method defines how two transcript objects will be considered as overlapping.
        It is used by the BronKerbosch method, and must be implemented at the class level for each child object.'''        
        raise NotImplementedError("The is_intersecting method should be defined for each child!")

    ###### Properties #######

    @property
    def stranded(self):
        '''This property determines whether a monosublocus will consider the strand for e.g. the in_locus method.
        By default, the parameter is set to True (i.e. the loci are strand-specific).
        At the moment, the only class which modifies the parameter is the superlocus class.'''
        return self.__stranded
    
    @stranded.setter
    def stranded(self, *args):
        if len(args)==0:
            args=[True]
        stranded=args[0]
        if type(stranded)!=bool:
            raise ValueError("The stranded attribute must be boolean!")
        self.__stranded=stranded
        
    @property
    def id(self):
        return "{0}:{1}{2}:{3}-{4}".format(
                                            self.__name__,
                                            self.chrom,
                                            self.strand,
                                            self.start,
                                            self.end)