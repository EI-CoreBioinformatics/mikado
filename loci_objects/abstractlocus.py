import os,sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import abc
import random
from copy import copy
import logging

class abstractlocus(metaclass=abc.ABCMeta):
    
    '''This abstract class defines the basic features of any locus-like object.
    It also defines methods/properties that are needed throughout the program,
    e.g. the Bron-Kerbosch algorithm for defining cliques, or the find_retained_introns method.'''

    
    __name__ = "abstractlocus"
    
    ###### Special methods #########
    
    @abc.abstractmethod
    def __init__(self, logger=None):
        self.transcripts = dict()
        self.introns, self.exons, self.splices = set(), set(), set()
        #Consider only the CDS part
        self.cds_introns, self.best_cds_introns = set(), set()
        self.start, self.end, self.strand = float("Inf"), float("-Inf"), None
        self.stranded=True
        self.initialized = False
        #raise NotImplementedError("This is an abstract class and should not be called directly!")
        if logger is None:
            self.logger = logging.Logger("{0}_logger".format(self.__name__))
        else:
            self.logger=logger
        
    
    @abc.abstractmethod
    def __str__(self):
        pass
    
    
    def __eq__(self, other):
        if type(self)!=type(other):
            return False
        for feature in ["chrom", "strand","start","end","exons","introns","splices","stranded"]:
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
    
    @staticmethod
    def evaluate(param, conf):
    
        '''This static method evaluates a single parameter using the requested operation from the JSON dict file.'''
        
        if conf["operator"]=="eq":
            return float(param)==float(conf["value"])
        elif conf["operator"]=="ne":
            return float(param)!=float(conf["value"])
        elif conf["operator"]=="gt":
            return float(param)>float(conf["value"])
        elif conf["operator"]=="lt":
            return float(param)<float(conf["value"])
        elif conf["operator"]=="ge":
            return float(param)>=float(conf["value"])
        elif conf["operator"]=="le":
            return float(param)<=float(conf["value"])
        elif conf["operator"]=="in":
            return param in conf["value"]
        elif conf["operator"]=="not in":
            return param not in conf["value"]
    
    
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
    def BronKerbosch(cls, clique, candidates, non_clique, original, inters = None, neighbours = None ):
        '''Implementation of the Bron-Kerbosch algorithm with pivot to define the subloci.
        We are using the class method "is_intersecting" to define the neighbours.
        Wiki: http://en.wikipedia.org/wiki/Bron%E2%80%93Kerbosch_algorithm
        The function takes four arguments and two keyword arguments:
            - clique            the cliques already found. This should be initialised as an empty set.
            - candidates        the elements to be analysed. This should be a set.
            - non_clique        Elements which have already been analysed and determined not to be cliques
            - original          A copy of the candidates set. Needed for the iteration.
            
        Keyword arguments:
             - inters=None      The intersection function to determine the cliques. It defaults to the class "is_intersecting" method.
             - neighbours=None  The function used to determine the neighbours of an element. It defaults to the class "neighbours" method.
        '''

        pool=set.union(candidates,non_clique)

        if not any((candidates, non_clique)) or len( pool )==0:
            yield clique
            return
        
        if inters is None:
            inters = cls.is_intersecting
        if neighbours is None:
            neighbours = cls.neighbours
        
        #Check the functions are actually functions
        assert hasattr(inters, "__call__") and hasattr(neighbours, "__call__")

        pivot = random.sample( pool, 1)[0]
        pivot_neighbours = neighbours(pivot, original, inters = inters)
        excluded = set.difference( candidates, pivot_neighbours)

        for vertex in excluded:
            vertex_neighbours = neighbours(vertex, original, inters = inters )
            clique_vertex = set.union(clique, set([vertex]))
            for result in cls.BronKerbosch(
                    clique_vertex,
                    set.intersection(candidates, vertex_neighbours),
                    set.intersection(non_clique, vertex_neighbours),
                    original,
                    neighbours = neighbours,
                    inters = inters
                    ):
                yield result
            candidates.remove(vertex)
            non_clique.add(vertex)

    @classmethod
    def neighbours( cls, vertex, graph, inters = None):
        if inters is None:
            inters = cls.is_intersecting
        
        '''Function to define the vertices which are near a given vertex in the graph.'''
        return set(filter(lambda x: inters(vertex, x), graph))
        
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
    def choose_best(cls, transcripts):
        '''Given a dictionary of metrics, this function will select the best according to the "score" item inside the dictionary itself.
        Form of the dictionary:
        
        dict( (<transcript_id>, dict( (key,val),.. (<"score">, val), ...)))
        
        It returns two items:
            transcript_id, score
        
        '''
        
        best_score,best_tid=float("-Inf"),[]
        if len(transcripts)==0:
            raise ValueError("Empty dictionary!")
        for tid in transcripts:
            score=transcripts[tid].score
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
        return best_tid


    ####### Class instance methods  #######


    def add_transcript_to_locus(self, transcript_instance, check_in_locus = True):
        '''This method checks that a transcript is contained within the superlocus (using the "in_superlocus" class method) and
        upon a successful check extends the superlocus with the new transcript.
        More precisely, it updates the boundaries (start and end), adds the transcript to the internal "transcripts" store,
        and extends the splices and introns with those found inside the transcript.'''
        transcript_instance.finalize()
        if self.initialized is True:
            if check_in_locus is False:
                pass
            elif not self.in_locus(self, transcript_instance):
                raise AssertionError("""Trying to merge a locus with an incompatible transcript!
                Locus: {lchrom}:{lstart}-{lend} {lstrand} [{stids}]
                Transcript: {tchrom}:{tstart}-{tend} {tstrand} {tid}
                """.format(
                           lchrom = self.chrom, lstart=self.start, lend = self.end, lstrand = self.strand,
                           tchrom = transcript_instance.chrom,
                           tstart =  transcript_instance.start,
                           tend =  transcript_instance.end,
                           tstrand = transcript_instance.strand,
                           tid = transcript_instance.id,
                           stids = ", ".join(list( self.transcripts.keys()) ),
                           
                           ))
        else:
            self.strand = transcript_instance.strand
            self.chrom = transcript_instance.chrom
#         if self.in_locus(self, transcript) is True:
#             if transcript.id in self.transcripts:
#                 raise KeyError("Trying to add transcript {0} to the monosublocus, but a different transcript with the same name is already present!".format(transcript.id))
        self.start = min(self.start, transcript_instance.start)
        self.end = max(self.end, transcript_instance.end)
        self.transcripts[transcript_instance.id]=copy(transcript_instance)
        self.splices.update(transcript_instance.splices)
        self.introns.update(transcript_instance.introns)
        
        self.cds_introns.update(transcript_instance.cds_introns)
        self.best_cds_introns.update(transcript_instance.best_cds_introns)
        
        self.exons.update(set(transcript_instance.exons))
        for transcript_id in self.transcripts:
            self.transcripts[transcript_id].parent=self.id

        if self.initialized is False:
            self.initialized = True
        self.source=transcript_instance.source
        return

    def remove_transcript_from_locus(self, tid):
        if tid not in self.transcripts:
            raise KeyError("Transcript {0} is not present in the locus.".format(tid))
        
        if len(self.transcripts)==1:
            self.transcripts = dict()
            self.introns, self.exons, self.splices = set(), set(), set()
            self.cds_introns, self.best_cds_introns = set(), set()
            self.start, self.end, self.strand = float("Inf"), float("-Inf"), None
            self.stranded=True
            self.initialized = False
        
        else:
            
            self.end=max(self.transcripts[t].end for t in self.transcripts if t!=tid)
            self.start=min(self.transcripts[t].start for t in self.transcripts if t!=tid)
            
            #Remove excess exons
            other_exons = [set(self.transcripts[otid].exons) for otid in self.transcripts if otid!=tid]
            other_exons = set.union( *other_exons )
            exons_to_remove = set.difference(set(self.transcripts[tid].exons), other_exons)
            self.exons.difference_update(exons_to_remove) 
            
            #Remove excess introns
            other_introns = set.union(*[ set(self.transcripts[otid].introns) for otid in self.transcripts if otid!=tid  ])
            introns_to_remove = set.difference(set(self.transcripts[tid].introns), other_introns)
            self.introns.difference_update(introns_to_remove) 

            #Remove excess cds introns
            other_cds_introns = set.union(*[ set(self.transcripts[otid].cds_introns) for otid in self.transcripts if otid!=tid  ])
            cds_introns_to_remove = set.difference(set(self.transcripts[tid].cds_introns), other_cds_introns)
            self.cds_introns.difference_update(cds_introns_to_remove) 

            #Remove excess best_cds_introns
            other_best_cds_introns = set.union(*[ set(self.transcripts[otid].best_cds_introns) for otid in self.transcripts if otid!=tid  ])
            best_cds_introns_to_remove = set.difference(set(self.transcripts[tid].best_cds_introns), other_best_cds_introns)
            self.best_cds_introns.difference_update(best_cds_introns_to_remove) 

            #Remove excess splices
            other_splices = set.union(*[ self.transcripts[otid].splices for otid in self.transcripts if otid!=tid  ])
            splices_to_remove = set.difference(self.transcripts[tid].splices, other_splices)
            self.splices.difference_update(splices_to_remove) 
                
            del self.transcripts[tid]
            for tid in self.transcripts:
                self.transcripts[tid].parent = self.id


    def find_retained_introns(self, transcript_instance):
         
        '''This method checks the number of exons that are possibly retained introns for a given transcript.
        To perform this operation, it checks for each non-CDS exon whether it exists a sublocus intron that
        is *completely* contained within a transcript exon.
        CDS exons are ignored because their retention might be perfectly valid.
        The results are stored inside the transcript instance, in the "retained_introns" tuple.'''
         
        transcript_instance.retained_introns=[]
        for exon in filter(lambda e: e not in transcript_instance.combined_cds, transcript_instance.exons):
            #Check that the overlap is at least as long as the minimum between the exon and the intron.
            if any(filter(
                          lambda junction: self.overlap(exon,junction)>=junction[1]-junction[0],
                          self.introns                          
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
        
    @property
    def name(self):
        return self.id