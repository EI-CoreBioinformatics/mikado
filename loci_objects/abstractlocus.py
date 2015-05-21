import os,sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import abc
import random
from copy import copy
import logging
from loci_objects.exceptions import NotInLocusError

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
        self.combined_cds_introns, self.selected_cds_introns = set(), set()
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
    
    def __repr__(self):
        return "\t".join([self.__name__,
                          self.chrom,
                          self.start,
                          self.end,
                          self.strand,
                          ",".join([t.id for t in self.transcripts]) if len(self.transcripts)>0 else "NA" ])
    
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
        a=sorted(a)
        b=sorted(b)
        
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
    def find_cliques(cls,objects, inters=None):
        '''Wrapper for the BronKerbosch algorithm, which returns the maximal cliques in the graph.
        It is the new interface for the BronKerbosch function, which is not called directly from outside this class any longer.
        The "inters" keyword provides the function used to determine whether two vertices are connected or not in the graph.
        '''
        if inters is None:
            inters = cls.is_intersecting
        assert hasattr(inters, "__call__")
                
        graph = dict()
        for obj in objects:
            graph[obj]=set( other_obj for other_obj in objects if (obj!=other_obj) and inters(obj,other_obj) is True )

        final_cliques = []
        candidates = set(graph.keys())
        non_clique=set()
        clique=set()  
        
        for vertex in graph:
            neighbours = graph[vertex]
            vertex_clique = set.union(clique, set([vertex])) # this is identical to neighbours
            vertex_candidates = set.intersection( candidates, neighbours )
            vertex_non_clique = set.intersection( non_clique, neighbours)
            final_cliques.extend(cls.BronKerbosch(graph, vertex_clique, vertex_candidates, vertex_non_clique, []))
            candidates.remove(vertex)
            non_clique.add(vertex)

        return final_cliques


    @classmethod    
    def BronKerbosch(cls, graph, clique, candidates, non_clique, final_cliques):
        '''Implementation of the Bron-Kerbosch algorithm with pivot to define the subloci.
        We are using the class method "is_intersecting" to define the neighbours.
        Wiki: http://en.wikipedia.org/wiki/Bron%E2%80%93Kerbosch_algorithm
        The function takes five arguments:
            - graph            The original graph. It must be a *dictionary* of the form dict[vertex]=set(intersecting vertices) - i.e. the neighbours 
            - clique            the current clique. Initialize to an empty set. 
            - candidates        the elements to be analysed. This should be a set.
            - non_clique        Elements which have already been analysed and determined not to be in the clique
            - final_clique    The list of the final cliques. It should be initialised to an empty list.
        '''

        pool=set.union(candidates,non_clique)
        if (len(candidates)==0 and len(non_clique)==0) or len(pool)==0:
        # if not any((candidates, non_clique)) or len( pool )==0:
            final_cliques.append(clique)
        else:
            pivot = random.sample( pool, 1)[0]
            excluded = set.difference( candidates, graph[pivot])

            for vertex in excluded:
                vertex_neighbours = graph[vertex]
                cls.BronKerbosch(graph, clique.union(set([vertex])),
                                 candidates.intersection(vertex_neighbours),
                                 non_clique.intersection(vertex_neighbours),
                                 final_cliques)
                
                candidates.remove(vertex)
                non_clique.add(vertex)
                
        return final_cliques

    @classmethod
    def merge_cliques(cls, cliques):
        '''This class method will merge together intersecting cliques found by the Bron-Kerbosch algorithm.
        It is therefore used to e.g. create the subloci.
        It is a somewhat naive implementation; it might be made better by looking for a more specific algorithm.
        Usually the method should be called as follows:
            - cliques = self.find_cliques(objects, inters=intersecting_function)
            - merged_cliques = self.merge_cliques(cliques) 
        '''
        merged_cliques = set()

        cliques=sorted(cliques, key=len, reverse=True)
#        print("# of cliques:", len(cliques))
        
        while len(cliques)>0:
            node=cliques[0]
            cliques.remove(node)
            new_node = set(node)
            intersecting=set()
            to_remove=[]
            for index in range(len(cliques)):
                other_clique=cliques[index]
                if set.intersection(other_clique,new_node) != set():
                    new_node.update(other_clique)
                    to_remove.append(index)
            cliques=[ v for i,v in enumerate(cliques) if i not in to_remove] #Remove cliques we have already assigned
            
            for merged_clique in merged_cliques:
                mc = set(merged_clique)
                if set.intersection(mc, new_node ) != set():
                    new_node.update(mc)
                    intersecting.add(merged_clique)
            for s in intersecting:
                merged_cliques.remove(s)
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
        if len(transcripts)==0:
            raise ValueError("No transcripts provided!")
        
        scores=dict((tid,transcripts[tid].score) for tid in transcripts  )
        
        selected_tid = list(filter( lambda tid: scores[tid]==max(scores.values()), scores ))
        if len(selected_tid)!=1:
            if len(selected_tid)==0:
                raise ValueError("Odd. I have not been to select the best transcript among these:\n{0}".format(
                                                                                                               "\n".join(["{0}\t{1}".format(key, scores[key]) for key in scores ]  )  
                                                                                                               ))
            else:
                selected_tid = random.sample(selected_tid,1)
        selected_tid=selected_tid[0]
        return selected_tid


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
                raise NotInLocusError("""Trying to merge a locus with an incompatible transcript!
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

        self.combined_cds_introns=set.union(self.combined_cds_introns,transcript_instance.combined_cds_introns)
        self.selected_cds_introns.update(transcript_instance.selected_cds_introns)

        self.exons.update(set(transcript_instance.exons))
        for transcript_id in self.transcripts:
            self.transcripts[transcript_id].parent=self.id

        if self.initialized is False:
            self.initialized = True
        self.source=transcript_instance.source
        assert transcript_instance.id in self.transcripts
        return

    def remove_transcript_from_locus(self, tid):
        if tid not in self.transcripts:
            raise KeyError("Transcript {0} is not present in the locus.".format(tid))
        
        if len(self.transcripts)==1:
            self.transcripts = dict()
            self.introns, self.exons, self.splices = set(), set(), set()
            self.cds_introns, self.selected_cds_introns = set(), set()
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
            other_cds_introns = set.union(*[ set(self.transcripts[otid].combined_cds_introns) for otid in self.transcripts if otid!=tid  ])
            cds_introns_to_remove = set.difference(set(self.transcripts[tid].combined_cds_introns), other_cds_introns)
            self.combined_cds_introns.difference_update(cds_introns_to_remove) 

            #Remove excess selected_cds_introns
            other_selected_cds_introns = set.union(*[ set(self.transcripts[otid].selected_cds_introns) for otid in self.transcripts if otid!=tid  ])
            selected_cds_introns_to_remove = set.difference(set(self.transcripts[tid].selected_cds_introns), other_selected_cds_introns)
            self.selected_cds_introns.difference_update(selected_cds_introns_to_remove) 

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
