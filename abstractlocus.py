import abc
import operator
import random
from copy import deepcopy as copy

class abstractlocus:
    
    __metaclass__  = abc.ABCMeta
    
    @abc.abstractmethod
    def __init__(self):
        raise NotImplementedError("This is an abstract class and should not be called directly!")
    
    @staticmethod
    def overlap(a, b, flank=0):
        '''This static method returns the overlap between two intervals. Values<=0 indicate no overlap.
        The optional "flank" argument (default 0) allows to expand a superlocus upstream and downstream.
        As a static method, it can be used also outside of any instance - "superlocus.overlap()" will function.
        Input: two 2-tuples of integers.
        '''
        
        #Removed the flexibility, as this will allow to transform the function cythonically
#         if hasattr(a, "start") and hasattr(b,"start"):
#             right_boundary=min(a.end+flank, b.end+flank)
#             left_boundary=max(a.start-flank, b.start-flank)
#         elif type(a)==type(b)==tuple and len(a)==len(b)==2:
        right_boundary=min(a[1]+flank, b[1]+flank)
        left_boundary=max(a[0]-flank, b[0]-flank)
        
        return right_boundary - left_boundary 

    def add_transcript_to_locus(self, transcript):
        '''This method checks that a transcript is contained within the superlocus (using the "in_superlocus" class method) and
        upon a successful check extends the superlocus with the new transcript.
        More precisely, it updates the boundaries (start and end), adds the transcript to the internal "transcripts" store,
        and extends the splices and junctions with those found inside the transcript.'''
        transcript.finalize()
        if self.in_locus(self, transcript) is True:
            self.start = min(self.start, transcript.start)
            self.end = max(self.end, transcript.end)
            self.transcripts.add(transcript)
            self.splices=set.union(self.splices, transcript.splices)
            for junction in transcript.junctions:
                if type(junction)!=tuple: raise TypeError(transcript.id,junction)
                self.junctions.add(junction) 
        return

    @classmethod
    def in_locus(cls, superlocus, transcript, flank=0):
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
            cls.overlap( (superlocus.start,superlocus.end), (transcript.start,transcript.end), flank=flank  ) > 0:
            return True
        return False 

    @abc.abstractmethod
    def is_intersecting(self):
        '''This class method defines how two transcript objects will be considered as overlapping.
        It is used by the BronKerbosch method, and must be implemented at the class level for each child object.'''        
        raise NotImplementedError("The is_intersecting method should be defined for each child!")
    
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
#        print(pivot_neighbours)                                                                                                                            
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
#        return clique

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
