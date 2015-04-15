import abc


class abstractlocus:
    
    @abc.abstractmethod
    def __init__(self):
        raise NotImplementedError("This is an abstract class and should not be called directly!")
    
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
            self.junctions=set.union(self.junctions, transcript.splices)
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
            cls.overlap( superlocus, transcript, flank=flank  ) > 0:
            return True
        return False 

    @abc.abstractclassmethod
    def is_intersecting(self):
        '''This class method defines how two transcript objects will be considered as overlapping.
        It is used by the BronKerbosch method, and must be implemented at the class level for each child object.'''        
        pass
    
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