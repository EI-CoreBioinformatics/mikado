#!/usr/bin/env python3

from abstractlocus import abstractlocus
import operator
from copy import copy
from sublocus import sublocus

class superlocus(abstractlocus):
    
    '''The superlocus class is used to define overlapping regions on the genome, and it receives as input
    transcript class instances.'''
    
    ####### Special methods ############
    
    def __init__(self, transcript, stranded=True):
        
        '''The superlocus class is instantiated from a transcript class, which it copies in its entirety.
        
        It will therefore have the following attributes:
        - chrom, strand, start, end
        - splices - a *set* which contains the position of each splice site
        - junctions - a *set* which contains the positions of each *splice junction* (registered as 2-tuples)
        - transcripts - a *set* which holds the transcripts added to the superlocus'''
        
        transcript.finalize()
        self.stranded=stranded
        self.feature="superlocus"
        self.__dict__.update(transcript.__dict__)
        self.splices = set(self.splices)
        self.junctions = set(self.junctions)
        self.transcripts = dict()
        super().add_transcript_to_locus(transcript)
        
        return

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
        counter=0
        
        lines=[superlocus_line]        
        for subl in iter(sorted(self.subloci, key=operator.attrgetter("start","end") )):
            counter+=1
            #subl.id = "{0}.{1}".format(superlocus_id, counter)
            subl.parent = superlocus_id
            #attr_field = "ID={0};Parent={1};".format(sublocus_id, superlocus_id)
            lines.append(str(subl).rstrip())

        try:
            return "\n".join(lines)
        except TypeError:
            raise TypeError(lines)    

    ############ Class instance methods ############

    def split_strands(self):
        '''This method will divide the superlocus on the basis of the strand.
        The rationale is to parse a GFF file without regard for the strand, in order to find all intersecting loci;
        and subsequently break the superlocus into the different components.
        '''
        
        if self.stranded is True:
            yield self
        
        else:
            plus, minus, nones = [], [], []
            for cdna_id in self.transcripts:
                cdna=self.transcripts[cdna_id]
                if cdna.strand == "+":
                    plus.append(cdna)
                elif cdna.strand == "-":
                    minus.append(cdna)
                elif cdna.strand is None:
                    nones.append(cdna)

            for strand in plus, minus, nones:
                if len(strand)>0:
                    new_locus = superlocus(strand[0], stranded=True)
                    for cdna in strand[1:]:
                        new_locus.add_transcript_to_locus(cdna)
                    yield new_locus
                    
    def define_subloci(self):
        '''This method will define all subloci inside the superlocus.
        Steps:
            - Call the BronKerbosch algorithm to define cliques
            - Call the "merge_cliques" algorithm the merge the cliques.
            - Create "sublocus" objects from the merged cliques and store them inside the instance store "subloci"       
        '''
        
        candidates = set(self.transcripts.values()) # This will order the transcripts based on their position
        if len(candidates)==0:
            raise ValueError("This superlocus has no transcripts in it!")
        
        
        original=copy(candidates)
        
        cliques = set( tuple(clique) for clique in self.BronKerbosch(set(), candidates, set(), original))
        
        subloci = self.merge_cliques(cliques)
        
        #Now we should define each sublocus and store it in a permanent structure of the class
        self.subloci = []
    
        for subl in subloci:
            if len(subl)==0:
                continue
            
            new_sublocus = sublocus(subl[0])
            if len(subl)>1:
                for ttt in subl[1:]:
                    new_sublocus.add_transcript_to_locus(ttt)
                    
            self.subloci.append(new_sublocus)

    ############# Class methods ###########
    
    @classmethod
    def is_intersecting(cls,transcript, other):
        '''When comparing two transcripts, for the definition of subloci inside superloci we follow these rules:
        If both are multiexonic, the function verifies whether there is at least one intron in common.
        If both are monoexonic, the function verifies whether there is some overlap between them.
        If one is monoexonic and the other is not,  the function will return False by definition.        
        '''
        
        if transcript.id==other.id: return False # We do not want intersection with oneself
        monoexonic_check = len( list(filter(lambda x: x.monoexonic is True, [transcript, other]   )  )   )
        
        flag=False
        if monoexonic_check==0: #Both multiexonic
            for junc in transcript.junctions:
                if junc in other.junctions:
                    flag=True
                    break
                    
            else:
                flag=False
#             return any( filter( lambda j: j in other.junctions, transcript.junctions ) )
        
        elif monoexonic_check==1: #One monoexonic, the other multiexonic: different subloci by definition
            flag=False
        
        elif monoexonic_check==2:
            if cls.overlap(
                           (transcript.start, transcript.end),
                           (other.start, other.end)
                           )>=0: #A simple overlap analysis will suffice
                flag=True
        return flag
    