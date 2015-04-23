#!/usr/bin/env python3

from abstractlocus import abstractlocus
#import operator
from copy import copy
from sublocus import sublocus
from monosublocus_holder import monosublocus_holder

class superlocus(abstractlocus):
    
    '''The superlocus class is used to define overlapping regions on the genome, and it receives as input
    transcript class instances.'''
    
    __name__ = "superlocus"
    
    ####### Special methods ############
    
    def __init__(self, transcript_instance, stranded=True):
        
        '''The superlocus class is instantiated from a transcript_instance class, which it copies in its entirety.
        
        It will therefore have the following attributes:
        - chrom, strand, start, end
        - splices - a *set* which contains the position of each splice site
        - junctions - a *set* which contains the positions of each *splice junction* (registered as 2-tuples)
        - transcripts - a *set* which holds the transcripts added to the superlocus'''
        
        super().__init__()
        self.stranded=stranded
        self.feature=self.__name__
        #self.__dict__.update(transcript_instance.__dict__)
        self.splices = set(self.splices)
        self.junctions = set(self.junctions)
        self.transcripts = dict()
        super().add_transcript_to_locus(transcript_instance)
        self.available_monolocus_metrics = []
        self.available_sublocus_metrics = []
        self.set_flags()
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

        superlocus_line = [self.chrom, "locus_pipeline", "superlocus", self.start, self.end, ".", strand, ".", "ID={0}".format(self.id) ]
        superlocus_line = "\t".join(str(s) for s in superlocus_line)
        lines=[superlocus_line]

        if self.loci_defined is True:
            for locus_instance in self.loci:
                lines.append(str(locus_instance).rstrip())
        elif self.monosubloci_defined is True:
            for monosublocus_instance in self.monosubloci:
                lines.append(str(monosublocus_instance).rstrip())
        else:
            self.define_subloci()
            for sublocus_instance in self.subloci:
                lines.append(str(sublocus_instance).rstrip())
        
        return "\n".join(lines)

    ############ Class instance methods ############

    def split_strands(self):
        '''This method will divide the superlocus on the basis of the strand.
        The rationale is to parse a GFF file without regard for the strand, in order to find all intersecting loci;
        and subsequently break the superlocus into the different components.
        Notice that each strand might generate more than one superlocus, if genes on a different strand link what are
        two different superloci.
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

            new_loci = []
            for strand in plus, minus, nones:
                if len(strand)>0:
                    strand = sorted(strand)
                    new_locus = superlocus(strand[0], stranded=True)
                    for cdna in strand[1:]:
                        if new_locus.in_locus(new_locus, cdna):
                            new_locus.add_transcript_to_locus(cdna)
                        else:
                            new_loci.append(new_locus)
                            new_locus = superlocus(cdna, stranded=True)
                            
                    new_loci.append(new_locus)
            for new_locus in iter(sorted(new_loci)):
                yield new_locus

    def set_flags(self):
        '''Method called by __init__ to set basic flags. These are used throughout the program to avoid unnecessary calculations.'''
        self.subloci_defined = False
        self.monosubloci_defined = False
        self.loci_defined = False
        self.monosubloci_metrics_calculated = False


    ###### Sublocus-related steps ######
                    
    def define_subloci(self):
        '''This method will define all subloci inside the superlocus.
        Steps:
            - Call the BronKerbosch algorithm to define cliques
            - Call the "merge_cliques" algorithm the merge the cliques.
            - Create "sublocus" objects from the merged cliques and store them inside the instance store "subloci"       
        '''
        
        if self.subloci_defined is True:
            return
        
        candidates = set(self.transcripts.values()) # This will order the transcripts based on their position
        if len(candidates)==0:
            raise ValueError("This superlocus has no transcripts in it!")
        
        
        original=copy(candidates)
        
        cliques = set( tuple(clique) for clique in self.BronKerbosch(set(), candidates, set(), original))
        
        subloci = self.merge_cliques(cliques)
        self.subloci = []
        #Now we should define each sublocus and store it in a permanent structure of the class
        for subl in subloci:
            if len(subl)==0:
                continue
            subl=sorted(subl)
            new_sublocus = sublocus(subl[0])
            for ttt in subl[1:]:
                new_sublocus.add_transcript_to_locus(ttt)
            new_sublocus.parent = self.id
            self.subloci.append(new_sublocus)
        self.subloci=sorted(self.subloci)
        self.subloci_defined = True

    def get_sublocus_metrics(self):
        '''Wrapper function to calculate the metrics inside each sublocus.'''
        
        self.define_subloci()
        self.sublocus_metrics = []
        for sublocus_instance in self.subloci:
            sublocus_instance.get_metrics(self)
            for metric in sublocus_instance.metrics:
                if metric not in self.sublocus_metrics:
                    self.sublocus_metrics.append(metric)

    def define_monosubloci(self):

        '''This is a wrapper method that defines the monosubloci for each sublocus.
        '''
        if self.monosubloci_defined is True:
            return
        
        self.define_subloci()
        self.monosubloci=[]
        #Extract the relevant transcripts
        for sublocus_instance in sorted(self.subloci):
            sublocus_instance.define_monosubloci()
            for ml in sublocus_instance.monosubloci:
                ml.parent = self.id
                self.monosubloci.append(ml)
            
        self.monosubloci = sorted(self.monosubloci)
        self.monosubloci_defined = True

    def print_monolocus_metrics(self, rower):
        '''Wrapper function to pass to a csv.DictWriter object the metrics of the transcripts in the monosubloci.'''
        
        raise NotImplementedError()

            
    def define_loci(self):
        '''This is the final method in the pipeline. It creates a container for all the monosubloci
        (an instance of the class monosublocus_holder) and retrieves the loci it calculates internally.'''
        
        if self.loci_defined is True:
            return
        
        self.monoholder = None
        
        for monosublocus_instance in sorted(self.monosubloci):
            if self.monoholder is None:
                self.monoholder = monosublocus_holder(monosublocus_instance)
            else:
                self.monoholder.add_monosublocus(monosublocus_instance)
                
            pass
            
        self.monoholder.define_loci()
        self.loci = []
        for locus_instance in self.monoholder.loci:
            locus_instance.parent = self.id
            self.loci.append(locus_instance)
            
        self.loci=sorted(self.loci)
        self.loci_defined = True
        
        return

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