#!/usr/bin/env python3

import sys,os.path
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from loci_objects.abstractlocus import abstractlocus
from copy import copy
from loci_objects.sublocus import sublocus
from loci_objects.monosublocus_holder import monosublocus_holder
from loci_objects.GFF import gffLine

class superlocus(abstractlocus):
    
    '''The superlocus class is used to define overlapping regions on the genome, and it receives as input
    transcript class instances.'''
    
    __name__ = "superlocus"
    available_metrics = sublocus.available_metrics
    
    ####### Special methods ############
    
    def __init__(self, transcript_instance, stranded=True, json_dict = None, purge=False ):
        
        '''The superlocus class is instantiated from a transcript_instance class, which it copies in its entirety.
        
        It will therefore have the following attributes:
        - chrom, strand, start, end
        - splices - a *set* which contains the position of each splice site
        - introns - a *set* which contains the positions of each *splice junction* (registered as 2-tuples)
        - transcripts - a *set* which holds the transcripts added to the superlocus
        
        The constructor method takes the following keyword arguments:
        - stranded    True if all transcripts inside the superlocus are required to be on the same strand
        - json_dict    Required. A dictionary with the coniguration necessary for scoring transcripts.
        - purge        Flag. If True, all loci holding only transcripts with a 0 score will be deleted from further consideration. 
        
        '''
        
        super().__init__()
        self.stranded=stranded
        self.feature=self.__name__
        if json_dict is None or type(json_dict) is not dict:
            raise ValueError("I am missing the configuration for prioritizing transcripts!")
        self.json_dict = json_dict
        
        #Dynamically load required modules
        if "modules" in self.json_dict:
            import importlib
            for mod in self.json_dict["modules"]:
                globals()[mod]=importlib.import_module(mod)
                
        #self.__dict__.update(transcript_instance.__dict__)
        self.splices = set(self.splices)
        self.introns = set(self.introns)
        self.transcripts = dict()
        self.purge = purge
        super().add_transcript_to_locus(transcript_instance)
        if self.stranded is True:
            self.strand = transcript_instance.strand
        self.available_monolocus_metrics = []
        self.available_sublocus_metrics = []
        self.set_flags()
        return

    def __str__(self, level=None):
        
        '''This function will return the desired level of children loci.
        The keyword level accepts the following four values:
        - "None" - print whatever is available.
        - "loci" - print the final loci
        - "monosubloci" - print the monosubloci
        - "subloci" - print the subloci.
        
        The function will then return the desired location in GFF3-compliant format. 
        '''

        superlocus_line=gffLine('')
        superlocus_line.chrom=self.chrom
        superlocus_line.feature=self.__name__
        superlocus_line.start,superlocus_line.end,superlocus_line.score=self.start, self.end, "."
        superlocus_line.strand=self.strand
        superlocus_line.phase, superlocus_line.score=None,None
        superlocus_line.id,superlocus_line.name=self.id, self.name

        lines=[]
        if level not in (None, "loci", "subloci", "monosubloci"):
            raise ValueError("Unrecognized level: {0}".format(level))
        
        if level=="loci" or (level is None and self.loci_defined is True):
            self.define_loci()
            if len(self.loci)>0:
                source="{0}_loci".format(self.source)
                superlocus_line.source=source
                lines.append(str(superlocus_line))
                for locus_instance in self.loci:
                    locus_instance.source=source
                    lines.append(str(locus_instance).rstrip())
        elif level=="monosubloci" or (level is None and self.monosubloci_defined is True):
            self.define_monosubloci()
            if len(self.monosubloci)>0:
                source="{0}_monosubloci".format(self.source)
                superlocus_line.source=source
                lines.append(str(superlocus_line))
                for monosublocus_instance in self.monosubloci:
                    monosublocus_instance.source=source
                    lines.append(str(monosublocus_instance).rstrip())
        elif level=="subloci" or (level is None and self.monosubloci_defined is False):
            source="{0}_subloci".format(self.source)
            superlocus_line.source=source
            lines=[str(superlocus_line)]
            self.define_subloci()
            for sublocus_instance in self.subloci:
                sublocus_instance.source=source
                lines.append(str(sublocus_instance).rstrip())
        
        if len(lines)>0:
            lines.append("###")
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
                    new_locus = superlocus(strand[0], stranded=True, json_dict=self.json_dict, purge=self.purge)
                    for cdna in strand[1:]:
                        if new_locus.in_locus(new_locus, cdna):
                            new_locus.add_transcript_to_locus(cdna)
                        else:
                            new_loci.append(new_locus)
                            new_locus = superlocus(cdna, stranded=True, json_dict=self.json_dict, purge=self.purge)
                            
                    new_loci.append(new_locus)
            for new_locus in iter(sorted(new_loci)):
                yield new_locus

    def set_flags(self):
        '''Method called by __init__ to set basic flags. These are used throughout the program to avoid unnecessary calculations.'''
        self.subloci_defined = False
        self.monosubloci_defined = False
        self.loci_defined = False
        self.monosubloci_metrics_calculated = False

    def load_cds(self, cds_dict, trust_strand=False):
        if cds_dict is None:
            return
        for tid in self.transcripts:
            self.transcripts[tid].load_cds(cds_dict, trust_strand = trust_strand)


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
            new_sublocus = sublocus(subl[0], json_dict=self.json_dict)
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
            sublocus_instance.get_metrics()

    def define_monosubloci(self):

        '''This is a wrapper method that defines the monosubloci for each sublocus.
        '''
        if self.monosubloci_defined is True:
            return
        
        self.define_subloci()
        self.monosubloci=[]
        #Extract the relevant transcripts
        for sublocus_instance in sorted(self.subloci):
            sublocus_instance.define_monosubloci(purge=self.purge)
            for ml in sublocus_instance.monosubloci:
                ml.parent = self.id
                self.monosubloci.append(ml)
            
        if self.purge is True:
            self.monosubloci = list(
                                    filter(lambda ms: ms.score > 0, self.monosubloci )
                                    )
            
        self.monosubloci = sorted(self.monosubloci)
        self.monosubloci_defined = True

    def print_monolocus_metrics(self, rower):
        '''Wrapper function to pass to a csv.DictWriter object the metrics of the transcripts in the monosubloci.'''
        
        raise NotImplementedError()

    def print_subloci_metrics(self ):
        
        '''Wrapper method to create a csv.DictWriter instance and call the sublocus.print_metrics method
        on it for each sublocus.'''
        
        self.get_sublocus_metrics()
        
        for slocus in self.subloci:
            for row in slocus.print_metrics():
                yield row

    def print_monoholder_metrics(self ):

        '''Wrapper method to create a csv.DictWriter instance and call the monosublocus_holder.print_metrics method
        on it.'''
        
        
        self.define_loci()

        #self.available_monolocus_metrics = set(self.monoholder.available_metrics)
        if self.monoholder is None:
            return ''
        for row in self.monoholder.print_metrics():
            yield row
            
    def define_loci(self):
        '''This is the final method in the pipeline. It creates a container for all the monosubloci
        (an instance of the class monosublocus_holder) and retrieves the loci it calculates internally.'''
        
        if self.loci_defined is True:
            return
        
        self.define_monosubloci()
        self.calculate_mono_metrics()

        self.loci = []        
        if self.monoholder is None:
            self.loci_defined = True
            return
        self.monoholder.define_loci(purge=self.purge)

        for locus_instance in self.monoholder.loci:
            locus_instance.parent = self.id
            self.loci.append(locus_instance)
            
        if self.purge is True:
            self.loci = list(
                             filter(lambda lc: lc.score>0, self.loci)
                             )
            
        self.loci=sorted(self.loci)
        self.loci_defined = True
        
        return
    
    def calculate_mono_metrics(self):
        '''Wrapper to calculate the metrics for the monosubloci'''
        self.monoholder = None
        
        for monosublocus_instance in sorted(self.monosubloci):
            if self.monoholder is None:
                self.monoholder = monosublocus_holder(monosublocus_instance, json_dict=self.json_dict)
            else:
                self.monoholder.add_monosublocus(monosublocus_instance)
                
        

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
        
        if monoexonic_check==0: #Both multiexonic
            for junc in transcript.introns:
                if junc in other.introns:
                    return True
        
        elif monoexonic_check==1: #One monoexonic, the other multiexonic: different subloci by definition
            return False
        
        elif monoexonic_check==2:
            if cls.overlap(
                           (transcript.start, transcript.end),
                           (other.start, other.end)
                           )>0: #A simple overlap analysis will suffice
                return True
        return False
    