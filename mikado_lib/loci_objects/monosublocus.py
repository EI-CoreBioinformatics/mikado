#!/usr/bin/env python3

import sys,os.path
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from mikado_lib.loci_objects.abstractlocus import Abstractlocus
from mikado_lib.parsers.GFF import GffLine


class Monosublocus(Abstractlocus):
    
    '''Very basic class which holds a single transcript.'''
    
    __name__ = "Monosublocus"
    
    ########### Special methods ############
    
    def __init__(self, transcript_instance, logger=None):
        
        self.counter = 0 # simple tag to avoid collisions
        super().__init__()
        self.monoexonic = transcript_instance.monoexonic # this must be defined straight away
        super().add_transcript_to_locus(transcript_instance)
        self.score = transcript_instance.score
#         self.__dict__.update(transcript_instance.__dict__)
        self.feature="Monosublocus"
        self.parent=None
        self.score = transcript_instance.score
        
#         self.source = "locus_pipeline"
        self.tid = transcript_instance.id
        self.set_logger(logger)
        self.attributes=dict()
        
    def __str__(self, print_cds=True):
        
        lines=[]

        self_line=GffLine('')
        for attr in ["chrom", 'feature','source','start','end','strand']:
            setattr(self_line,attr, getattr(self,attr))
        self_line.phase,self_line.score=None,self.score
        self_line.id="{0}_{1}".format(self.source,self.id)
        self_line.name=self.name
        self_line.parent=self.parent
        self_line.attributes.update(self.attributes)
        if "is_fragment" in self.attributes and self.attributes["is_fragment"] is False:
            del self_line.attributes["is_fragment"]
        self_line.attributes["multiexonic"]=(not self.monoexonic)
        lines.append(str(self_line))
            
        for tid in self.transcripts:
            transcript_instance=self.transcripts[tid]
            transcript_instance.source=self.source
            transcript_instance.parent=self_line.id
            self.logger.debug(self.attributes)
            for attribute in self.attributes:
                if attribute not in transcript_instance.attributes:
                    if attribute=="is_fragment" and self.attributes[attribute] is False: continue
                    transcript_instance.attributes[attribute] = self.attributes[attribute]
#                  
            lines.append(transcript_instance.__str__(print_cds=print_cds).rstrip())
            
        return "\n".join(lines)
        
        
    ########### Class instance methods ##############
        
 
    def add_transcript_to_locus(self, transcript):
        '''For this basic class, this method raises a NotImplementedError -
        as this container should hold only one transcript.'''
        
        raise NotImplementedError("In a Monosublocus there should be one and only one transcript!")


    def is_intersecting(self):
        '''Not implemented: this function makes no sense for a single-transcript container.'''
        raise NotImplementedError("Monosubloci hold a single transcript, so intersections are not calculated.")

    ########## Properties ############

    @property
    def source(self):
        return self.__source
    
    @source.setter
    def source(self,*args):
        if len(args)==0:
            args=["locus_pipeline"]
        assert len(args)==1 and type(args[0]) is str
        self.__source = args[0]
        
    @property
    def id(self):
        if self.monoexonic is True:
            addendum = "mono"
        else:
            addendum = "multi"
        if self.counter>0:
            addendum="{0}.{1}".format(addendum,self.counter)
        return "{0}.{1}".format(super().id, addendum)
        