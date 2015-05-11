#!/usr/bin/env python3

import sys,os.path
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from loci_objects.abstractlocus import abstractlocus
from loci_objects.GFF import gffLine
#from builtins import str
#from transcript import transcript
#import random, sys

class monosublocus(abstractlocus):
    
    '''Very basic class which holds a single transcript.'''
    
    __name__ = "monosublocus"
    
    ########### Special methods ############
    
    def __init__(self, transcript_instance):
        
        super().__init__()
        self.monoexonic = transcript_instance.monoexonic # this must be defined straight away
        super().add_transcript_to_locus(transcript_instance)
        self.score = transcript_instance.score
#         self.__dict__.update(transcript_instance.__dict__)
        self.feature="monosublocus"
        self.parent=None
        self.score = transcript_instance.score
#         self.source = "locus_pipeline"
        self.tid = transcript_instance.id
        
    def __str__(self, print_cds=True):
        
        lines=[]

        self_line=gffLine('')
        for attr in ["chrom", 'feature','source','start','end','strand']:
            setattr(self_line,attr, getattr(self,attr))
        self_line.phase,self_line.score=None,self.score
        self_line.id="{0}_{1}".format(self.source,self.id)
        self_line.name=self.name
        self_line.parent=self.parent
        self_line.attributes["multiexonic"]=(not self.monoexonic)
        lines.append(str(self_line))
            
        for tid in self.transcripts:
            transcript_instance=self.transcripts[tid]
            transcript_instance.source=self.source
            transcript_instance.parent=self_line.id 
            lines.append(transcript_instance.__str__(print_cds=print_cds).rstrip())
            
        return "\n".join(lines)
        
        
    ########### Class instance methods ##############
        
 
    def add_transcript_to_locus(self, transcript):
        '''For this basic class, this method raises a NotImplementedError -
        as this container should hold only one transcript.'''
        
        raise NotImplementedError("In a monosublocus there should be one and only one transcript!")


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
        
        return "{0}.{1}".format(super().id, addendum)
        