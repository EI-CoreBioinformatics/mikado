#!/usr/bin/env python3

from loci_objects.abstractlocus import abstractlocus
from builtins import str
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
#         self.__dict__.update(transcript_instance.__dict__)
        self.feature="monosublocus"
        self.parent=None
        self.source = "locus_pipeline"
        self.tid = transcript_instance.id
        
    def __str__(self):
        
        lines=[]
        if self.strand is not None:
            strand=self.strand
        else:
            strand="."
            
        attr_field="ID={0};{1}{2}".format(
                                        self.id,
                                        "{0};".format(self.parent) if self.parent is not None else "",
                                        "multiexonic={0}".format(not self.monoexonic)
                                        )
        self_line = [ self.chrom, self.source, self.feature, self.start, self.end,
                             ".", strand, ".", attr_field]
        lines.append("\t".join([str(s) for s in self_line]))
        for tid in self.transcripts:
            transcript_instance=self.transcripts[tid]

            lines.append(str(transcript_instance).rstrip())
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
        