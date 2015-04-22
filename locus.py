#!/usr/bin/env python3

from abstractlocus import abstractlocus
from builtins import str
#from transcript import transcript
#import random, sys

class locus(abstractlocus):
    
    '''Very basic class which holds a single transcript.'''
    
    ########### Special methods ############
    
    def __init__(self, transcript):
        
        transcript.finalize()
        self.__dict__.update(transcript.__dict__)
        self.feature="locus"
        self.splices = set(self.splices)
        self.junctions = set(self.junctions)
        self.transcripts = set()
        self.transcripts.add(transcript)
        self.id=None
        self.parent=None
        self.source = "locus_pipeline"
        
    def __str__(self):
        
        lines=[]
        if self.strand is not None:
            strand=self.strand
        else:
            strand="."
        
        if self.id is None:
            self.id = "locus:{0}{3}:{1}-{2}".format(self.chrom, self.start, self.end, strand)
            
        attr_field="ID={0};{1}".format(
                                        self.id,
                                        "{0};".format(self.parent) if self.parent is not None else ""
                                        )
        self_line = [ self.chrom, self.source, self.feature, self.start, self.end,
                             ".", strand, ".", attr_field]
        lines.append("\t".join([str(s) for s in self_line]))
        transcript=list(self.transcripts)[0]
        transcript.parent=self.id
        lines.append(str(transcript).rstrip())
        return "\n".join(lines)
        
        
    ########### Class instance methods ##############
        
 
    def add_transcript_to_locus(self, transcript):
        '''For this basic class, this method raises a NotImplementedError -
        as this container should hold only one transcript.'''
        
        raise NotImplementedError("In a locus there should be one and only one transcript!")

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
