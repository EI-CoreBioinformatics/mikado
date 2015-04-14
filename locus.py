#!/usr/bin/env python3

from superlocus import superlocus
import random, sys

class locus(superlocus):
    
    '''This class '''
    
    def __init__(self, gffline):
        if gffline.feature!="sublocus":
            raise TypeError("This is not a sublocus line!\n{0}".format(gffline))
        
        for key in ["chrom","start","end","strand","id","parent" ]:
            self.__dict__[key]=gffline.__dict__[key]
        self.transcripts = set()
        self.junctions = set()
        self.splices = set()
        self.scores = None
        
    def add_transcript(self, transcript):
        '''Simple aliasing of the add_to_superlocus method.'''
        assert transcript.parent==self.id #check the transcript actually comes from this sublocus .. 
        super(locus, self).add_to_superlocus(self,transcript)
        return
    
    @classmethod
    def is_intersecting(cls, transcript, other):
        '''Override of the analogous method in the superlocus class.
        Previously, we were checking for *intronic* overlaps. Now, instead,
        we will check for *exonic* overlaps between two transcripts.        
        '''
        
        if transcript.id==other.id: return False # We do not want intersection with oneself
        for exon in transcript.exons:
            if any(filter(lambda oexon: cls.overlap( exon, oexon )>0, other.exons )) is True:
                return True
            else:
                return False
    
    def load_scores(self, score_file):
        '''Simple mock function to load scores for the transcripts from a tab-delimited file.
        This implementation will have to be rewritten for the final pipeline.'''
        
        transcripts=dict()
        
        for line in open(score_file):
            tid, score = line.rstrip().split()
            score=float(score)
            tid_in_locus = list(filter(lambda tid: tid.id==tid, self.transcripts))
            if tid_in_locus!=[]:
                transcripts[tid]=dict(
                                      ("score", score),
                                      ("transcript", tid_in_locus[0]))
        if len(transcripts)!=len(self.transcripts):
            raise ValueError("I have not been able to find a score for all of the transcripts!")
    
    
    def split_locus(self):
        '''This method will split the locus if the transcript chosen for the '''
        
        if self.scores is None:
            raise ValueError("Transcript scores must be loaded before this method can be used!")
        
        best_score = max( [self.scores[tid]["score"] for tid in self.scores ]  )
        best_tid = list(filter(lambda tid: self.scores[tid]["score"]==best_score, self.scores.keys()  ))
        if len(best_tid)!=1:
            if len(best_tid)==0:
                raise ValueError("Odd. I have not been able to find the transcript with the best score: {0}".format(best_score))
            else:
                print("WARNING: multiple choices with the same score. I will choose one randomly.", file=sys.stderr)
                best_tid=random.sample(best_tid, 1)
        best_tid=best_tid[0]
        
        raise NotImplementedError()
        