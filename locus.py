#!/usr/bin/env python3

from abstractlocus import abstractlocus
#from transcript import transcript
#import random, sys

class locus(abstractlocus):
    
    '''Very basic class which holds a single transcript.'''
    
    def __init__(self, transcript):

        self.transcripts = set()
        self.junctions = set()
        self.splices = set()
        self.id=None
        self.parent=None
        
        super().add_transcript_to_locus(transcript)
        
 
    def add_transcript_to_locus(self, transcript):
        raise NotImplementedError("In a locus there should be one and only one transcript!")
    
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
        self_line = [ self.chrom, "locus_pipeline", "sublocus", self.start, self.end,
                             ".", strand, ".", attr_field]
        lines.append("\t".join(self_line))
        transcript=list(self.transcripts)[0]
        lines.append(str(transcript).rstrip())
        return lines
    
#     @property
#     def splitted(self):
#         return self.__splitted
#     
#     @splitted.setter
#     def splitted(self,verified):
#         if type(verified)!=bool:
#             raise TypeError()
#         self.__splitted=verified
#     
#     def init_from_line(self, gffline):
#         if gffline.feature!="sublocus":
#             raise TypeError("This is not a sublocus line!\n{0}".format(gffline))
#         
#         for key in ["chrom","start","end","strand","id","parent" ]:
#             self.__dict__[key]=gffline.__dict__[key]
#         
#         
#     def add_transcript_to_locus(self, transcript_instance):
#         '''
#         Override of the parent method. We are not adding the transcript to a generic set
#         but rather to a dictionary. 
#         '''
#         
#         assert transcript_instance.parent==self.id #check the transcript actually comes from this sublocus ..
#         super().add_transcript_to_locus(transcript_instance)
#         self.metrics[transcript_instance.id]=dict()
#         
#         return
#     
#     @classmethod
#     def is_intersecting(cls, transcript, other):
#         '''
#         Implementation of the is_intersecting method. For the loci, it should be 
#         '''
#          
#         if transcript.id==other.id: return False # We do not want intersection with oneself
#         for exon in transcript.exons:
#             if any(
#                    filter(
#                           #Check that at least one couple of exons are overlapping
#                           lambda oexon: cls.overlap( (exon.start,exon.end), (oexon.start,oexon.end) )>=0, 
#                           other.exons
#                           )
#                    ) is True:
#                 return True
#             else:
#                 return False
#     
#     def load_scores(self, score_file):
#         '''Simple mock function to load scores for the transcripts from a tab-delimited file.
#         This implementation will have to be rewritten for the final pipeline.
#         Moreover, this is potentially a disk-killer - each sublocus will have to parse the whole shebang!
#         It will probably be much better to write down everything in a SQLite DB, or parse it separately once
#         and write a "connector" to the global dictionary. Or something.
#         '''
#         
#         for line in open(score_file):
#             tid, score = line.rstrip().split()
#             score=float(score)
#             
#             if tid not in self.metrics: continue 
#             self.metrics["score"]=score
#         
#         #Using the _ as variable because it is ignored by checks for used variables 
#         if sum([1 for _ in filter(lambda tid: "score" in self.metrics[tid], self.metrics  )])!=len(self.metrics):
#             raise ValueError("I have not been able to find a score for all of the transcripts!")
#      
#     @classmethod
#     def choose_best(cls, metrics):
#         best_score,best_tid=float("-Inf"),[None]
#         for tid in metrics:
#             score=metrics[tid]["score"]
#             if score>best_score:
#                 best_score,best_tid=score,[tid]
#             elif score==best_score:
#                 best_tid.append(tid)
#         if len(best_tid)!=1:
#             if len(best_tid)==0:
#                 raise ValueError("Odd. I have not been able to find the transcript with the best score: {0}".format(best_score))
#             else:
#                 print("WARNING: multiple transcripts with the same score. I will choose one randomly.", file=sys.stderr)
#                 best_tid=random.sample(best_tid, 1) #this returns a list
#         best_tid=best_tid[0]
#         return best_tid
#         
#     def split_locus(self, min_score=float("-Inf")):
#         '''This method will split the locus if the transcript chosen for the '''
#         
#         if self.splitted is True:
#             return
#         
#         raise NotImplementedError("Still in development.")
#         
#         if any(filter(lambda x: "score" not in self.metrics[x], self.metrics)):
#             raise ValueError("Transcript scores must be loaded before this method can be used!")
#         
#         
#         #The logic here is different compared to when I was creating the subloci.
#         #Earlier, I wanted to keep *all* the transcripts of a clique
#         #Now, instead, I am retrieving all hits of a clique, 
#         best_tid=self.choose_best(self.metrics)
#         self.loci=[]
#         if self.metrics[best_tid]["score"]<min_score:
#             print("No transcript found passing the minimum score.")
#             return
#                     
#         candidates=set(self.transcripts)
#         best_transcript=set([t.id==best_tid for t in self.transcripts]  ) #Retrieve the "best" transcript
#         self.loci.append(best_transcript)
#         initial_clique=self.BronKerbosch(best_transcript, candidates, set())
#         candidates=set.difference(candidates, initial_clique)
#         while len(candidates)>0:
#             result=self.BronKerbosch(set(), candidates, set())
#             metrics=dict((tid, self.metrics[tid] ) for tid in [trans.id for trans in result]  )
#             self.loci.append(self.choose_best(metrics))
#             candidates=set.difference(candidates, result)
#             
#         self.splitted = True
#         return
#         
#         
#         
#         
#     
#     def has_retained_introns(self, transcript_instance):
#         
#         '''This method checks the number of exons that are possibly retained introns for a given transcript.
#         To perform this operation, it checks for each exon whether it exists a sublocus intron that
#         is *completely* contained within a transcript exon.'''
#         
#         retained_introns=0
#         for exon in transcript_instance.exons:
#            
#             if any(filter(
#                           lambda junction: self.overlap(exon,junction)==junction[1]-junction[0],
#                           self.junctions                          
#                           )) is True:
#                 retained_introns+=1
#         return retained_introns
#     
