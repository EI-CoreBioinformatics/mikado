from abstractlocus import abstractlocus
import random,sys,operator
from locus import locus

class sublocus(abstractlocus):
    
    '''
    The sublocus class is created either by the superlocus class during the subloci definition, or directly using a G(T|F)line-like object.
    It is used to define the final loci.
    '''
    
    def __init__(self, span):
        
        '''This '''
        
        self.transcripts = set()
        self.fixedSize=True if span.feature=="sublocus" else False
        self.feature="sublocus"
        self.splices = set()
        self.junctions = set()
        self.metrics=dict()
        self.splitted=False

        #Copy attributes into the current dictionary
        for key in ["parent", "id", "start", "end", "chrom", "strand", "attributes"]:
            setattr(self, key, getattr(span, key))
        
        setattr( self, "monoexonic", getattr(span, "monoexonic", None)  )
        
        if span.feature=="transcript" or "RNA" in span.feature.upper():
            self.add_transcript_to_locus(span)
        
    
    def add_transcript_to_locus(self, transcript):
        transcript.finalize()
        if self.monoexonic is None:
            self.monoexonic = transcript.monoexonic
        elif self.monoexonic!=transcript.monoexonic:
            raise ValueError("Sublocus and transcript are not compatible!\n{0}\t{1}\t{2}\t{3}\t{4}\n{5}".format(self.chrom,
                                                                                                           self.start,
                                                                                                           self.end,
                                                                                                           self.strand,
                                                                                                           self.monoexonic,
                                                                                                           transcript))
        super().add_transcript_to_locus(transcript)
        self.metrics[transcript.id]=dict()
    
    @property
    def splitted(self):
        return self.__splitted
    
    @splitted.setter
    def splitted(self,verified):
        if type(verified)!=bool:
            raise TypeError()
        self.__splitted=verified
    
#     @classmethod
#     def is_intersecting(cls, first, second):
#                 
#         
#         ''' For the definition of the subloci we are interested in exonic overlaps, not intronic overlaps.'''
#         
#         for exon in first.exons:
#             if any(filter( lambda oexon: cls.overlap(exon, oexon), second.exons  ) ):
#                 return True
#         return False
# #        raise NotImplementedError()

    @classmethod
    def is_intersecting(cls, transcript, other):
        '''
        Implementation of the is_intersecting method. For the loci, it should be 
        '''
          
        if transcript.id==other.id: return False # We do not want intersection with oneself
        for exon in transcript.exons:
            if any(
                   filter(
                          #Check that at least one couple of exons are overlapping
                          lambda oexon: cls.overlap( (exon.start,exon.end), (oexon.start,oexon.end) )>=0, 
                          other.exons
                          )
                   ) is True:
                return True
        
        return False
        
    
    def define_loci(self):
         
        self.loci=[]
        best_tid,best_score=self.choose_best(self.metrics)
        best_transcript=next(filter(lambda x: x.id==best_tid, self.transcripts))
        best_transcript.score=best_score
        remaining = self.transcripts.copy()
        remaining.remove(best_transcript)
        not_intersecting=set(filter(lambda x: not self.is_intersecting(best_transcript, x), remaining ))
        
        while len(not_intersecting)>0:
            metrics=dict((tid, self.metrics[tid]) for tid in self.metrics if tid in [ trans.id for trans in not_intersecting ]  )
            best_tid=self.choose_best(metrics)
            best_transcript=next(filter(lambda x: x.id==best_tid, not_intersecting))
            new_locus = locus(best_transcript)
            self.loci.append(new_locus)
            remaining = not_intersecting.copy()
            remaining.remove(best_transcript)
            not_intersecting=set(filter(lambda x: not self.is_intersecting(best_transcript, x), remaining ))
    
        self.splitted=True
        return
    
    def load_scores(self, scores):
        '''Simple mock function to load scores for the transcripts from a tab-delimited file.
        This implementation will have to be rewritten for the final pipeline.
        Moreover, this is potentially a disk-killer - each sublocus will have to parse the whole shebang!
        It will probably be much better to write down everything in a SQLite DB, or parse it separately once
        and write a "connector" to the global dictionary. Or something.
        '''
        
        for tid in filter(lambda tid: tid in self.metrics, scores):
            self.metrics[tid]["score"]=scores[tid]
        
        #Using the _ as variable because it is ignored by checks for used variables
        try: 
            if sum([1 for _ in filter(lambda tid: "score" in self.metrics[tid], self.metrics  )])!=len(self.metrics):
                raise ValueError("I have not been able to find a score for all of the transcripts!")
        except TypeError as err:
            raise TypeError("{0}\n{1}".format(err, self.metrics) )
            
    
    @classmethod
    def choose_best(cls, metrics):
        best_score,best_tid=float("-Inf"),[None]
        for tid in metrics:
            score=metrics[tid]["score"]
            if score>best_score:
                best_score,best_tid=score,[tid]
            elif score==best_score:
                best_tid.append(tid)
        if len(best_tid)!=1:
            if len(best_tid)==0:
                raise ValueError("Odd. I have not been able to find the transcript with the best score: {0}".format(best_score))
            else:
                print("WARNING: multiple transcripts with the same score. I will choose one randomly.", file=sys.stderr)
                best_tid=random.sample(best_tid, 1) #this returns a list
        best_tid=best_tid[0]
        return best_tid, best_score
    
    def __str__(self):
        
        if self.strand is not None:
            strand=self.strand
        else:
            strand="."
        
        lines=[]
        if self.splitted is False:
        
            lines=[]
            attr_field="ID={0};Parent={1};multiexonic={2}".format(
                                                                  self.id,
                                                                  self.parent,
                                                                  str(not self.monoexonic)
                                                                  )
                            
            self_line = [ self.chrom, "locus_pipeline", "sublocus", self.start, self.end,
                             ".", strand, ".", attr_field]
            lines.append("\t".join([str(s) for s in self_line]))

        
            transcripts = iter(sorted(self.transcripts, key=operator.attrgetter("start", "end")))
            for transcript in transcripts:
                transcript.parent=self.id
                lines.append(str(transcript).rstrip())
                
        else:
            for slocus in sorted(self.loci): #this should function ... I have implemented the sorting in the class ...
                lines.append(str(slocus).rstrip())
                
        return "\n".join(lines)
        
        
        
        
    