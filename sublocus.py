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
        self.exons=set()

        #Copy attributes into the current dictionary
        for key in ["parent", "id", "start", "end", "chrom", "strand", "attributes"]:
            setattr(self, key, getattr(span, key))
        
        setattr( self, "monoexonic", getattr(span, "monoexonic", None)  )
        
        if span.feature=="transcript" or "RNA" in span.feature.upper():
            self.add_transcript_to_locus(span)
        
    
    def add_transcript_to_locus(self, transcript):
        if transcript is None: return
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
        self.exons = set.union(self.exons, transcript.exons)
        self.metrics[transcript.id]=dict()
    
    @property
    def splitted(self):
        '''The splitted flag indicates whether a sublocus has already been processed to produce the necessary loci.
        It must be set as a boolean flag (hence why it is coded as a property)'''
        return self.__splitted
    
    @splitted.setter
    def splitted(self,verified):
        if type(verified)!=bool:
            raise TypeError()
        self.__splitted=verified
    
    @classmethod
    def is_intersecting(cls, transcript, other):
        '''
        Implementation of the is_intersecting method. Here at the level of the sublocus,
        the intersection is seen as overlap between exons. 
        '''
          
        if transcript.id==other.id: return False # We do not want intersection with oneself
        for exon in transcript.exons:
            if any(
                   filter(
                          #Check that at least one couple of exons are overlapping
                          lambda oexon: cls.overlap( exon, oexon )>=0, 
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
        new_locus = locus(best_transcript)
        self.loci.append(new_locus)
        remaining = self.transcripts.copy()
        remaining.remove(best_transcript)
        not_intersecting=set(filter(lambda x: not self.is_intersecting(best_transcript, x), remaining ))
        
        while len(not_intersecting)>0:
            metrics=dict((tid, self.metrics[tid]) for tid in self.metrics if tid in [ trans.id for trans in not_intersecting ]  )
            best_tid,best_score=self.choose_best(metrics)
            best_transcript=next(filter(lambda x: x.id==best_tid, not_intersecting))
            new_locus = locus(best_transcript)
            self.loci.append(new_locus)
            remaining = not_intersecting.copy()
            remaining.remove(best_transcript)
            not_intersecting=set(filter(lambda x: not self.is_intersecting(best_transcript, x), remaining ))
    
        self.splitted=True
        return
    
    
    def calculate_score(self, tid):
        '''This function will try to calculate a score for a transcript various inputs.
        The scoring function must consider the following factors:
        - No. of exons
        - % of exons on the total of the exons of the sublocus
        - % of introns on the total of the intronts of the sublocus
        - % of retained introns (see has_retained_introns)
        - (if portcullis-like data available) % of high/low confidence introns
        - No. and % of CDS exons
        - CDS length
        - multiple ORFS (negative)
        - CDS exons present in the longest ORF
        - Top CDS length
        - Fraction of the transcript which is coding
        - Total CDS length
        - Difference between maximum and total CDS length (for transcripts with multiple ORFs)        
        '''
        
        transcript_instance = next(filter( lambda t: t.id==tid, self.transcripts))
        
        self.metrics[tid]["exons"]=len(transcript_instance.exons)
        self.metrics[tid]["exon_frac"] = len(set.intersection( self.exons,transcript_instance.exons   ))/len(self.exons)
        self.metrics[tid]["intron_frac"] = len(set.intersection( self.junctions,transcript_instance.junctions   ))/len(self.junctions)
        self.metrics[tid]["retained_introns"] = self.count_retained_introns(transcript_instance)
        
        self.metrics[tid]["cds"] = len(transcript_instance.cds)
        self.metrics[tid]["cds_length"] = transcript_instance.cds_length
        self.metrics["cds_fraction"] = transcript_instance.length
        
        pass
        
    
    
    def count_retained_introns(self, transcript_instance):
         
        '''This method checks the number of exons that are possibly retained introns for a given transcript.
        To perform this operation, it checks for each exon whether it exists a sublocus intron that
        is *completely* contained within a transcript exon.'''
         
        retained_introns=0
        for exon in transcript_instance.exons:
            
            if any(filter(
                          lambda junction: self.overlap(exon,junction)==junction[1]-junction[0],
                          self.junctions                          
                          )) is True:
                retained_introns+=1
        return retained_introns
    
    
    
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
                print("WARNING: multiple transcripts with the same score ({0}). I will choose one randomly.".format(", ".join(best_tid)
                                                                                                                    ) , file=sys.stderr)
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
#            print(len(self.loci))
            for slocus in sorted(self.loci): #this should function ... I have implemented the sorting in the class ...
                lines.append(str(slocus).rstrip())
                
        return "\n".join(lines)
        
        
        
        
    