from abstractlocus import abstractlocus
#import random
from copy import copy
from monosublocus import monosublocus
from transcript import transcript

class sublocus(abstractlocus):
    
    '''
    The sublocus class is created either by the superlocus class during the subloci definition, or directly using a G(T|F)line-like object.
    It is used to define the final monosubloci.
    '''
    
    __name__ = "sublocus"
    
    ################ Class special methods ##############
    
    def __init__(self, span, source=None):
        
        '''This class takes as input a "span" feature - e.g. a gffLine or a transcript_instance. 
        The span instance should therefore have such attributes as chrom, strand, start, end, attributes. '''
        super().__init__()
        
        self.fixedSize=True if span.feature=="sublocus" else False
        self.feature="sublocus"
        if source is not None:
            self.source = source
        else:
            self.source = "locus_pipeline"
        self.metrics=dict()
        self.splitted=False
        self.metrics_calculated=False #Flag to indicate that we have not calculated the metrics for the transcripts
        self.available_metrics=[] # List to retrieve the available metrics. Calculated at runtime.
        setattr( self, "monoexonic", getattr(span, "monoexonic", None)  )
        assert hasattr(self, "monoexonic")

        if type(span) is transcript:
            self.add_transcript_to_locus(span)
        else:
            for key in ["parent", "start", "end", "chrom", "strand", "attributes"]:
                setattr(self, key, getattr(span, key))
        
        
    def __str__(self):
        #print(self.id, len(self.transcripts), file=sys.stderr)
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
                            
            self_line = [ self.chrom, self.source, "sublocus", self.start, self.end,
                     ".", strand, ".", attr_field]
            lines.append("\t".join([str(s) for s in self_line]))
        
            for tid in sorted(self.transcripts, key=lambda tid: self.transcripts[tid]):
                lines.append(str(self.transcripts[tid]).rstrip())
                
        else:
            for slocus in sorted(self.monosubloci): #this should function ... I have implemented the sorting in the class ...
                lines.append(str(slocus).rstrip())
                
        return "\n".join(lines)

    
    ########### Class instance methods #####################
    
    def add_transcript_to_locus(self, transcript_instance):
        
        '''This is an override of the original method, as at the sublocus stage we need to accomplish a couple of things more:
        - check that transcripts added to the sublocus are either all monoexonic or all multiexonic
        - change the id of the transcripts to  
        '''
        
        if transcript_instance is None: return
        assert hasattr(self, "monoexonic")
        if self.initialized is False:
            self.monoexonic = transcript_instance.monoexonic
        elif self.monoexonic!=transcript_instance.monoexonic:
            raise ValueError("Sublocus and transcript are not compatible!\n{0}\t{1}\t{2}\t{3}\t{4}\n{5}".format(self.chrom,
                                                                                                           self.start,
                                                                                                           self.end,
                                                                                                           self.strand,
                                                                                                           self.monoexonic,
                                                                                                           transcript))

        super().add_transcript_to_locus(transcript_instance)
        
        #Update the id
        self.metrics[transcript_instance.id]=dict()

    def define_monosubloci(self):
        '''This function retrieves the best non-overlapping transcripts inside the sublocus, according to the score
        calculated by calculate_scores (explicitly called inside the method).'''
        
        self.monosubloci=[]
        self.calculate_scores()
        remaining = self.transcripts.copy()
        
        while len(remaining)>0:
            metrics = dict( 
                           (tid, self.metrics[tid]) for tid in self.metrics if tid in remaining
                            )
            best_tid,best_score=self.choose_best(metrics)
            best_transcript = remaining[best_tid]
            best_transcript.score = best_score
            new_locus = monosublocus(best_transcript)
            self.monosubloci.append(new_locus)
            remaining = remaining.copy()
            del remaining[best_tid]
            for tid in list(remaining.keys()):
                if self.is_intersecting(best_transcript, remaining[tid] ):
                    del remaining[tid]
    
        self.splitted=True
        return
    
    def calculate_metrics(self, tid):
        '''This function will calculate the metrics which will be used to derive a score for a transcript.
        The different attributes of the transcript will be stored inside the transcript class itself,
        to be more precise, into an internal dictionary ("metrics").
         
        The scoring function must consider the following factors:
        - "exons":              No. of exons 
        - "exon_frac":          % of exons on the total of the exons of the sublocus
        - "intron_frac":        % of introns on the total of the intronts of the sublocus
        - "retained_introns":   no. of retained introns (see has_retained_introns)
        - "retained_frac":      % of cdna_length that is in retained introns 
        - "cds_length":         length of the CDS
        - "cds_fraction":       length of the CDS/length of the cDNA
        - "internal_cds_num":   number of internal CDSs. 1 is top, 0 is worst, each number over 1 is negative.             
        - "max_internal_cds_exon_num":   number of CDS exons present in the longest ORF.
        - "max_internal_cds_length":     length of the greatest CDS
        - "max_internal_cds_fraction":   fraction of the cDNA which is in the maximal CDS
        - "cds_not_maximal":            length of CDS *not* in the maximal ORF
        - "cds_not_maximal_fraction"    fraction of CDS *not* in the maximal ORF 
        '''
    
        transcript_instance = self.transcripts[tid]
        transcript_instance.finalize() # The transcript must be finalized before we can calculate the score.
        
        transcript_instance.metrics["exons"]=len(transcript_instance.exons)
        transcript_instance.metrics["exon_frac"] = len(set.intersection( self.exons,transcript_instance.exons   ))/len(self.exons)
        if len(self.junctions)>0:
            transcript_instance.metrics["intron_frac"] = len(set.intersection( self.junctions,transcript_instance.junctions   ))/len(self.junctions)
        else:
            transcript_instance.metrics["intron_frac"]=0
        self.find_retained_introns(transcript_instance)
        transcript_instance.metrics["retained_introns"] = len(transcript_instance.retained_introns)
        transcript_instance.metrics["retained_frac"]=sum(e[1]-e[0]+1 for e in transcript_instance.retained_introns)/transcript_instance.cdna_length
        transcript_instance.metrics["cds_exons"] = len(transcript_instance.cds)
        transcript_instance.metrics["cds_length"] = transcript_instance.cds_length
        transcript_instance.metrics["cds_fraction"] = transcript_instance.cds_length/transcript_instance.cdna_length
        transcript_instance.metrics["utr_fraction"] = 1-transcript_instance.metrics["cds_fraction"] 
        transcript_instance.metrics["internal_cds_num"] = transcript_instance.internal_cds_num
        transcript_instance.metrics["max_internal_cds_num"] = len(transcript_instance.max_internal_cds)
        transcript_instance.metrics["max_internal_cds_length"] = transcript_instance.max_internal_cds_length
        if transcript_instance.cds_length>0:
            transcript_instance.metrics["max_internal_cds_fraction"] = transcript_instance.max_internal_cds_length/transcript_instance.cds_length
        else:
            transcript_instance.metrics["max_internal_cds_fraction"] = 0
        transcript_instance.metrics["cds_not_maximal"] = transcript_instance.cds_length - transcript_instance.max_internal_cds_length
        if transcript_instance.max_internal_cds_length>0:
            transcript_instance.metrics["cds_not_maximal_fraction"] = (transcript_instance.cds_length - transcript_instance.max_internal_cds_length)/transcript_instance.cds_length
        else:
            transcript_instance.metrics["cds_not_maximal_fraction"] = 0
        transcript_instance.metrics["cdna_length"] = transcript_instance.cdna_length
        for metric in transcript_instance.metrics:
            if metric not in self.available_metrics:
                self.available_metrics.append(metric)

        self.transcripts[tid]=transcript_instance
        
    def find_retained_introns(self, transcript_instance):
         
        '''This method checks the number of exons that are possibly retained introns for a given transcript.
        To perform this operation, it checks for each non-CDS exon whether it exists a sublocus intron that
        is *completely* contained within a transcript exon.
        CDS exons are ignored because their retention might be perfectly valid.
        The results are stored inside the transcript instance, in the "retained_introns" tuple.'''
         
        transcript_instance.retained_introns=[]
        for exon in filter(lambda e: e not in transcript_instance.cds, transcript_instance.exons):
            #Check that the overlap is at least as long as the minimum between the exon and the intron.
            if any(filter(
                          lambda junction: self.overlap(exon,junction)>=junction[1]-junction[0],
                          self.junctions                          
                          )) is True:
                    transcript_instance.retained_introns.append(exon)
#                     print("Retained:", transcript_instance.id, exon, file=sys.stderr)
        transcript_instance.retained_introns=tuple(transcript_instance.retained_introns)
            
    def load_scores(self, scores):
        '''Simple mock function to load scores for the transcripts from an external dictionary.
        '''
        
        for tid in filter(lambda tid: tid in self.metrics, scores):
            self.metrics[tid]["score"]=scores[tid]
        
        #Using the _ as variable because it is ignored by checks for used variables
        try: 
            if sum([1 for _ in filter(lambda tid: "score" in self.metrics[tid], self.metrics  )])!=len(self.metrics):
                raise ValueError("I have not been able to find a score for all of the transcripts!")
        except TypeError as err:
            raise TypeError("{0}\n{1}".format(err, self.metrics) )
            
            
    def calculate_scores(self, maximum_utr=3, minimum_cds_fraction=0.2, order=("exons", "cdna_length", "intron_frac", "utr_fraction", )):
        '''
        Function to calculate a score for each transcript, given the metrics derived
         with the calculate_metrics method.
         Keyword arguments:
         
         - maximum_utr                Maximum number of UTR exons
         - minimum_cds_fraction       How much of the cdna should be CDS, at a minimum
         '''
        
        self.get_metrics()
        
        def keyfunction(transcript_instance, order=[]   ):
            if len(order)==0:
                raise ValueError("No information on how to sort!")
            return tuple([transcript_instance.metrics[x] for x in order])
        
        current_score=0
        
        for tid in sorted( self.transcripts, key = lambda tid: keyfunction(self.transcripts[tid], order=order), reverse=True  ):
            transcript_instance = self.transcripts[tid]
            if len(transcript_instance.utr)>maximum_utr or self.metrics[transcript_instance.id]["cds_fraction"]<minimum_cds_fraction:
                score=float("-Inf")
            else:
                score=current_score
                current_score+=1
            self.metrics[tid]["score"]=score
            self.transcripts[tid].score=score

        self.scores_calculated=True
         
        return
    
    
    def print_metrics(self, rower):
        
        '''This class takes as input a csv.DictWriter class, which it uses to print out a table of each transcript metrics.'''
        
        #Check that rower is an instance of the csv.DictWriter class
        if not hasattr(rower, "fieldnames") or not hasattr(rower, "writerow") or not hasattr(rower, "writer"):
            raise AttributeError("Invalid rower provided: should be a DictWriter instance, instead is {0}".format(type(rower)))
        
        self.get_metrics() 
        
        #The rower is an instance of the DictWriter class from the standard CSV module
        
        for tid in sorted(self.transcripts, key=lambda tid: self.transcripts[tid] ):
            row=dict()
            for key in rower.fieldnames:
                if not key in self.metrics[tid]:
                    row[key]=getattr(self.transcripts[tid], key, "NA")
                else:
                    row[key]=self.metrics[tid][key]
                if type(row[key]) is float:
                    row[key] = round(row[key],2)
                elif row[key] is None:
                    row[key]="NA"
            rower.writerow(row)
        return
    
    def get_metrics(self):
        
        '''Quick wrapper to calculated the metrics for all the transcripts.'''
        
        if self.metrics_calculated is True:
            return
        
        self.metrics=dict()
        
        for tid in self.transcripts:
            self.calculate_metrics(tid)
            self.metrics[tid]=copy(self.transcripts[tid].metrics)

        self.metrics_calculated = True
        return


    ############### Class methods ################

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
    
    ########### Properties
    
    @property
    def splitted(self):
        '''The splitted flag indicates whether a sublocus has already been processed to produce the necessary monosubloci.
        It must be set as a boolean flag (hence why it is coded as a property)'''
        return self.__splitted
    
    @splitted.setter
    def splitted(self,verified):
        if type(verified)!=bool:
            raise TypeError()
        self.__splitted=verified
        
    @property
    def id(self):
        if self.monoexonic is True:
            addendum = "mono"
        else:
            addendum = "multi"
        
        return "{0}.{1}".format(super().id, addendum)