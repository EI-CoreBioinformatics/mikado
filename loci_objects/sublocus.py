from loci_objects.abstractlocus import abstractlocus
#import random
#from copy import copy
import re
from loci_objects.monosublocus import monosublocus
from loci_objects.transcript import transcript

class sublocus(abstractlocus):
    
    '''
    The sublocus class is created either by the superlocus class during the subloci definition, or directly using a G(T|F)line-like object.
    It is used to define the final monosubloci.
    '''
    
    __name__ = "sublocus"
    available_metrics = [
                "tid",
                "parent",
                "score",
                "exon_num",
                "exon_fraction",
                "intron_fraction",
                "retained_intron_num",
                "retained_fraction",
                "combined_cds_length",
                "combined_cds_num",
                "combined_cds_num_fraction",
                "combined_cds_fraction",
                "combined_utr_length",
                "cdna_length",
                "cds_length",
                "utr_num",
                "five_utr_length",
                "five_utr_num",
                "three_utr_length",
                "three_utr_num",
                "number_internal_orfs",
                "cds_fraction",
                "best_cds_number",
                "best_cds_number_fraction",
                "cds_num",
                "cds_fraction",
                "cds_not_maximal",
                "cds_not_maximal_fraction", 
                "has_start",
                "has_stop",
                "is_complete",
                'utr_fraction',
                ]
    
    am=available_metrics[:3]
    am.extend(sorted(available_metrics[3:]))
    available_metrics=am[:]
    del am
    
    ################ Class special methods ##############
    
    def __init__(self, span, source=None, json_dict = None ):
        
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
        setattr( self, "monoexonic", getattr(span, "monoexonic", None)  )
        if json_dict is None or type(json_dict) is not dict:
            raise ValueError("I am missing the configuration for prioritizing transcripts!")
        self.json_dict = json_dict
        
        #This part is necessary to import modules
        if "modules" in self.json_dict:
            import importlib
            for mod in self.json_dict["modules"]:
                globals()[mod]=importlib.import_module(mod)
        
        assert hasattr(self, "monoexonic")

        if type(span) is transcript:
            self.add_transcript_to_locus(span)
        else:
            for key in ["parent", "start", "end", "chrom", "strand", "attributes"]:
                setattr(self, key, getattr(span, key))
        
        
    def __str__(self):
        if self.strand is not None:
            strand=self.strand
        else:
            strand="."
        
        lines=[]
        
        if self.splitted is False:
            attr_field="ID={0};Name={1};Parent={2};multiexonic={3}".format(
                                                              self.id,
                                                              self.name,
                                                              self.parent,
                                                              str(not self.monoexonic)
                                                            )
                            
            self_line = [ self.chrom, self.source, "sublocus", self.start, self.end,
                     ".", strand, ".", attr_field]
            lines.append("\t".join([str(s) for s in self_line]))
        
            for tid in sorted(self.transcripts, key=lambda tid: self.transcripts[tid]):
                lines.append(str(self.transcripts[tid]).rstrip())
        else:
            for slocus in sorted(self.monosubloci): 
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
        '''This function will calculate the metrics for a transcript which are relative in nature
        i.e. that depend on the other transcripts in the sublocus. Examples include the fraction
        of introns or exons in the sublocus, or the number/fraction of retained introns.  
        '''
    
        transcript_instance = self.transcripts[tid]
        self.transcripts[tid].finalize() # The transcript must be finalized before we can calculate the score.
        
        
        self.transcripts[tid].exon_fraction = len(set.intersection( self.exons,self.transcripts[tid].exons   ))/len(self.exons)
        if len(self.junctions)>0:
            transcript_instance.intron_fraction = len(set.intersection( self.junctions,transcript_instance.junctions   ))/len(self.junctions)
        else:
            transcript_instance.intron_fraction = 0
        self.find_retained_introns(transcript_instance)
        transcript_instance.retained_fraction=sum(e[1]-e[0]+1 for e in transcript_instance.retained_introns)/transcript_instance.cdna_length
        self.transcripts[tid]=transcript_instance
        
    def find_retained_introns(self, transcript_instance):
         
        '''This method checks the number of exons that are possibly retained introns for a given transcript.
        To perform this operation, it checks for each non-CDS exon whether it exists a sublocus intron that
        is *completely* contained within a transcript exon.
        CDS exons are ignored because their retention might be perfectly valid.
        The results are stored inside the transcript instance, in the "retained_introns" tuple.'''
         
        transcript_instance.retained_introns=[]
        
        for exon in filter(lambda e: e not in transcript_instance.non_overlapping_cds, transcript_instance.exons):
            #Check that the overlap is at least as long as the minimum between the exon and the intron.
            if any(filter(
                          lambda junction: self.overlap(exon,junction)>=junction[1]-junction[0],
                          self.junctions                          
                          )) is True:
                    transcript_instance.retained_introns.append(exon)
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
            
            
    def calculate_scores(self):
        '''
        Function to calculate a score for each transcript, given the metrics derived
        with the calculate_metrics method.
        The exact ordering must be provided from outside using a JSON file.
         
         '''
        
        self.get_metrics()
        self.metrics=dict()
        
        for tid in self.transcripts:
            values = []
            for param in self.json_dict["parameters"]:
                val = getattr(self.transcripts[tid], param)
                if "operation" in self.json_dict["parameters"][param]:
                    try:
                        val = eval( re.sub("param", str(val),  self.json_dict["parameters"][param]["operation"] ) )
                    except ZeroDivisionError:
                        val = 0
                if "multiplier" in self.json_dict["parameters"][param]:
                    val*=self.json_dict["parameters"][param]["multiplier"]
                elif "bool_value" in self.json_dict["parameters"][param]:
                    assert type(val) is bool
                    if val is True:
                        val = self.json_dict["parameters"][param]["bool_value"][0]
                    else:
                        val = self.json_dict["parameters"][param]["bool_value"][1]
                assert type(val) in (int,float), (param, val)
                values.append(val)
            score=sum(values)
            self.metrics[tid]=dict()
            self.transcripts[tid].score=self.metrics[tid]["score"]=score
        
        if "requirements" in self.json_dict:
            for key in self.json_dict["requirements"]:
                key, conf = key, self.json_dict["requirements"][key]
                if conf["type"]=="min":
                    oper=">="
                elif conf["type"]=="eq":
                    oper="=="
                elif conf["type"]=="max":
                    oper="<="
                elif conf["type"]=="uneq":
                    oper="!="
                else:
                    raise TypeError("Cannot recognize this type: {0}".format(conf["type"]))
                validator_str = "{score} {oper} {ref_val}"
                if hasattr(self.transcripts[tid], key):
                    val=getattr(self.transcripts[tid], key)
                else:
                    val=self.transcripts[tid].metrics[key]
                not_passing = set(filter( lambda tid: eval(validator_str.format(
                                                                              score=val,
                                                                              oper=oper,
                                                                              ref_val=conf["value"]
                                                                              )) is False, self.metrics
                                           ))
                if len(not_passing)==len(self.metrics): #all transcripts in the locus fail to pass the filter
                    continue
                else:
                    for tid in not_passing:
                        self.transcripts[tid].score=self.metrics[tid]["score"]=0
                      
    
    def print_metrics(self):
        
        '''This class yields dictionary "rows" that will be given to a csv.DictWriter class.'''
        
        #Check that rower is an instance of the csv.DictWriter class
#        self.get_metrics()
        self.calculate_scores() 
        
        #The rower is an instance of the DictWriter class from the standard CSV module
        
        for tid in sorted(self.transcripts.keys(), key=lambda tid: self.transcripts[tid] ):
            row={}
            for key in self.available_metrics:
                if key.lower() in ("id", "tid"):
                    row[key]=tid
                elif key.lower()=="parent":
                    row[key]=self.id
                elif key=="score":
                    row[key]=self.transcripts[tid].score
                else:
                    row[key]=getattr(self.transcripts[tid], key, "NA")
                if type(row[key]) is float:
                    row[key] = round(row[key],2)
                elif row[key] is None or row[key]=="":
                    row[key]="NA"
            yield row
        return
    
    def get_metrics(self):
        
        '''Quick wrapper to calculate the metrics for all the transcripts.'''
        
        if self.metrics_calculated is True:
            return
        
        for tid in self.transcripts:
            self.calculate_metrics(tid)

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