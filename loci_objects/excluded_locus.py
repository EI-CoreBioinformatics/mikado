import sys,os.path
from loci_objects.get_metrics_name import get_metrics_name
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from loci_objects.abstractlocus import abstractlocus

#Resolution order is important here!
class excluded_locus(abstractlocus):
    
    '''This is a container that groups together the transcripts surviving the selection for the monosublocus.
    The class inherits from both sublocus and abstractlocus (the main abstract class) in order to be able to reuse
    some of the code present in the former.
    Internally, the most important method is define_loci - which will select the best transcript(s) and remove all the overlapping ones.
    The intersection function for this object is quite laxer than in previous stages, and so are the requirements for the inclusion.
    '''
    @staticmethod
    def get_available_metrics(filename=None):
        '''Wrapper for the "get_metrics_name" function to retrieve the necessary metrics names.'''
        
        return get_metrics_name(filename=filename)
    
    __name__ = "excluded_transcripts"
    available_metrics = []
    if available_metrics == []:
        available_metrics = get_available_metrics.__func__()

    def __init__(self, monosublocus_instance, json_dict=None, metrics=None):
        
        abstractlocus.__init__(self)
        self.splitted=False
        self.metrics_calculated = False
        self.json_dict = json_dict
        #Add the transcript to the locus
        self.add_monosublocus(monosublocus_instance)

    def add_transcript_to_locus(self, transcript_instance):
        '''Override of the sublocus method, and reversal to the original method in the abstractlocus class.'''
        abstractlocus.add_transcript_to_locus(self, transcript_instance, check_in_locus=False)
            
    def add_monosublocus(self, monosublocus_instance):
        '''Wrapper to extract the transcript from the monosubloci and pass it to the constructor.'''
        assert len(monosublocus_instance.transcripts)==1
        for tid in monosublocus_instance.transcripts:
            self.add_transcript_to_locus(monosublocus_instance.transcripts[tid])
            
    def __str__(self):
        '''This special method is explicitly *not* implemented; this locus object is not meant for printing, only for computation!'''
        raise NotImplementedError("This is a container used for computational purposes only, it should not be printed out directly!")

    def print_metrics(self):
        
        '''This class yields dictionary "rows" that will be given to a csv.DictWriter class.'''
        
        #Check that rower is an instance of the csv.DictWriter class
#        self.get_metrics()
#         self.calculate_scores() 
        
        #The rower is an instance of the DictWriter class from the standard CSV module
        
        for tid in sorted(self.transcripts.keys(), key=lambda tid: self.transcripts[tid] ):
            row={}
            for key in self.available_metrics:
                if key.lower() in ("id", "tid"):
                    row[key]=tid
                elif key.lower()=="parent":
                    row[key]=self.id
                else:
                    row[key]=getattr(self.transcripts[tid], key, "NA")
                if type(row[key]) is float:
                    row[key] = round(row[key],2)
                elif row[key] is None or row[key]=="":
                    row[key]="NA"
            yield row
        return


    @classmethod
    def is_intersecting(cls):
        raise NotImplementedError()