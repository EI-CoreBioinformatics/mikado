from loci_objects.abstractlocus import abstractlocus
from loci_objects.sublocus import sublocus
from loci_objects.locus import locus

#Resolution order is important here!
class monosublocus_holder(sublocus,abstractlocus):
    
    '''This is a container that groups together the transcripts surviving the selection for the monosublocus.
    The class inherits from both sublocus and abstractlocus (the main abstract class) in order to be able to reuse
    some of the code present in the former.
    Internally, the most important method is define_loci - which will select the best transcript(s) and remove all the overlapping ones.
    The intersection function for this object is quite laxer than in previous stages, and so are the requirements for the inclusion.
    '''

    __name__ = "monosubloci_holder"

    def __init__(self, monosublocus_instance, json_dict=None):
        
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
        
    def calculate_metrics(self, tid):
        '''This function will recalculate the *relative* metrics for a transcript,
        which have varied since now we are considering only a subset of the original transcripts.
        Relative metrics include, at this stage, only the fraction of retained introns.
        '''
    
        if self.metrics_calculated is True:
            return
        
        transcript_instance = self.transcripts[tid]
        #Check that metrics had already been calculated
        assert transcript_instance.finalized is True

        self.find_retained_introns(transcript_instance)
        transcript_instance.parent=self.id
        self.transcripts[tid]=transcript_instance

        self.metrics_calculated = True
        return

    def define_monosubloci(self):
        '''Overriden and set to NotImplemented to avoid cross-calling it when inappropriate.'''
        raise NotImplementedError("Monosubloci are the input of this object, not the output.")
    
    def define_loci(self, purge=False):
        '''This is the main function of the class. It is analogous to the define_subloci class defined
        for sublocus objects, but it returns "locus" objects (not "monosublocus").'''
        if self.splitted is True:
            return
        self.calculate_scores()
        
        self.loci=[]
        remaining = self.transcripts.copy()
        
        while len(remaining)>0:
            best_tid=self.choose_best(remaining.copy())
            best_transcript = remaining[best_tid]
            new_remaining = remaining.copy()
            del new_remaining[best_tid]
            if best_transcript.score==0 and purge is True:
                pass
            else:
                new_locus = locus(best_transcript)
                self.loci.append(new_locus)
            for tid in remaining:
                if tid==best_tid: continue
                if self.is_intersecting(best_transcript, new_remaining[tid]):
                    del new_remaining[tid]
            remaining=new_remaining.copy()
    
        self.splitted = True
        return


    @classmethod
    def is_intersecting(cls, transcript_instance, other):
        '''
        Implementation of the is_intersecting method. Now that we are comparing transcripts that
        by definition span multiple subloci, we have to be less strict in our definition of what
        counts as an intersection.
        Criteria:
        - 1 splice site in common (splice, not junction)
        - 2+ overlapping exons or 1 exon for monoexonic
        
        The overlap check must be done both ways, because doing it only in one sense might miss cases like this:
        
        A ++++++++++++++++++++++--------+++++
        B ++++------+++++---------++++
        
        If we do it only one way, is_intersecting(A,B)!=is_intersecting(B,A), which is
        NOT a desirable outcome (the result is that they are indeed intersecting).
        '''
        if transcript_instance==other:
            return False # We do not want intersection with oneself

        if cls.overlap((transcript_instance.start,transcript_instance.end), (other.start,other.end) )<=0: return False
        if len(set.intersection( set(transcript_instance.junctions), set(other.junctions)))>0:
            return True
        
        overlapping_exons = set()
        other_overlapping_exons = set()
        
        
        for exon in transcript_instance.exons:
            for oexon in other.exons:
                
                if cls.overlap(exon, oexon) > 0:
                    other_overlapping_exons.add(oexon)
                    overlapping_exons.add(exon)
                if max(len(other_overlapping_exons), len(overlapping_exons))>1 or \
                    (len(other_overlapping_exons)==1 and other.monoexonic is True) or \
                    (len(overlapping_exons)==1 and transcript_instance.monoexonic is True):
                    return True
                    
        
        return False

    @property
    def id(self):
        return abstractlocus.id.fget(self)  # @UndefinedVariable