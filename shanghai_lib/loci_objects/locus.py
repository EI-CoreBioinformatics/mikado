import sys,os.path
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from shanghai_lib.loci_objects.monosublocus import monosublocus
from shanghai_lib.loci_objects.abstractlocus import abstractlocus

class locus(monosublocus,abstractlocus):
    
    '''Minimal which inherits from monosublocus. It will define the final loci.''' 
    
    __name__ = "locus"
    
    def __init__(self,transcript_instance):
        self.counter=0
        super().__init__(transcript_instance)
    
    def __str__(self, print_cds=True):
          
        self.feature=self.__name__
        return super().__str__(print_cds=print_cds)
    
    def other_is_fragment(self,other):
        '''This function checks whether another *monoexonic* locus on the opposite strand* is a fragment, by checking that:
            - the fragment is completely contained inside the coordinates of the other transcript
            - the transcript has at least some exonic overlap and overlaps at most one exon.
        '''
        
        if self.monoexonic is True or other.monoexonic is False or other.strand==self.strand:
            return False
        if type(self)!=type(other):
            raise TypeError("I can compare only loci.")

        other_transcript=other.transcripts[list(other.transcripts.keys())[0]]
        if self.overlap( (other.start, other.end), (self.start,self.end))<other_transcript.cdna_length-1:
            return False
        overlapping_exons=0
        for exon in self.exons:
            if self.overlap(exon, (other.start,other.end))>0:
                overlapping_exons+=1
        if overlapping_exons!=1: return False
        return True
    
    @property
    def id(self):
        Id = abstractlocus.id.fget(self)  # @UndefinedVariable
        if self.counter>0:
            Id = "{0}.{1}".format(Id, self.counter)
        return Id
        