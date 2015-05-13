import sys,os.path
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from loci_objects.monosublocus import monosublocus
from loci_objects.abstractlocus import abstractlocus

class locus(monosublocus,abstractlocus):
    
    '''Minimal which inherits from monosublocus. It will define the final loci.''' 
    
    __name__ = "locus"
    
    def __init__(self,transcript_instance):
        self.counter=0
        super().__init__(transcript_instance)
    
    def __str__(self, print_cds=True):
          
        self.feature=self.__name__
        return super().__str__(print_cds=print_cds)
    
    def other_is_fragment(self,other, percentage=1):
        '''This function checks that another *monoexonic* locus on the opposite strand* does not verify one of the following:
            - it is contained for more than (exon_length)*percentage inside one of the locus exons
            - it is not partially contained inside an intron
        This should get rid of monoexonic fragments that plague the output of RNA-Seq reconstruction programs.
        '''
        
        if type(percentage) not in (float,int) or not 0<percentage<=1:
            raise ValueError("Invalid percentage, it should be between 0 and 1. Received: {0}".format(percentage))
        
        if self.monoexonic is True or other.monoexonic is False or other.strand==self.strand:
            return False
        if type(self)!=type(other):
            raise TypeError("I can compare only loci.")

        other_exon = list(other.exons)[0]
        threshold = (other_exon[1]+1-other_exon[0])*percentage
        for exon in self.exons:
            if self.overlap( other_exon, exon  )>=threshold:
                return True
            
        for intron in self.introns:
            if self.overlap(other_exon,intron)>=threshold: 
                return True
            
        return False
    
    @property
    def id(self):
        Id = abstractlocus.id.fget(self)  # @UndefinedVariable
        if self.counter>0:
            Id = "{0}.{1}".format(Id, self.counter)
        return Id
        