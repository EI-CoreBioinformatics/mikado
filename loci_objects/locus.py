import sys,os.path
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from loci_objects.monosublocus import monosublocus
from loci_objects.abstractlocus import abstractlocus

class locus(monosublocus,abstractlocus):
    
    '''Minimal which inherits from monosublocus. It will define the final loci.''' 
    
    __name__ = "locus"
    
    def __str__(self):
          
        self.feature=self.__name__
        return super().__str__()
    
    def other_is_fragment(self,other, percentage=1):
        '''This function checks that another monoexonic locus is not completely contained within an exon of the instance,
        only on the opposite strand.
        This should get rid of monoexonic fragments that plague the output of RNA-Seq reconstruction programs.'''
        
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
        return False
    
    @property
    def id(self):
        return abstractlocus.id.fget(self)  # @UndefinedVariable