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
          
    @property
    def id(self):
        return abstractlocus.id.fget(self)  # @UndefinedVariable