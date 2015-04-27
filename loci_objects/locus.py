from loci_objects.monosublocus import monosublocus
from loci_objects.abstractlocus import abstractlocus

class locus(monosublocus,abstractlocus):
    
    '''Minimal which inherits from monosublocus. It will define the final loci.''' 
    
    __name__ = "locus"
    
    def __str__(self):
        
        lines=[]
        if self.strand is not None:
            strand=self.strand
        else:
            strand="."
            
        attr_field="ID={0};{1}".format(
                                        self.id,
                                        "{0};".format(self.parent) if self.parent is not None else "",
                                        )
        self_line = [ self.chrom, self.source, self.__name__, self.start, self.end,
                             ".", strand, ".", attr_field]
        lines.append("\t".join([str(s) for s in self_line]))
        for tid in self.transcripts:
            transcript_instance=self.transcripts[tid]
            lines.append(str(transcript_instance).rstrip())
        return "\n".join(lines)
    
    
    @property
    def id(self):
        return abstractlocus.id.fget(self)  # @UndefinedVariable