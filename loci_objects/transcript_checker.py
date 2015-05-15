import sys,os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from loci_objects.transcript import transcript


class IncorrectStrandError(Exception):
    pass


class transcript_checker(transcript):
    
    '''This is a subclass of the generic transcript class. Its purpose is to compare
    the information of the transcript instance with the information contained in a
    genomic FASTA file, to verify some of the information.
    At the moment, the class implements only a check on the strandedness made by extracting
    the FASTA sequence of the splice sites and verifying that they are concordant with the annotated strand.
    Keyword arguments:
        - strand_specific: if set, monoexonic transcripts are not set to "unknown" strand.
        - lenient: boolean. If set to True, a transcript with mixed splices will not throw an exception
        but rather will just report the number of splices supporting each strand.
    '''
    
    def __init__(self, gffLine, fasta_index, strand_specific=False, lenient=False):
        if fasta_index is None:
            raise ValueError()
        super().__init__(gffLine)
        self.original_strand = gffLine.strand
        assert self.original_strand==self.strand
        self.parent = gffLine.parent 
        self.fasta_index = fasta_index
        self.strand_specific=strand_specific
        self.checked = False
        self.lenient=lenient
        self.mixed_splices = False
                
    @property
    def strand_specific(self):
        return self.__strand_specific
    
    @strand_specific.setter
    def strand_specific(self,*args):
        if type(args[0]) is not bool:
            raise TypeError("Invalid value for boolean property: {0}".format(args[0]))
        self.__strand_specific=args[0]
                
    def __str__(self, print_cds=True):
        
        self.check_strand()
        if self.mixed_splices is True:
            self.attributes["mixed_splices"]=self.mixed_attribute
        
        return super().__str__(print_cds=print_cds)
    
    def check_strand(self):
        
        self.finalize()
        if self.checked is True:
            return
        
        canonical_splices = [
                             ("GT","AG"),
                             ("GC","AG"),
                             ("AT","AC") 
                             ]
        
        if self.strand_specific is False and self.monoexonic is True:
            self.strand=None
            return
        
        elif self.monoexonic is False:
            canonical_counter=dict()
            for strand in ("+","-",None):
                canonical_counter[strand]=0
            
            assert len(self.introns)>0
            
            
            for intron in self.introns:
                splice_donor = self.fasta_index[self.chrom][intron[0]-1:intron[0]+1]
                splice_acceptor = self.fasta_index[self.chrom][intron[1]-2:intron[1]]
                if (str(splice_donor),str(splice_acceptor)) in canonical_splices:
                    if self.strand == "+":
                        canonical_counter["+"]+=1
                    elif self.strand == "-":
                        canonical_counter["-"]+=1
                else:
                    rsa = splice_donor.reverse_complement()
                    rsd = splice_acceptor.reverse_complement()
                    splice_acceptor, splice_donor = rsa, rsd
                    if (str(splice_donor),str(splice_acceptor)) in canonical_splices:
                        if self.strand=="-":
                            canonical_counter["+"]+=1
                        else:
                            canonical_counter["-"]+=1
                    else:
                        canonical_counter[None]+=1

            if canonical_counter["+"]>0 and canonical_counter["-"]>0:
                if self.lenient is False:
                    raise IncorrectStrandError("Transcript {0} has {1} positive and {2} negative splice junctions. Aborting.".format(
                                                                                                                                 self.id,
                                                                                                                                 canonical_counter["+"],
                                                                                                                                 canonical_counter["-"]
                                                                                                                                 )
                                           )
                else:
                    self.mixed_splices = True

                if canonical_counter["+"]>=canonical_counter["-"]:
                    self.mixed_attribute="{0}concordant,{1}discordant".format(canonical_counter["+"],canonical_counter["-"])
                else:
                    self.reverse_strand()
                    self.mixed_attribute="{0}concordant,{1}discordant".format(canonical_counter["-"],canonical_counter["+"])
                return
            
            if canonical_counter["-"]>0:
                self.reverse_strand()

        self.checked = True
