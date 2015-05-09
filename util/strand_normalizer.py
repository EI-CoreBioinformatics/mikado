import sys,os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from loci_objects.GFF import GFF3
from loci_objects.GTF import GTF
from loci_objects.transcript import transcript
from Bio import SeqIO
import argparse

class transcript_checker(transcript):
    
    def __init__(self, gffLine, fasta_index, strand_specific=False):
        if fasta_index is None:
            raise ValueError()
        super().__init__(gffLine)
        self.parent = gffLine.parent 
        self.fasta_index = fasta_index
        self.strand_specific=strand_specific
                
    @property
    def strand_specific(self):
        return self.__strand_specific
    
    @strand_specific.setter
    def strand_specific(self,*args):
        if type(args[0]) is not bool:
            raise TypeError("Invalid value for boolean property: {0}".format(args[0]))
        self.__strand_specific=args[0]
                
    def __str__(self):
        
        self.check_strand()
        return super().__str__()
    
    def check_strand(self):
        
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
            
            for intron in self.introns:
                splice_donor = self.fasta_index[self.chrom][intron[0]-1:intron[0]+1]
                assert len(splice_donor)==2
                splice_acceptor = self.fasta_index[self.chrom][intron[0]-2:intron[0]]
                assert len(splice_acceptor)==2
                if self.strand == "-":
                    splice_acceptor = splice_donor.reverse_complement()
                    splice_donor = splice_acceptor.reverse_complement()
                    
                if (splice_donor,splice_acceptor) in canonical_splices:
                    canonical_counter["+"]+=1
                else:
                    reversed_splice_donor = splice_acceptor.reverse_complement()
                    reversed_splice_acceptor = splice_donor.reverse_complement()
                    if (reversed_splice_donor,reversed_splice_acceptor) in canonical_splices:
                        canonical_counter["-"]+=1
                    else:
                        canonical_counter[None]+=1

            if canonical_counter["+"]>canonical_counter["-"]+canonical_counter[None]:
                pass
            elif canonical_counter["-"]==len(self.introns):
                self.reverse_strand()
            elif len(canonical_counter[None])==len(self.introns):
                self.strand=None

def main():
    
    def to_gff(string):
        if not os.path.exists(string) or not os.path.isfile(string) or not os.stat(string).st_size>0:
            raise ValueError("Invalid input file.")
        if  string[-4:]==".gtf":
            gff_function=GTF
        else:
            gff_function=GFF3 
        
        gff=gff_function(string)
        record=next(gff)
        for record in gff:
            if record.header is False:
                gff.close()
                return gff_function(string)
        raise ValueError("Empty GFF file provided.")
    
    def to_seqio(string):
        if not os.path.exists(string) or not os.path.isfile(string) or not os.stat(string).st_size>0:
            raise ValueError("Invalid input file.")
        return SeqIO.index(string, "fasta")
    
    parser=argparse.ArgumentParser("Script to check the correctness of the strand for aligned/assembled transcripts.")
    parser.add_argument("--fasta", type=to_seqio, required=True,
                        help="Genome FASTA file. Required." )
    parser.add_argument("gff", type=to_gff, help="Input GFF3/GTF file.")
    parser.add_argument("out", default=sys.stdout, nargs='?', type=argparse.FileType('w'),
                        help="Output file. Default: STDOUT.")
    args=parser.parse_args()
    
    currentTranscript = None
    
    if args.gff.name[-4:]=="gtf":
        is_gff=False
    else:
        is_gff=True
        
    for record in args.gff:
        if record.is_parent:
            continue
        elif record.header:
            print(record, file=args.out)
        elif record.is_transcript:
            if currentTranscript is not None:
                print(currentTranscript, file=args.out)
            currentTranscript = transcript_checker(record, args.fasta)
        elif record.is_exon:
            currentTranscript.addExon(record)
            
    if currentTranscript is not None:
        print(currentTranscript, file=args.out)
            
            
    
if __name__=='__main__': main()




                    
                    