import sys,os
from copy import copy
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
        self.original_strand = gffLine.strand
        assert self.original_strand==self.strand
        self.parent = gffLine.parent 
        self.fasta_index = fasta_index
        self.strand_specific=strand_specific
        self.checked = False
                
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

            if canonical_counter["+"]>canonical_counter["-"]+canonical_counter[None]:
                pass
            elif canonical_counter["-"]==len(self.introns):
                self.reverse_strand()
            elif canonical_counter[None]==len(self.introns):
                self.strand=None

        self.checked = True

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
        print("Loading reference")
#         seqdict = SeqIO.index(string, "fasta")
        seqdict = SeqIO.to_dict(SeqIO.parse(open(string), 'fasta'))
        print("Reference loaded") 
        return seqdict
    
    parser=argparse.ArgumentParser("Script to check the correctness of the strand for aligned/assembled transcripts.")
    parser.add_argument("--fasta", type=to_seqio, required=True,
                        help="Genome FASTA file. Required." )
    parser.add_argument("gff", type=to_gff, help="Input GFF3/GTF file.")
    parser.add_argument("out", default=sys.stdout, nargs='?', type=argparse.FileType('w'),
                        help="Output file. Default: STDOUT.")
    args=parser.parse_args()
    
    currentTranscripts = []
    currentParent = None
    
    if args.gff.name[-4:]=="gtf":
        is_gff=False
    else:
        is_gff=True
    
    
    currentSeq = None
    
    for record in args.gff:
        if record.header is False:
            assert record.chrom in args.fasta
            sequence = args.fasta[record.chrom]
            if currentSeq is None or list(currentSeq.keys())[0]!=record.chrom:
                print("Loading sequence {0}".format(record.chrom))
                currentSeq=dict()
                currentSeq[record.chrom] = sequence.seq
                print("Loaded sequence {0}".format(record.chrom))
        
        if record.is_parent:
            if currentParent is not None:
                if len(currentTranscripts)==0:
                    print(currentParent, file=args.out)
                else:
                    for tran in currentTranscripts:
                        tran.check_strand()
                    strands = set([t.strand for t in currentTranscripts])
                    if len(strands)==1:
                        strand = strands.pop()
                        currentParent.strand = strand
                        print(currentParent, file=args.out)
                        for tran in currentTranscripts:
                            print(tran, file=args.out)
                    else:
                        original_strand = currentParent.strand
                        print(currentParent, file=args.out)
                        for tran in filter(lambda t: t==original_strand, currentTranscripts):
                            print(tran, file=args.out)
                        for strand in filter(lambda s: s!=original_strand, strands):
                            newPar = copy(currentParent)
                            new_id = "{0}:{1}".format(currentParent.id, strand)
                            newPar.id = new_id
                            newPar.strand = strand
                            print(newPar, file=args.out)
                            for tran in filter(lambda t: t==original_strand, currentTranscripts):
                                tran.parent = new_id
                                print(tran, file=args.out)
            currentParent=record
            currentTranscripts=[]
        elif record.header:
            print(record, file=args.out)
        elif record.is_transcript:
            if is_gff is False:
                print(currentTranscripts[0], file=args.out)
            new_tran = transcript_checker(record, currentSeq)
            currentTranscripts.append(new_tran)
        elif record.is_exon:
            for tran in currentTranscripts:
                if tran.id==record.parent:
                    tran.addExon(record)


    if is_gff is False and len(currentTranscripts)>0:
        print(currentTranscripts[0], file=args.out)
    elif is_gff is True:
        if currentParent is not None:
            if len(currentTranscripts)==0:
                print(currentParent, file=args.out)
            else:
                for tran in currentTranscripts:
                    tran.check_strand()
                strands = set([t.strand for t in currentTranscripts])
                if len(strands)==1:
                    currentParent.strand = strands.pop()
                    print(currentParent, file=args.out)
                    for tran in currentTranscripts:
                        print(tran, file=args.out)
                else:
                    original_strand = currentParent.strand
                    print(currentParent, file=args.out)
                    for tran in filter(lambda t: t==original_strand, currentTranscripts):
                        print(tran, file=args.out)
                    for strand in filter(lambda s: s!=original_strand, strands):
                        newPar = copy(currentParent)
                        new_id = "{0}:{1}".format(currentParent.id, strand)
                        newPar.id = new_id
                        print(newPar, file=args.out)
                        for tran in filter(lambda t: t==original_strand, currentTranscripts):
                            tran.parent = new_id
                            print(tran, file=args.out)
            
            
    
if __name__=='__main__': main()




                    
                    