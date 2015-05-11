import sys,os
from copy import copy
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from loci_objects.GFF import GFF3
from loci_objects.GTF import GTF
from loci_objects.transcript_checker import transcript_checker
from Bio import SeqIO
import argparse

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
    parser.add_argument("-s", "--strand-specific", dest="strand_specific", action="store_true", default=False,
                        help="Flag. If set, monoexonic transcripts will be left on their strand rather than being moved to the unknown strand.")
    parser.add_argument("-l", "--lenient", action="store_true", default=False,
                        help="Flag. If set, transcripts with mixed +/- splices will not cause exceptions but rather be annotated as problematic.")
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
            new_tran = transcript_checker(record, currentSeq, lenient=args.lenient, strand_specific=args.strand_specific)
            currentTranscripts.append(new_tran)
        elif record.is_exon:
            for tran in currentTranscripts:
                if tran.id not in record.parent:
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




                    
                    