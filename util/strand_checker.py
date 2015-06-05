import sys,os
from copy import deepcopy
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import loci_objects
from Bio import SeqIO
import argparse

def main():
    
    def to_gff(string):
        if not os.path.exists(string) or not os.path.isfile(string) or not os.stat(string).st_size>0:
            raise ValueError("Invalid input file.")
        if  string[-4:]==".gtf":
            gff_function=loci_objects.GTF.GTF
        else:
            gff_function=loci_objects.GFF.GFF3 
        
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
                        help="""Flag. If set, transcripts with mixed +/- splices will not cause exceptions but rather be annotated as problematic.
If set, the output will be GFF3, regardless of the input format.""")
    parser.add_argument("gff", type=to_gff, help="Input GFF3/GTF file.")
    parser.add_argument("out", default=sys.stdout, nargs='?', type=argparse.FileType('w'),
                        help="Output file. Default: STDOUT.")
    args=parser.parse_args()
    
    currentTranscripts = dict()
    currentParent = None
    
    if args.gff.name[-3:]=="gtf":
        is_gff=False
    else:
        is_gff=True

#   print(is_gff, file=sys.stderr)
    
    currentSeq = None
    
    reversed_transcripts=0
    
    for record in args.gff:
        
        if currentSeq is None or record.chrom not in currentSeq:
            currentSeq=dict()
            currentSeq[record.chrom]=args.fasta[record.chrom]
        if record.header is True:
            continue
        if record.is_parent is True:
            if is_gff is True and record.is_transcript is False:
                to_delete=[]
                for tid,tr in currentTranscripts.items():
                        try:
                            tr.check_strand()
                            if tr.reversed is True:
                                reversed_transcripts+=1
                        except loci_objects.exceptions.IncorrectStrandError:
                            to_delete.append(tid)
                        except loci_objects.exceptions.InvalidTranscript:
                            to_delete.append(tid)
                for tid in to_delete:
                    del currentTranscripts[tid]
                if len(currentTranscripts)>0:
                    strands = dict()
                    for tid, tr in currentTranscripts.items():
                        if tr.reversed is True: reversed_transcripts+=1
                        if tr.strand not in strands:
                            strands[tr.strand]=[]
                        strands[tr.strand].append(tr)
                    for strand in strands:
                        if strand == currentParent.strand:
                            print(currentParent, file=args.out)
                            parent_id=currentParent.id
                        else:
                            newParent = deepcopy(currentParent)
                            newParent.strand = strand
                            newParent.id = "{0}.strand{1}".format(currentParent.id, strand)
                            print(newParent, file=args.out)
                            parent_id = newParent.id
                        for tr in strands[strand]:
                            tr.parent = parent_id
                            print(tr, file=args.out)
                currentParent=record
                currentTranscripts=dict()

            else:
                tr = currentTranscripts[currentTranscripts.keys()[0]]
                try:
                    tr.check_strand()
                    print(tr, file=args.out)
                except loci_objects.exceptions.IncorrectStrandError:
                    pass
                currentTranscripts=dict()
                currentTranscript = loci_objects.transcript_checker.transcript_checker(
                                                                                       record,
                                                                                       currentSeq,
                                                                                       strand_specific=args.strand_specific,
                                                                                       lenient=args.lenient
                                                                                       )

                currentTranscripts[currentTranscript.id]=currentTranscript

        elif record.is_transcript is True and is_gff is True:
            try:
                assert currentParent.id in record.parent, record
            except TypeError:
                raise TypeError(str(currentParent), str(record))
            currentTranscript = loci_objects.transcript_checker.transcript_checker(
                                                                               record,
                                                                               currentSeq,
                                                                               strand_specific=args.strand_specific,
                                                                               lenient=args.lenient
                                                                            )            
            currentTranscripts[currentTranscript.id]=currentTranscript

        elif record.is_exon is True:
            assert len(currentTranscripts)>0, (str(currentParent),currentTranscripts)
            for parent in record.parent:
                try:
                    currentTranscripts[parent].addExon(record)
                except KeyError:
                    print(currentTranscripts)
                    raise

    to_delete = []         
    for tid,tr in currentTranscripts.items():
        try:
            tr.check_strand()
            if tr.reversed is True:
                reversed_transcripts+=1
        except loci_objects.exceptions.IncorrectStrandError:
            to_delete.append(tid)
    for tid in to_delete:
        del currentTranscripts[tid]
    if len(currentTranscripts)>0:
        strands = dict()
        for tid, tr in currentTranscripts.items():
            if tr.strand not in strands:
                strands[tr.strand]=[]
            strands[tr.strand].append(tr)
        for strand in strands:
            parent_id = None
            if is_gff is True:
                if strand == currentParent.strand:
                    print(currentParent, file=args.out)
                    parent_id=currentParent.id
                else:   
                    newParent = deepcopy(currentParent)
                    newParent.strand = strand
                    newParent.id = "{0}.strand{1}".format(currentParent.id, strand)
                    print(newParent, file=args.out)
                    parent_id = newParent.id
            for tr in strands[strand]:
                if parent_id is not None: tr.parent = parent_id
                print(tr, file=args.out)

if __name__=='__main__': main()




                    
                    
