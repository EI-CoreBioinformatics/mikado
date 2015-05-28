import sys,os,argparse
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from loci_objects import exceptions
from loci_objects.transcript_checker import transcript_checker
from loci_objects import GTF
from Bio import SeqIO

def main():
    
    def to_gff(string):
        if  string[-4:]!=".gtf":
            raise TypeError("This script takes as input only GTF files.")
        
        gff=GTF.GTF(string)
        record=next(gff)
        for record in gff:
            if record.header is False:
                gff.close()
                return GTF.GTF(string)
        raise ValueError("Empty GTF file provided.")
    
    def to_seqio(string):
        if not os.path.exists(string) or not os.path.isfile(string) or not os.stat(string).st_size>0:
            raise ValueError("Invalid input file.")
        print("Loading reference")
#         seqdict = SeqIO.index(string, "fasta")
        seqdict = SeqIO.to_dict(SeqIO.parse(open(string), 'fasta'))
        print("Reference loaded") 
        return seqdict
    
    parser=argparse.ArgumentParser("""Script to prepare a GTF for the pipeline; it will perform the following operations:
    1- add the "transcript" feature
    2- sort by coordinates
    3- check the strand""")
    parser.add_argument("--fasta", type=to_seqio, required=True,
                        help="Genome FASTA file. Required." )
    parser.add_argument("-s", "--strand-specific", dest="strand_specific", action="store_true", default=False,
                        help="Flag. If set, monoexonic transcripts will be left on their strand rather than being moved to the unknown strand.")
    parser.add_argument("-l", "--lenient", action="store_true", default=False,
                        help="Flag. If set, transcripts with mixed +/- splices will not cause exceptions but rather be annotated as problematic.")
    parser.add_argument("gff", type=to_gff, help="Input GTF file.")
    parser.add_argument("out", default=sys.stdout, nargs='?', type=argparse.FileType('w'),
                        help="Output file. Default: STDOUT.")
    args=parser.parse_args()
    
    exon_lines=dict()
    
    for row in args.gff:
        if row.feature!="exon": continue
        if row.transcript not in exon_lines:
            exon_lines[row.transcript]=[]
        exon_lines[row.transcript].append(row)
    
    transcripts=[]
    
    for tid in exon_lines:
        lines=exon_lines[tid]
        transcript_line = lines[0]
        transcript_line.feature="transcript"
        transcript_line.start=min(r.start for r in lines)
        transcript_line.end=max(r.end for r in lines)
        transcript_object=transcript_checker(transcript_line, args.fasta, lenient=args.lenient, strand_specific=args.strand_specific)
        for line in lines:
            transcript_object.addExon(line)
        try:
            transcript_object.finalize()
            transcript_object.check_strand()
        except exceptions.IncorrectStrandError:
            continue
        except exceptions.InvalidTranscript:
            continue
        transcripts.append(transcript_object)
        
    for obj in sorted(transcripts):
        print(obj.__str__(to_gtf=True), file=args.out)
    
if __name__=="__main__": main()