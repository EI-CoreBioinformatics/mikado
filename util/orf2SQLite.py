import sys,os,argparse
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from loci_objects import orf
from Bio import SeqIO

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--fasta", default=None)
    parser.add_argument("bed12")
    parser.add_argument("db")
    args=parser.parse_args()
    args.fasta=SeqIO.index(args.fasta, "fasta")
    
    serializer=orf.orfSerializer(args.bed12, args.db, fasta_index=args.fasta)
    serializer.serialize()
    
if __name__=="__main__": main()