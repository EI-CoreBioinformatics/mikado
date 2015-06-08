import sys,os,argparse
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from shanghai_lib import orf
from Bio import SeqIO

def main():
    parser=argparse.ArgumentParser("Simple script to serialise ORF BED files into the SQLite DB.")
    parser.add_argument("--fasta", default=None)
    parser.add_argument("bed12")
    parser.add_argument("db")
    args=parser.parse_args()
    if args.fasta is not None:
        args.fasta=SeqIO.index(args.fasta, "fasta")
    
    serializer=orf.orfSerializer(args.bed12, args.db, fasta_index=args.fasta)
    serializer.serialize()
    
if __name__=="__main__": main()