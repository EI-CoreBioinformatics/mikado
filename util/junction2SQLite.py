import sys,os,argparse
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from loci_objects import junction

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--fai", default=None)
    parser.add_argument("bed12")
    parser.add_argument("db")
    args=parser.parse_args()
    
    serializer=junction.junctionSerializer(args.bed12, args.db, fai=args.fai)
    serializer.serialize()
    
if __name__=="__main__": main()