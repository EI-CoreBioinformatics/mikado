import sys,os,argparse
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from shanghai_lib.serializers import junction

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--fai", default=None)
    parser.add_argument("--max-objects", dest="max_objects", default=10**5, type=int)
    parser.add_argument("bed12")
    parser.add_argument("db")
    args=parser.parse_args()
    
    serializer=junction.junctionSerializer(args.bed12, args.db, fai=args.fai, maxobjects=args.max_objects)
    serializer()
    
if __name__=="__main__": main()