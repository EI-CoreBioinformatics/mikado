#!/usr/bin/env python3


import sys,argparse,os
sys.path.append(
                os.path.dirname(os.path.dirname(__file__))
                )
from shanghai_lib.parsers import GTF


# genes=set()
# for mrna in filter(lambda x: "-mRNA-" in x, ids):
#     genes.add(re.sub("-mRNA.*","",mrna))

def main():
    parser=argparse.ArgumentParser("Quick utility to retrieve (or exclude) records from a GTF file.")
    parser.add_argument('-v', dest='reverse', action="store_true", default=False,
                        help="Flag. If selected, it excludes the selected ids.")
    parser.add_argument('ids', type=argparse.FileType('r'), help="The file with the ids.")
    parser.add_argument('gtf', type=argparse.FileType('r'), help="The GTF file to analyze.")
    parser.add_argument("out", nargs="?", default=sys.stdout,
                        type=argparse.FileType("w"),help="Output file. Default: STDOUT")
    args=parser.parse_args()

    ids=set([f.rstrip() for f in args.ids])

    for record in GTF.GTF(args.gtf):
        if not record: continue
        if record.transcript in ids:           
            if not args.reverse: print(record, file=args.out)
        elif args.reverse: print(record, file=args.out)
        
    # if record.feature=="gene" and record.id in genes: print(record)

if __name__=='__main__': main()
