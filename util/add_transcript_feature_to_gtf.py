#!/usr/bin/env python3


import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from loci_objects.GTF import GTF
from copy import deepcopy
import operator
import argparse

class obj(object): pass

current=obj()
current.transcript=None

rows=[]

parser=argparse.ArgumentParser("Script to add a transcript feature to e.g. Cufflinks GTFs")
parser.add_argument("gtf", type=argparse.FileType("r"),
                    help="Input GTF")
parser.add_argument("out", default=sys.stdout, nargs="?",
                    type=argparse.FileType("w"),
                    help="Output file. Default: stdout.")
args=parser.parse_args()

args.gtf.close()

for record in GTF(args.gtf.name):
        if current.transcript!=record.transcript:
                if current.transcript!=None:
                        print(current, file=args.out)
                        exon_no=0
                        for row in filter(lambda x: x.feature=="exon", sorted(rows, key=operator.attrgetter("start"))):
                                exon_no+=1
                                row.info["exon_number"]=exon_no
                                print(row, file=args.out)
                        exon_no=0
                        for row in filter(lambda x: x.feature=="CDS", sorted(rows, key=operator.attrgetter("start"))):
                                exon_no+=1
                                row.info["exon_number"]=exon_no
                                print(row, file=args.out)
                rows=[record]
                current=deepcopy(record)
                current.feature="transcript"
                
        else:
                current.stop=max(current.stop,record.stop)
                current.start=min(current.start, record.start)
                rows.append(record)

print(current, file=args.out)
for row in rows: print(row, file=args.out)
