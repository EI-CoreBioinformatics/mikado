#!/usr/bin/env python3

#This block is necessary to add to the PATH the parent folder, where the various modules are located
import sys
import os
path=os.path.dirname( os.path.dirname(os.path.abspath(__file__) ))
sys.path.insert(0,path)

#from loci_objects.sublocus import sublocus
from loci_objects.transcript import transcript
import re

'''Quick script to automate the generation of metrics definition from the files.'''

def main():
    metrics = transcript.get_available_metrics()

    print()

    for metric in ["tid", "parent", "score"]+sorted(filter(lambda x: x not in ("tid", "parent", "score", metrics), metrics)):
        docstring=getattr(transcript,metric).__doc__ 
        if docstring is None:
            docstring=''
        print( "- *{0}*:".format(metric), re.sub(" +", " ", re.sub("\n", " ", docstring)), sep="\t")
       
    print() 

if __name__=='__main__': main()