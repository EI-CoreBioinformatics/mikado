#!/usr/bin/env python3


#from loci_objects.sublocus import sublocus
from mikado_lib.loci_objects.transcript import Transcript
import re

'''Quick script to automate the generation of metrics definition from the files.'''

def main():
    metrics = Transcript.get_available_metrics()

    print()

    for metric in ["tid", "parent", "score"]+sorted(filter(lambda x: x not in ("tid", "parent", "score", metrics), metrics)):
        docstring=getattr(Transcript,metric).__doc__
        if docstring is None:
            docstring=''
        print( "- *{0}*:".format(metric), re.sub(" +", " ", re.sub("\n", " ", docstring)), sep="\t")
       
    print() 

if __name__=='__main__': main()