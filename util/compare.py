import sys,argparse,os
import collections
import operator
import random
sys.path.append(
                os.path.dirname(
                                os.path.dirname(__file__)
                                ))
import csv
from shanghai_lib.loci_objects.transcript import transcript
import shanghai_lib.exceptions
from shanghai_lib.loci_objects.abstractlocus import abstractlocus
from shanghai_lib.parsers.GTF import GTF
from shanghai_lib.parsers.GFF import GFF3
import bisect # Needed for efficient research

'''This is still an embryo. Ideally, this program would perform the following functions:

1- Define precision/recall for the annotation
2- Use a "flexible" mode to determine accuracy
3- Detect gene *fusions* as well as gene *splits*  

'''

def get_best(positions, indexer, tr, args):
    
    keys = list(indexer[tr.chrom])
    indexed = bisect(keys, (tr.start,tr.end) )
    
    found = []
    #Search left
    left_index=right_index=indexed
    while positions[tr.chrom][keys[left_index]].start+args.distance>=tr.start:
        found.append( keys[left_index] )
        left_index-=1
        
    right_index+=1
    while positions[tr.chrom][keys[right_index]].end-args.distance<=tr.end:
        found.append(keys[right_index])
        right_index+=1
        
    if len(found)==0:
        ccode = "u"
        match = None
        return ccode, match

    distances = sorted([ (key, max(0, max(tr.start,key[0] ) - min(tr.end, key[1])   )) for key in found   ], key = operator.itemgetter(1))  
    
    #Polymerase run-on
    
    
    if distances[0][1]>0:
        match = random.choice( positions[tr.chrom][key]  ).id
        ccode = "p"
        return ccode, match
    
    genes = 
    


def main():
    
    def to_gtf(string):
        '''Function to recognize the input file type and create the parser.'''
        
        if string.endswith(".gtf"):
            return GTF(string)
#         elif string.endswith('.gff') or string.endswith('.gff3'):
#             return GFF3(string)
        else:
            raise ValueError('Unrecognized file format.')
        
    
    parser=argparse.ArgumentParser('Tool to define the spec/sens of predictions vs. references.')
    input_files=parser.add_argument_group('Prediction and annotation files.')
    input_files.add_argument('-r', '--reference', type=to_gtf, help='Reference annotation file.', required=True)
    input_files.add_argument('-p', '--prediction', type=to_gtf, help='Prediction annotation file.', required=True)
    parser.add_argument('--distance', type=int, default=2000, 
                        help='Maximum distance for a transcript to be considered a polymerase run-on. Default: %(default)s')
    
    args=parser.parse_args()

    transcripts = dict()
    positions = collections.defaultdict(dict)
    if type(args.reference) is GTF:
        ref_gtf = True
    else:
        ref_gtf = False
    
    for row in args.reference:
        #Assume we are going to use GTF for the moment
        if row.is_transcript is True:
            transcripts[row.transcript]=transcript(row)
        else:
            transcripts[row.transcript].addExon(row)
        
    for tr in transcripts:
        tr = transcripts[tr].finalize()
        key = (tr.start,tr.end)
        if key not in positions[tr.chrom]: 
            positions[tr.chrom][key]=[]
        positions[tr.chrom][key].append(tr)
    
    indexer = collections.defaultdict(list).fromkeys(positions)
    for chrom in positions:
        indexer[chrom]=sorted(positions[chrom].keys())

    currentTranscript = None
    for row in args.prediction:
        if row.is_transcript is True:
            if currentTranscript is not None:
                ccode, match = get_best(positions, indexer, currentTranscript, args)
            currentTranscript=transcript(row)
        else:
            currentTranscript.addExon(row)

    get_best(positions, currentTranscript)

if __name__=='__main__': main()