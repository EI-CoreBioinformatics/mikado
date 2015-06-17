import sys,argparse,os
import collections
import operator
import random
sys.path.append(
                os.path.dirname(
                                os.path.dirname(__file__)
                                ))
import itertools
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

def get_best(positions:dict, indexer:dict, tr:transcript, args:argparse.Namespace):
    
    keys = indexer[tr.chrom]
    if len(keys)==0:
        ccode = "u"
        match = None
        return tr.id, ccode, match, 0
        
    indexed = bisect.bisect(keys, (tr.start,tr.end) )
    
    found = []

    search_right = True
    search_left = True
    

    left_index=max(0,min(indexed, len(keys)-1)) #Must be a valid list index
    if left_index==0: search_left = False
    
    right_index=min(len(keys)-1, left_index+1)
    if right_index==left_index: search_right=False
    
    while search_left is True and keys[left_index][1]+args.distance>=tr.start:
        found.append( keys[left_index] )
        left_index-=1
        if left_index<0: break
    
    while search_right is True and keys[right_index][0]-args.distance<=tr.end:
        found.append(keys[right_index])
        right_index+=1
        if right_index>=len(keys): break
        

    distances = []
    for key in found:
#         print(key)
        distances.append( (key, max( 0, max(tr.start,key[0] ) - min(tr.end, key[1])   )   ) )

    distances = sorted(distances, key = operator.itemgetter(1))  
#     print(distances)
    #Polymerase run-on
    if len(found)==0:
        ccode = "u"
        match = None
        return tr.id, ccode, match, indexed
    
    if distances[0][1]>0:
        match = random.choice( positions[tr.chrom][key]  ).id
        ccode = "p"
        return tr.id, ccode, match, indexed, distances[0][1]
#         else:
#             match=None
#             ccode="u"
#             return tr.id, ccode, match, indexed, distances[0][1]
    else:
        matches = []
        for d in distances:
            if d[1]==0: matches.append(d)
            else: break
            
        if len(matches)>1:
            match=",".join(str(positions[tr.chrom][key[0]][0].id) for key in matches )
            ccode = "f"
        else:
            match = random.choice( positions[tr.chrom][key]  )
            ccode = "m"
            print(tr.id, match.id, ccode)
            
            for tra in match:
                print( tr.id, tra.id, calc_compare(tr, tra)  )
                
            match=match.id
            
        return tr.id, ccode, match, indexed, distances[0][1]
    

def calc_compare(tr:transcript, other:transcript) ->  [float, float, float, float, float, float, str]:
    '''Function to compare two transcripts and determine a ccode.''' 
    
    nucl_overlap = len(set.intersection(
                                        set( itertools.chain( range(x[0], x[1]+1) for x in tr.exons  ) ),
                                        set( itertools.chain( range(x[0], x[1]+1) for x in other.exons  ) ),
                                    
                                    )
                       )
    
    nucl_recall = nucl_overlap/other.cdna_length   # Sensitivity
    nucl_precision = nucl_overlap/tr.cdna_length 
    if max(nucl_recall, nucl_precision) == 0:
        nucl_f1 = 0
    else:
        nucl_f1 = 2*(nucl_recall*nucl_precision)/(nucl_recall + nucl_precision) 

    if min(tr.exon_num, other.exon_num)>1:
        assert min(len(tr.splices), len(other.splices))>0, (tr.introns, tr.splices)
        junction_overlap = len( set.intersection(set(tr.splices), set(other.splices))   )
        junction_recall = junction_overlap / len(other.splices)
        junction_precision = junction_overlap / len(tr.splices)
        if max(junction_recall, junction_precision)>0:
            junction_f1 = 2*(junction_recall*junction_precision)/(junction_recall+junction_precision)
        else:
            junction_f1 = 0

    else:
        junction_overlap=junction_f1=junction_precision=junction_recall=0
    
    if junction_f1 == 1:
        ccode = "=" #We have recovered all the junctions
    elif tr.start>other.end or tr.end<other.start: # Outside the transcript - polymerase run-on
        ccode = "p" 
 
    elif min(tr.exon_num, other.exon_num)>1:
        if junction_precision == 1 and junction_recall < 1:
            ccode = "c" # all the junctions are correct, but we are missing some .. contained
        elif junction_recall == 1 and junction_precision<1:
            ccode = "n" # we have recovered all the junctions AND added some other junctions of our own
        elif junction_recall>0 and junction_precision>0:
            ccode = "j"
        elif (junction_recall==0 and junction_precision==0):
            if nucl_f1>0:
                ccode = "o" 
            else:
                if nucl_overlap == 0:
                    if other.start<tr.start<other.end: #The only explanation for no nucleotide overlap and no nucl overlap is that it is inside an intron
                        ccode = "I"
                    elif tr.start<other.start<tr.end:
                        ccode="K" #reverse intron retention
    else:
        if tr.exon_num==1 and other.exon_num>1:
            if nucl_overlap>0:
                ccode = "e"
            elif  other.start < tr.start < other.end:
                ccode = "i" #Monoexonic fragment inside an intron
        elif tr.exon_num>1 and other.exon_num==1:
            if nucl_overlap>0:
                ccode = "h" #Extension
            else:
                ccode = "k" #Reverse intron retention - it's the annotated model which is inside the intron!
        else:
            ccode="o" #just a generic exon overlap

    return nucl_precision, nucl_recall, nucl_f1, junction_precision, junction_recall, junction_f1, ccode


class gene:
    
    def __init__(self, tr:transcript, gid=None):
        
        self.chrom, self.start, self.end, self.strand = tr.chrom, tr.start, tr.end, tr.strand
        self.id = gid
        self.transcripts = dict()
        self.transcripts[tr.id]=tr
        
    def add(self, tr:transcript):
        self.start=min(self.start, tr.start)
        self.end = max(self.end, tr.end)
        self.transcripts[tr.id]=tr
        
    def __getitem__(self, tid:str) -> transcript:
        return self.transcripts[tid]
    
    def finalize(self):
        to_remove=set()
        for tid in self.transcripts:
            try:
                self.transcripts[tid].finalize()
            except:
                to_remove.add(tid)
                
        for k in to_remove: del self.transcripts[k]
    
    def __str__(self):
        return " ".join(self.transcripts.keys())
    
    def __iter__(self) -> transcript:
        '''Iterate over the transcripts attached to the gene.'''
        return iter(self.transcripts.values())
    

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
            tr = transcript(row)
            if row.gene not in transcripts:
                transcripts[row.gene] = gene(tr, gid=row.gene)
            transcripts[row.gene].add(tr)
            assert tr.id in transcripts[row.gene].transcripts
        else:
            
            transcripts[row.gene][row.transcript].addExon(row)
        
    for tr in transcripts:
        transcripts[tr].finalize()
        key = (transcripts[tr].start,transcripts[tr].end)
        if key not in positions[transcripts[tr].chrom]: 
            positions[transcripts[tr].chrom][key]=[]
        positions[transcripts[tr].chrom][key].append(transcripts[tr])
    
    indexer = collections.defaultdict(list).fromkeys(positions)
    for chrom in positions:
        indexer[chrom]=sorted(positions[chrom].keys())

    currentTranscript = None
    for row in args.prediction:
        if row.header is True:
            continue
        if row.is_transcript is True:
            if currentTranscript is not None:
                currentTranscript.finalize()
                print(*get_best(positions, indexer, currentTranscript, args))
            currentTranscript=transcript(row)
        else:
            currentTranscript.addExon(row)

    currentTranscript.finalize()
    print(*get_best(positions, indexer, currentTranscript, args))

if __name__=='__main__': main()