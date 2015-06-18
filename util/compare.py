import sys,argparse,os
import collections
import operator
import random
import asyncio
import queue
import threading
import concurrent.futures
import time
sys.path.append(
                os.path.dirname(
                                os.path.dirname(__file__)
                                ))
import csv
from shanghai_lib.loci_objects.transcript import transcript
# import shanghai_lib.exceptions
# from shanghai_lib.loci_objects.abstractlocus import abstractlocus
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
        result = args.formatter( tr.id, "NA", ccode, *[0]*6  )
        args.queue.put(result)
        return
        
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
        args.queue.put( args.formatter( tr.id, "NA", ccode, *[0]*6   ) )
        return

    if distances[0][1]>0:
        match = random.choice( positions[tr.chrom][key]  ).id
        ccode = "p"
        args.queue.put(args.formatter( tr.id, match, ccode, *[0]*6))
        return 
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
            matches = [positions[tr.chrom][x[0]][0] for x in matches]
            res = []
            for match in matches:
                m_res = sorted([calc_compare(tr, tra, args.formatter) for tra in match], reverse=True, key=operator.attrgetter( "j_f1", "n_f1" )  )
                res.append(m_res[0])
                
            fields = [tr.id]
            fields.append( ",".join( x[1] for x in res  )  )
            ccode = ",".join(["f"] + [x.ccode for x in res]) 
            fields.append(ccode)
            for field in args.formatter._fields[3:]:
                fields.append( ",".join( str(getattr(x,field)) for x in res  ) )
            
            result = args.formatter(*fields)
            args.queue.put(result)
            return 
            
        else:
            match =  positions[tr.chrom][matches[0][0]][0]
            res = sorted([calc_compare(tr, tra, args.formatter) for tra in match], reverse=True, key=operator.attrgetter( "j_f1", "n_f1" )  )
            best = res[0]
            args.queue.put(best)

            return 
    
    
def calc_compare(tr:transcript, other:transcript, formatter:collections.namedtuple) ->  collections.namedtuple:

    '''Function to compare two transcripts and determine a ccode.
    Available ccodes (from Cufflinks documentation):
    
    - =    Complete match
    - c    Contained
    - j    Potentially novel isoform (fragment): at least one splice junction is shared with a reference transcript
    - e    Single exon transfrag overlapping a reference exon and at least 10 bp of a reference intron, indicating a possible pre-mRNA fragment.
    - i    A *monoexonic* transfrag falling entirely within a reference intron
    - o    Generic exonic overlap with a reference transcript
    - p    Possible polymerase run-on fragment (within 2Kbases of a reference transcript)
    - u    Unknown, intergenic transcript
    - x    Exonic overlap with reference on the opposite strand
    - s    An intron of the transfrag overlaps a reference intron on the opposite strand (likely due to read mapping errors)

    Please note that the description for i is changed from Cufflinks.

    We also provide the following additional classifications:
    
    - f    gene fusion
    - n    Potentially novel isoform, where all the known junctions have been confirmed and we have added others as well
    - I    *multiexonic* transcript falling completely inside a known transcript
    - h    the transcript is multiexonic and extends a monoexonic reference transcript
    - O    Reverse generic overlap - the reference is monoexonic and overlaps the prediction
    - K    Reverse intron retention - the annotated gene model retains an intron compared to the prediction

    ''' 
    
    tr_nucls = list()
    for x in tr.exons:
        tr_nucls.extend(range(x[0], x[1]+1))
    other_nucls = list()
    for x in other.exons:
        other_nucls.extend(range(x[0], x[1]+1))
    
    
    nucl_overlap = len(set.intersection(
                                        set(other_nucls), set(tr_nucls)
                                    )
                       )
                       
    assert nucl_overlap<=min(other.cdna_length,tr.cdna_length), (tr.id, tr.cdna_length, other.id, other.cdna_length, nucl_overlap)
    
    nucl_recall = nucl_overlap/other.cdna_length   # Sensitivity
    nucl_precision = nucl_overlap/tr.cdna_length 
    if max(nucl_recall, nucl_precision) == 0:
        nucl_f1 = 0
    else:
        nucl_f1 = 2*(nucl_recall*nucl_precision)/(nucl_recall + nucl_precision) 

    if min(tr.exon_num, other.exon_num)>1:
        assert min(len(tr.splices), len(other.splices))>0, (tr.introns, tr.splices)
        one_junction_confirmed = any(intron in other.introns for intron in tr.introns)
        junction_overlap = len( set.intersection(set(tr.splices), set(other.splices))   )
        junction_recall = junction_overlap / len(other.splices)
        junction_precision = junction_overlap / len(tr.splices)
        if max(junction_recall, junction_precision)>0:
            junction_f1 = 2*(junction_recall*junction_precision)/(junction_recall+junction_precision)
        else:
            junction_f1 = 0

    else:
        junction_overlap=junction_f1=junction_precision=junction_recall=0
    
    ccode = None
    
    if junction_f1 == 1  or (tr.exon_num==other.exon_num==1 and tr.start==other.start and tr.end==other.end):
        if tr.strand==other.strand or tr.strand is None:
            ccode = "=" #We have recovered all the junctions
        else:
            ccode = "c"
        
    elif tr.start>other.end or tr.end<other.start: # Outside the transcript - polymerase run-on
        ccode = "p" 
 
    elif nucl_precision==1:
        if tr.exon_num==1 or (tr.exon_num>1 and junction_precision==1):
            ccode="c"
        
 
    if ccode is None:
        if min(tr.exon_num, other.exon_num)>1:
            if junction_recall == 1 and junction_precision<1:
                ccode = "n" # we have recovered all the junctions AND added some other junctions of our own
            elif junction_recall>0 and 0<junction_precision<1:
                if one_junction_confirmed is True:
                    ccode = "j"
                else:
                    ccode = "o"
            elif junction_precision==1:
                ccode = "c"
                if nucl_precision<1:
                    for intron in other.introns:
                        if intron in tr.introns: continue
                        if intron[1]<tr.start: continue
                        elif intron[0]>tr.end: continue
                        if tr.start<intron[0] and intron[1]<tr.end:
                            ccode="j"
                            break
                    
                    
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
                if nucl_precision<1 and nucl_overlap>0:
                    outside = max( other.start-tr.start,0  )+max(tr.end-other.end,0) #Fraction outside
                    if tr.cdna_length-nucl_overlap-outside>10:
                        ccode = "e"
                    else:
                        ccode = "o"
                elif nucl_overlap>0:
                    ccode = "o"
                elif nucl_recall==0 and  other.start < tr.start < other.end:
                    ccode = "i" #Monoexonic fragment inside an intron
            elif tr.exon_num>1 and other.exon_num==1:
                if nucl_recall==1:
                    ccode = "h" #Extension
                else:
                    ccode = "O" #Reverse generic overlap
            else:
                if nucl_precision==1:
                    ccode="c" #just a generic exon overlap
                else:
                    ccode="o"
    
    if ccode in ("e","o","c") and tr.strand is not None and other.strand is not None and tr.strand!=other.strand:
        ccode="x"

    result = formatter(tr.id, other.id, ccode,
                     round(nucl_precision*100,2), round(100*nucl_recall,2),round(100*nucl_f1,2),
                     round(junction_precision*100,2), round(100*junction_recall,2), round(100*junction_f1,2),
                     ) 

    if ccode is None:
        raise ValueError(result)

    return result


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

#@asyncio.coroutine
def printer( args ):
   
    rower=csv.DictWriter(args.out, args.formatter._fields, delimiter="\t"  )
    rower.writeheader()
    print("Starting to listening")
    while True:
        try:
            res = args.queue.get(block=False)
            if res=="EXIT": return
            rower.writerow(res._asdict())
        except queue.Empty:
            time.sleep(0.1)
            continue


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
    parser.add_argument("-t", "--threads", default=1, type=int)
    parser.add_argument("-o","--out", default=sys.stdout, type = argparse.FileType("w") )

    
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
        elif row.header is True:
            continue
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
    print("Finished parsing reference")

    formatter=collections.namedtuple("compare",
                             [ "TID", "RefId", "ccode", "n_prec", "n_recall", "n_f1",
                              "j_prec", "j_recall", "j_f1",
                              ]
                             , verbose=False)
    args.formatter = formatter

    args.queue = queue.Queue()
#     loop = asyncio.get_event_loop()
    pp=threading.Thread(target=printer, args=(args,), name="printing_thread")
    
    print("Initialising printer")
    print("Initialised printer")
    pool = concurrent.futures.ProcessPoolExecutor(args.threads)
    
    currentTranscript = None
    for row in args.prediction:
        if row.header is True:
            continue
        if row.is_transcript is True:
            if currentTranscript is not None:
                currentTranscript.finalize()
                pool.submit(get_best(positions, indexer, currentTranscript, args))
            currentTranscript=transcript(row)
        else:
            currentTranscript.addExon(row)

    currentTranscript.finalize()
    pool.submit(get_best(positions, indexer, currentTranscript, args))
    print("Finished parsing")
    pool.shutdown(wait=True)

    args.queue.put("EXIT")
    
    pp.start()
#     loop.run_until_complete(asyncio.async(printer(args)))
    pp.join()

if __name__=='__main__': main()