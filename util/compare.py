import sys,argparse,os
from collections import namedtuple
import itertools
sys.path.append(os.path.dirname(os.path.dirname( os.path.abspath(__file__)  )))
import collections
import operator
#import random
import threading
import multiprocessing
import shanghai_lib.exceptions
import csv
from shanghai_lib.loci_objects.transcript import transcript
from shanghai_lib.result_storer import result_storer
# import shanghai_lib.exceptions
# from shanghai_lib.loci_objects.abstractlocus import abstractlocus
from shanghai_lib.parsers.GTF import GTF
from shanghai_lib.parsers.GFF import GFF3
import logging
from logging import handlers as log_handlers
import bisect # Needed for efficient research


'''This is still an embryo. Ideally, this program would perform the following functions:

1- Define precision/recall for the annotation
2- Use a "flexible" mode to determine accuracy
3- Detect gene *fusions* as well as gene *splits*  

'''

def get_best(positions:dict, indexer:dict, tr:transcript, args:argparse.Namespace):
    
    queue_handler = log_handlers.QueueHandler(args.log_queue)
    logger = logging.getLogger("selector")
    logger.addHandler(queue_handler)
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    logger.propagate=False
    logger.debug("Started with {0}".format(tr.id))
    
    tr.set_logger(logger)
    try:
        tr.finalize()
    except shanghai_lib.exceptions.InvalidTranscript:
        args.queue.put_nowait("mock")
        logger.warn("Invalid transcript: {0}.".format(tr.id))
        logger.removeHandler(queue_handler)
        queue_handler.close()
        return
    
    if args.protein_coding is True and tr.combined_cds_length == 0:
        args.queue.put_nowait("mock")
        logger.debug("No CDS for {0}. Ignoring.".format(tr.id))
        logger.removeHandler(queue_handler)
        queue_handler.close()
        return

    if tr.chrom in indexer:    
        keys = indexer[tr.chrom]
    else:
        keys=[]        
    indexed = bisect.bisect(keys, (tr.start,tr.end) )
    
    found = []

    search_right = True
    search_left = True
    

    left_index=max(0,min(indexed, len(keys)-1)) #Must be a valid list index
    if left_index==0: search_left = False
    
    right_index=min(len(keys)-1, left_index+1)
    if right_index==left_index: search_right=False
    
    while search_left is True:
        if  keys[left_index][1]+args.distance<tr.start:
            break
        found.append( keys[left_index] )
        left_index-=1
        if left_index<0: break
    
    while search_right is True:
        if keys[right_index][0]-args.distance>tr.end:
            break
        found.append(keys[right_index])
        right_index+=1
        if right_index>=len(keys): break

    distances = []
    for key in found:
        distances.append( (key, max( 0, max(tr.start,key[0] ) - min(tr.end, key[1])   )   ) )

    distances = sorted(distances, key = operator.itemgetter(1))

    #Unknown transcript
    if len(found)==0 or distances[0][1]>args.distance:
        ccode = "u"
        result = result_storer( "-", "-", ccode, tr.id, ",".join(tr.parent), *[0]*6+["-"] )

    #Polymerase run-on
    else:
        if distances[0][1]>0:
            match = positions[tr.chrom][distances[0][0]][0]
            m_distances = []
            for other in match:
                dist = max(tr.start-other.end, other.start-tr.end) 
                m_distances.append((other.id, dist))
            mmatch = sorted( m_distances, key=operator.itemgetter(1)  )[0][0] #Get the transcript with the least distance
            mmatch = match[mmatch]
            
            ccode = "p"
            result = calc_compare(tr, mmatch, result_storer)
        else:
            matches=list(filter(lambda x: x[1]==0, distances))
            
            if len(matches)>1:
                matches = [positions[tr.chrom][x[0]][0] for x in matches]
                
                strands = set(x.strand for x in matches)
                if len(strands)>1 and tr.strand in strands:
                    matches = list(filter(lambda match: match.strand == tr.strand, matches))
                if len(matches)==0:
                    raise ValueError("I filtered out all matches. This is wrong!")
                
                res = []
                for match in matches:
                    m_res = sorted([calc_compare(tr, tra, result_storer) for tra in match], reverse=True, key=operator.attrgetter( "j_f1", "n_f1" )  )
                    res.append(m_res[0])
                    
                values = []
                for key in result_storer.__slots__:
                    if key in ["GID", "TID", "distance"]:
                        values.append(getattr(res[0],key))
                    elif key=="ccode":
                        if len(res)>1:
                            values.append(tuple( ["f"]+[ getattr(x,key)[0]  for x in res])   ) #Add the "f" ccode for fusions
                        else:
                            values.append(tuple( getattr(res[0],key))   )
                    else:
                        val = tuple([ getattr(x,key)[0]  for x in res])   
                        assert len(val)==len(res)
                        assert type(val[0]) is not tuple, val
                        values.append(val)
                    
                result = result_storer(*values)
            else:
                match =  positions[tr.chrom][matches[0][0]][0]
                res = sorted([calc_compare(tr, tra, result_storer) for tra in match], reverse=True, key=operator.attrgetter( "j_f1", "n_f1" )  )
                result = res[0]
    
    args.queue.put_nowait(result)
    args.refmap_queue.put_nowait(result)
    args.stats_queue.put_nowait((tr,result))
        
            
    logger.debug("Finished with {0}".format(tr.id))
    logger.removeHandler(queue_handler)
    queue_handler.close()
    return
    
def calc_compare(tr:transcript, other:transcript, formatter:result_storer) ->  result_storer:

    '''Function to compare two transcripts and determine a ccode.
    Available ccodes (from Cufflinks documentation):
    
    - =    Complete intron chain match
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
    
    - f    gene fusion - in this case, this ccode will be followed by the ccodes of the matches for each gene, separated by comma
    - _    Complete match, for monoexonic transcripts
    - n    Potentially novel isoform, where all the known junctions have been confirmed and we have added others as well
    - I    *multiexonic* transcript falling completely inside a known transcript
    - h    the transcript is multiexonic and extends a monoexonic reference transcript
    - O    Reverse generic overlap - the reference is monoexonic and overlaps the prediction
    - K    Reverse intron retention - the annotated gene model retains an intron compared to the prediction
    - P    Possible polymerase run-on fragment (within 2Kbases of a reference transcript), on the opposite strand

    ''' 
    
    tr_nucls = set(itertools.chain(*[range(x[0], x[1]+1) for x in tr.exons]))
    other_nucls = set(itertools.chain(*[range(x[0], x[1]+1) for x in other.exons]))
    
    
    nucl_overlap = len(set.intersection(
                                        other_nucls, tr_nucls
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
    distance = 0
    if junction_f1 == 1:
        if tr.strand==other.strand or tr.strand is None:
            ccode = "=" #We have recovered all the junctions
        else:
            ccode = "c" #We will set this to x at the end of the function
        
    elif (tr.exon_num==other.exon_num==1 and tr.start==other.start and tr.end==other.end):
        junction_f1 = junction_precision = junction_precision = 1 #Set to one
        if tr.strand==other.strand or tr.strand is None:
            ccode = "_" #We have recovered all the junctions
        else:
            ccode = "x"
        
    elif tr.start>other.end or tr.end<other.start: # Outside the transcript - polymerase run-on
        distance = max(tr.start-other.end, other.start-tr.end)
        if other.strand == tr.strand:
            ccode = "p"
        else:
            ccode = "P" 
        distance = max(tr.start-other.end, other.start-tr.end) 
 
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
                    ccode = "n" #Extension
                else:
                    ccode = "O" #Reverse generic overlap
            elif tr.exon_num == other.exon_num ==1:
                junction_f1 = junction_precision = junction_precision = 1 #Set to one
                if nucl_precision==1:
                    ccode="c" #just a generic exon overlap
                else:
                    ccode="o"
    
    if ccode in ("e","o","c") and tr.strand is not None and other.strand is not None and tr.strand!=other.strand:
        ccode="x"

    result = formatter(other.id, ",".join(other.parent), ccode, tr.id, ",".join(tr.parent), 
                     round(nucl_precision*100,2), round(100*nucl_recall,2),round(100*nucl_f1,2),
                     round(junction_precision*100,2), round(100*junction_recall,2), round(100*junction_f1,2),
                     distance
                     )
    if ccode is None:
        raise ValueError("Ccode is null;\n{0}".format(  repr(result)))
 

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
        assert self.strand == tr.strand
        
    def __getitem__(self, tid:str) -> transcript:
        return self.transcripts[tid]
    
    def finalize(self):
        to_remove=set()
        for tid in self.transcripts:
            try:
                self.transcripts[tid].finalize()
            except shanghai_lib.exceptions.InvalidTranscript:
                to_remove.add(tid)
                
        for k in to_remove:
            del self.transcripts[k]
    
    def remove(self, tid:str):
        del self.transcripts[tid]
        if len(self.transcripts)==0:
            self.end=None
            self.start=None
            self.chrom=None
        self.start = min(self.transcripts[tid].start for tid in self.transcripts)
        self.end = max(self.transcripts[tid].end for tid in self.transcripts)
    
    def __str__(self):
        return " ".join(self.transcripts.keys())
    
    def __iter__(self) -> transcript:
        '''Iterate over the transcripts attached to the gene.'''
        return iter(self.transcripts.values())

#@asyncio.coroutine
def printer( args ):

    logger = logging.getLogger("printer")
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
         
    queue_handler = log_handlers.QueueHandler(args.log_queue)
    logger.addHandler(queue_handler)
    logger.propagate=False
    logger.info("Started the printer process")

   
    with open("{0}.tmap".format(args.out),'wt' ) as out:
   
        rower=csv.DictWriter(out, result_storer.__slots__, delimiter="\t"  )
        rower.writeheader()

        done=0
        while True:
            res = args.queue.get()
            if res=="EXIT":
                logger.debug("Receveid EXIT signal")
                logger.info("Finished {0} transcripts".format(done))
                logger.removeHandler(queue_handler)
                queue_handler.close()
                break
            done+=1
            if done % 10000 == 0:
                logger.info("Done {0} transcripts".format(done))
            elif done % 1000 == 0:
                logger.debug("Done {0} transcripts".format(done))
            if res!="mock":
                rower.writerow(res._asdict())
            args.queue.task_done()
    
    logger.removeHandler(queue_handler)
    queue_handler.close()
    
    return

def refmap_printer(args, genes):

    '''Function to print out the best match for each gene.'''
    
    gene_matches = dict()
    for gid in genes:
        gene_matches[gid]=dict()
        for tid in genes[gid]:
            gene_matches[gid][tid.id]=[]
    
    queue_handler = log_handlers.QueueHandler(args.log_queue)
    logger = logging.getLogger("selector")
    logger.addHandler(queue_handler)
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    logger.propagate=False
    logger.debug("Started with refmapping")

    while True:
        try:
            curr_match = args.refmap_queue.get()
        except Exception as err:
            logger.exception(err)
            return
        if curr_match == "EXIT":
            break
        else:
            if curr_match.ccode == ("u",):
                continue
            #This is necessary for fusion genes
            elif len(curr_match.RefGene)==1:
                gene_matches[ curr_match.RefGene[0] ][ curr_match.RefId[0] ].append(curr_match)
            else: #Fusion gene
                for index,(gid,tid) in enumerate(zip(curr_match.RefGene, curr_match.RefId)):
                    
                    result = result_storer( gid, tid, ["f",curr_match.ccode[index]], curr_match.TID, curr_match.GID,
                                            curr_match.n_prec[index], curr_match.n_recall[index], curr_match.n_f1[index],
                                            curr_match.j_prec[index], curr_match.j_recall[index], curr_match.j_f1[index],
                                              0
                                              )
                    gene_matches[gid][tid].append(result)
        args.refmap_queue.task_done()
                
        
    with open( "{0}.refmap".format(args.out), 'wt' ) as out:
        fields = ["RefId", "ccode", "TID", "GID",  "RefGene", "best_ccode", "best_TID", "best_GID"  ]
        out_tuple = namedtuple("refmap", fields)
        
        rower=csv.DictWriter(out, fields, delimiter="\t"  )
        rower.writeheader()
        
        for gid in sorted(gene_matches.keys()):
            best_picks = []
            rows=[]
            for tid in sorted(gene_matches[gid].keys()):
                if len(gene_matches[gid][tid])==0:
                    row=(tid, gid, "NA", "NA", "NA")
                else:
                    try:
                        best = sorted(gene_matches[gid][tid], key=operator.attrgetter( "j_f1", "n_f1" ), reverse=True)[0]
                        best_picks.append(best)
                    except TypeError as err:
                        logger.exception(gene_matches[gid][tid])
                        logger.exception(err)
                        raise
                    if len(best.ccode)==1:
                        row=tuple([tid, gid, ",".join(best.ccode), best.TID, best.GID])
                    else:
                        row=tuple([tid, gid, ",".join(best.ccode), best.TID, best.GID])
                rows.append(row)
                
            if len(best_picks)>0:
                best_pick = sorted( best_picks,  key=operator.attrgetter( "j_f1", "n_f1" ), reverse=True)[0]
            else:
                best_pick = None
            for row in rows:
                if best_pick is not None:
                    row = out_tuple( row[0], row[2], row[3], row[4], row[1], ",".join(best_pick.ccode), best_pick.TID, best_pick.GID   )
                else:
                    row = out_tuple( row[0], "NA", "NA", "NA", row[1], "NA", "NA", "NA"  )
                
                rower.writerow(row._asdict())
        pass
    
    logger.removeHandler(queue_handler)
    queue_handler.close()
    return


def stat_printer(genes, args):
    '''This function will determine the global statistics'''
    
    queue_handler = log_handlers.QueueHandler(args.log_queue)
    logger = logging.getLogger("stat_logger")
    logger.addHandler(queue_handler)
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
    logger.propagate=False
    logger.debug("Started with stat printing")
    
    ref_exons = dict()
    ref_introns = dict()
     
    queue = args.stats_queue
    ref_transcript_num = 0
    for gene in genes:
        for tr in genes[gene]:
            ref_transcript_num+=1
            if tr.chrom not in ref_exons:
                ref_exons[tr.chrom]=dict([ ("+", set()), ("-", set())  ]  )
                ref_introns[tr.chrom]=dict([ ("+", set()), ("-", set())  ]  )
            if tr.strand is None: s="+"
            else: s=tr.strand
            ref_exons[tr.chrom][s].update(tr.exons)
            ref_introns[tr.chrom][s].update(tr.introns)
 

    pred_exons = dict()
    pred_introns = dict()
    pred_transcripts = set()
    pred_genes = set()
    found_ref_transcripts = set()
    found_pred_transcripts = set()
    matching_chains=set()
    new_transcripts = set()
    new_genes = set()
    found_ref_genes = set()
    found_pred_genes = set()
 
    while True:
        delivery = queue.get()
        if delivery =="EXIT":
            break
        tr, result = delivery
        
        if tr.strand is None: s="+"
        else: s=tr.strand
        if tr.chrom not in pred_exons:
            pred_exons[tr.chrom]=dict([ ("+", set()), ("-", set()) ] )
            pred_introns[tr.chrom]=dict([ ("+", set()), ("-", set()) ] )
        pred_exons[tr.chrom][s].update(tr.exons)
        if tr.exon_num>1:
            pred_introns[tr.chrom][s].update(tr.introns)
        pred_transcripts.add(tr.id)
        pred_genes.add(tr.parent[0])
        if result.ccode==("_",):
            found_ref_transcripts.add(result.RefId[0])
            found_ref_genes.add(result.RefGene[0])
            found_pred_transcripts.add(tr.id)
            found_pred_genes.add(tr.parent[0])
        elif result.ccode==("=",):
            matching_chains.add(result.RefId[0])
            ref_tr = genes[result.RefGene[0]][result.RefId[0]] 
            if ref_tr.start==tr.start and ref_tr.end == tr.end:
                found_ref_transcripts.add(result.RefId[0])
                found_ref_genes.add(result.RefGene[0])
                found_pred_transcripts.add(tr.id)
                found_pred_genes.add(tr.parent[0])
        elif result.ccode==("u",):
            new_transcripts.add(result.TID)
            new_genes.add(result.GID)
        queue.task_done()
    
    found_exons = dict()
    new_exons = dict()
    missed_exons = dict()

    found_introns = dict()
    new_introns = dict()
    missed_introns = dict()


    for chrom in filter(lambda key: key not in pred_exons, ref_exons.keys()):
        missed_exons[chrom]=ref_exons[chrom]
        missed_introns[chrom]=ref_introns[chrom]

    for chrom in filter(lambda key: key not in ref_exons, pred_exons.keys()):
        new_exons[chrom]=pred_exons[chrom]
        new_introns[chrom]=pred_introns[chrom]

    for chrom in filter(lambda key: key in ref_exons, pred_exons.keys()):
        found_exons[chrom]=dict([ ("+", set()), ("-", set())  ])
        new_exons[chrom]=dict([ ("+", set()), ("-", set())  ])
        missed_exons[chrom]=dict([ ("+", set()), ("-", set())  ])
        
        found_introns[chrom]=dict([ ("+", set()), ("-", set())  ])
        new_introns[chrom]=dict([ ("+", set()), ("-", set())  ])
        missed_introns[chrom]=dict([ ("+", set()), ("-", set())  ])
        
        for strand in ("+","-"):
            found_exons[chrom][strand] = set.intersection(pred_exons[chrom][strand], ref_exons[chrom][strand])
            new_exons[chrom][strand] = set.difference(pred_exons[chrom][strand], ref_exons[chrom][strand])
            missed_exons[chrom][strand] = set.difference(ref_exons[chrom][strand], pred_exons[chrom][strand])
            
            found_introns[chrom][strand] = set.intersection(pred_introns[chrom][strand], ref_introns[chrom][strand])
            logger.debug( "Comparing introns for {0}{1}; prediction: {2}; reference: {3}".format(
                                                                                                chrom,
                                                                                                strand,
                                                                                                len(pred_introns[chrom][strand]),
                                                                                                len(ref_introns[chrom][strand])
                                                                                                 ) )
            new_introns[chrom][strand] = set.difference(pred_introns[chrom][strand], ref_introns[chrom][strand])
            missed_introns[chrom][strand] = set.difference(ref_introns[chrom][strand], pred_introns[chrom][strand])
    
    prediction_found_num = 0
    for chrom in found_exons:
        for strand in found_exons[chrom]:
            prediction_found_num+=len(found_exons[chrom][strand])
    
    
    prediction_new_num = 0
    for chrom in new_exons:
        for strand in new_exons[chrom]:
            prediction_new_num+= len(new_exons[chrom][strand])
    
    reference_missed_num = 0
    for chrom in missed_exons:
        for strand in missed_exons[chrom]:
            reference_missed_num += len(missed_exons[chrom][strand])
    
    if prediction_found_num+prediction_new_num>0:
        exon_prec = prediction_found_num/(prediction_found_num+prediction_new_num)
    else:
        exon_prec=0
    exon_recall = prediction_found_num / (prediction_found_num+reference_missed_num)
    if max(exon_prec,exon_recall)>0:
        exon_f1 = 2*(exon_prec*exon_recall)/(exon_prec+exon_recall)
    else:
        exon_f1 = 0

    prediction_introns_found_num=0
    prediction_introns_new_num=0
    prediction_introns_missed_num=0
        
    for chrom in found_introns:
        for strand in found_introns[chrom]:
            logger.debug("{0}{1} found introns: {2}".format(chrom, strand, len(found_introns[chrom][strand])))
            prediction_introns_found_num+=len(found_introns[chrom][strand])
    for chrom in new_introns:
        for strand in new_introns[chrom]:
            logger.debug("{0}{1} new introns: {2}".format(chrom, strand, len(new_introns[chrom][strand])))
            prediction_introns_new_num+=len(new_introns[chrom][strand])

    for chrom in missed_introns:
        for strand in missed_introns[chrom]:
            logger.debug("{0}{1} missed introns: {2}".format(chrom, strand, len(missed_introns[chrom][strand])))
            prediction_introns_missed_num+=len(missed_introns[chrom][strand])

    if prediction_found_num+prediction_new_num>0:
        intron_prec = prediction_introns_found_num/(prediction_introns_found_num+prediction_introns_new_num)
    else:
        intron_prec=0
    if prediction_introns_found_num+prediction_introns_missed_num==0:
        intron_recall=0
    else:
        intron_recall = prediction_introns_found_num / (prediction_introns_found_num+prediction_introns_missed_num)
    if max(intron_prec,intron_recall)>0:
        intron_f1 = 2*(intron_prec*intron_recall)/(intron_prec+intron_recall)
    else:
        intron_f1 = 0


    found_bases_num = 0
    new_bases_num = 0
    missing_bases_num = 0
    
    
    #For each chromosome, generate a set of integer values for both reference and prediction
    #These are the bases annotated in ref and pred, and using set.intersection/difference
    #we are able to derive the precision/recall stats
    for chrom in set.union(set(ref_exons.keys()), set(pred_exons.keys())  ):
        if chrom not in pred_exons:
            for strand in ref_exons[chrom]:
                missing_bases_num += len(set(itertools.chain(*[range(x[0], x[1]+1) for x in ref_exons[chrom][strand]]  )))
        elif chrom not in ref_exons:
            for strand in pred_exons[chrom]:
                new_bases_num += len(set(itertools.chain(*[range(x[0], x[1]+1) for x in pred_exons[chrom][strand]]  )))
        else:
            for strand in pred_exons[chrom]:
                ref_bases = set(itertools.chain(*[range(x[0], x[1]+1) for x in ref_exons[chrom][strand]]  ))
                pred_bases = set(itertools.chain(*[range(x[0], x[1]+1) for x in pred_exons[chrom][strand]]  ))
                found_bases_num += len(set.intersection(pred_bases,ref_bases )  )
                new_bases_num +=  len(set.difference(pred_bases,ref_bases )  )
                missing_bases_num +=  len(set.difference(ref_bases, pred_bases )  )
                
    logger.debug("Matching bases: {0:,}".format(found_bases_num))
    logger.debug("New bases: {0:,}".format(new_bases_num))
    logger.debug("Missed bases: {0:,}".format(missing_bases_num))
    
    if found_bases_num+new_bases_num>0:
        bases_prec = found_bases_num/(found_bases_num+new_bases_num)
    else:
        bases_prec = 0
    bases_recall = found_bases_num/(found_bases_num+missing_bases_num)
    if max(bases_recall,bases_prec)>0:
        bases_f1 = 2*(bases_recall*bases_prec)/(bases_recall+bases_prec)
    else:
        bases_f1 = 0
    
    transcript_recall = len(found_ref_transcripts)/ref_transcript_num
    if len(pred_transcripts)>0:
        transcript_prec = len(found_pred_transcripts)/len(pred_transcripts)
    else:
        transcript_prec = 0
    if max(transcript_prec,transcript_recall)>0:
        transcript_f1 = 2*(transcript_recall*transcript_prec)/(transcript_recall+transcript_prec)
    else:
        transcript_f1 = 0
    
    gene_recall = len(found_ref_genes)/len(genes)
    if len(pred_genes)>0:
        gene_prec = len(found_pred_genes)/len(pred_genes)
    else:
        gene_prec = 0
    if max(gene_prec, gene_recall)>0:
        gene_f1 = 2*(gene_recall*gene_prec)/(gene_recall+gene_prec)
    else:
        gene_f1 = 0
    
    ref_introns_num=0
    for chrom in ref_introns:
        for strand in ref_introns[chrom]:
            ref_introns_num+=len(ref_introns[chrom][strand])
    
    
    with open("{0}.stats".format(args.out),'wt') as out:
        
        print(ref_transcript_num, "reference RNAs in", len(genes), "genes", file=out )
        print( len(pred_transcripts) , "predicted RNAs in ", len(pred_genes), "genes", file=out  )
        
        print("-"*20, "|   Sn |   Sp |   F1 |", file=out  )
        print("           {0} {1:.2f}  {2:.2f}  {3:.2f}".format("Base level:",  bases_recall*100, bases_prec*100, bases_f1*100 ) , file=out   )
        print("           {0} {1:.2f}  {2:.2f}  {3:.2f}".format("Exon level:",  exon_recall*100, exon_prec*100, exon_f1*100 )  , file=out  )
        print("         {0} {1:.2f}  {2:.2f}  {3:.2f}".format("Intron level:",  intron_recall*100, intron_prec*100, intron_f1*100 )  , file=out  )
        print("     {0} {1:.2f}  {2:.2f}  {3:.2f}".format("Transcript level:",  transcript_recall*100, transcript_prec*100, transcript_f1*100 )  , file=out  )
        print("           {0} {1:.2f}  {2:.2f}  {3:.2f}".format("Gene level:",  gene_recall*100, gene_prec*100, gene_f1*100 )  , file=out  )
        print(file=out)
        print(" Matching intron chains: {0}".format(len(matching_chains)), file=out)
        print("           {0} {1}/{2}  ({3:.2f}%)".format("Missed exons:",
                                                   reference_missed_num,
                                                   prediction_found_num+reference_missed_num,
                                                   100*reference_missed_num/(prediction_found_num+reference_missed_num) ),
              file=out)
        print("            {0} {1}/{2}  ({3:.2f}%)".format("Novel exons:",
                                                   prediction_new_num,
                                                   prediction_found_num+prediction_new_num,
                                                   100*prediction_new_num/(prediction_found_num+prediction_new_num) if prediction_found_num+prediction_new_num>0 else 0  ),
              file=out)
        print("          {0} {1}/{2}  ({3:.2f}%)".format("Novel introns:",
                                                   prediction_introns_new_num,
                                                   prediction_introns_new_num+prediction_introns_found_num,
                                                   100*prediction_introns_new_num/(prediction_introns_new_num+prediction_introns_found_num) if prediction_introns_new_num+prediction_introns_found_num>0 else 0  ),
              file=out)
        print("         {0} {1}/{2}  ({3:.2f}%)".format("Missed introns:",
                                                   prediction_introns_missed_num,
                                                   prediction_introns_found_num+prediction_introns_missed_num,
                                                   100*prediction_introns_missed_num/(prediction_introns_found_num+prediction_introns_missed_num) ),
              file=out)
        
        
    logger.removeHandler(queue_handler)
    queue_handler.close()
    return

def main():
    
    def to_gtf(string):
        '''Function to recognize the input file type and create the parser.'''
        
        if string.endswith(".gtf"):
            return GTF(string)
        elif string.endswith('.gff') or string.endswith('.gff3'):
            return GFF3(string)
        else:
            raise ValueError('Unrecognized file format.')
        
    
    parser=argparse.ArgumentParser('Tool to define the spec/sens of predictions vs. references.')
    input_files=parser.add_argument_group('Prediction and annotation files.')
    input_files.add_argument('-r', '--reference', type=to_gtf, help='Reference annotation file.', required=True)
    input_files.add_argument('-p', '--prediction', type=to_gtf, help='Prediction annotation file.', required=True)
    parser.add_argument('--distance', type=int, default=2000, 
                        help='Maximum distance for a transcript to be considered a polymerase run-on. Default: %(default)s')
    parser.add_argument('-pc', '--protein-coding', dest="protein_coding", action="store_true", default=False,
                        help="Flag. If set, only transcripts with a CDS (both in reference and prediction) will be considered.")
#     parser.add_argument("-t", "--threads", default=1, type=int)
    parser.add_argument("-o","--out", default="shangai_compare", type = str,
                        help = "Prefix for the output files. Default: %(default)s" )
    parser.add_argument("-l","--log", default=None, type = str)
    parser.add_argument("-v", "--verbose", action="store_true", default=False)

    
    args=parser.parse_args()

    #Flags for the parsing
    if type(args.reference) is GFF3:
        ref_gff = True
    else:
        ref_gff = False
    
    if os.path.dirname(args.out)!='' and os.path.dirname(args.out)!=os.path.dirname(os.path.abspath(".")):
        dirname = os.path.dirname(args.out)
        if os.path.exists(dirname):
            assert os.path.isdir(dirname)
        else:
            os.makedirs(dirname)
        
    
    context = multiprocessing.get_context() #@UndefinedVariable
    manager = context.Manager()
    args.queue = manager.Queue(-1)

    logger = logging.getLogger("main")
    formatter = logging.Formatter("{asctime} - {name} - {levelname} - {message}", style="{")
    args.log_queue = manager.Queue()
    args.queue_handler = log_handlers.QueueHandler(args.log_queue)
    log_queue_listener = log_handlers.QueueListener(args.log_queue, logger)
    log_queue_listener.propagate=False
    log_queue_listener.start()

    if args.log is None:
        handler = logging.StreamHandler()
    else:
        if os.path.exists(args.log):
            os.remove(args.log)
        
        handler = logging.FileHandler(args.log, mode='a')
    handler.setFormatter(formatter)
     
    if args.verbose is False:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.DEBUG)
    logger.addHandler(handler)
    
    
    queue_logger = logging.getLogger("main_queue")
    queue_logger.setLevel(logging.INFO)
    main_queue_handler = log_handlers.QueueHandler(args.log_queue)
    queue_logger.addHandler(main_queue_handler)
    
    queue_logger.propagate=False
    queue_logger.info("Start")
    queue_logger.info("Command line: {0}".format(" ".join(sys.argv)))
    
    refmap_queue = manager.Queue(-1)
    
    pp=threading.Thread(target=printer, args=(args,), name="printing_thread", daemon=True)
    pp.start()
    args.refmap_queue = refmap_queue

    queue_logger.info("Starting parsing the reference")

    genes = dict()
    positions = collections.defaultdict(dict)

    transcript2gene = dict()
    logger.handlers[0].flush()
    
    for row in args.reference:
        #Assume we are going to use GTF for the moment
        if row.header is True:
            continue
#         logger.debug(str(row))
        if row.is_transcript is True:
            queue_logger.debug("Transcript\n{0}".format(str(row)))
            tr = transcript(row)
            transcript2gene[row.id]=row.gene
            if row.gene not in genes:
                genes[row.gene] = gene(tr, gid=row.gene)
            genes[row.gene].add(tr)
            assert tr.id in genes[row.gene].transcripts
        elif row.is_exon is True:
#             logger.debug(str(row))
#             assert type(row.transcript) is list
#             logger.debug("Exon found: {0}, {1}".format(row.transcript, row.parent))
            if ref_gff is True:
                for tr in row.transcript:
#                     logger.debug(tr)
                    gid = transcript2gene[tr]
                    genes[gid][tr].addExon(row)
            else:
#                 logger.debug(row.transcript)
                try:
                    genes[row.gene][row.transcript].addExon(row)
                except KeyError as exc:
                    assert row.gene in genes
                    queue_logger.exception(exc)
                    queue_logger.exception( "Keys for {0}: {1}".format(row.gene,  genes[row.gene].transcripts.keys() ))
                    raise
        else:
            continue
#     logger.info("Finished parsing the reference")
    genes_to_remove = set()
    for gid in genes:
        genes[gid].finalize()
        if args.protein_coding is True:
            to_remove=[]
            for tid in genes[gid].transcripts:
                if genes[gid].transcripts[tid].combined_cds_length==0:
                    to_remove.append(tid)
                    logger.debug("No CDS for {0}".format(tid))
            if len(to_remove)==len(genes[gid].transcripts):
                genes_to_remove.add(gid)
                logger.debug("Noncoding gene: {0}".format(tid))
                continue
            elif len(to_remove)>0:
                for tid in to_remove: genes[gid].remove(tid) 
        key = (genes[gid].start,genes[gid].end)
        if key not in positions[genes[gid].chrom]: 
            positions[genes[gid].chrom][key]=[]
        positions[genes[gid].chrom][key].append(genes[gid])
    
    for gid in genes_to_remove: del genes[gid]
    
    indexer = collections.defaultdict(list).fromkeys(positions)
    for chrom in positions:
        indexer[chrom]=sorted(positions[chrom].keys())
    queue_logger.info("Finished preparation")

    args.stats_queue = manager.Queue(-1)

    stat_thread = threading.Thread(target=stat_printer, args=(genes, args), name="stat_thread", daemon=True)
    stat_thread.start()

    queue_logger.debug("Initialised printer")    
#     logger.debug("Initialised worker pool")

    refmap_proc = threading.Thread(target=refmap_printer, args=(args, genes ), name="refmap_printer", daemon=True)
    refmap_proc.start()
    
    def reduce_args(args):
        new_args = args.__dict__.copy()
        name = argparse.Namespace()
        del new_args["out"]
        del new_args["prediction"]
        del new_args["reference"]
        name.__dict__.update(new_args)
        return name
    
    cargs = reduce_args(args)
    
    currentTranscript = None
    for row in args.prediction:
        if row.header is True:
            continue
#         queue_logger.debug("Row:\n{0:>20}".format(str(row)))
        if row.is_transcript is True:
            queue_logger.debug("Transcript row:\n{0}".format(str(row)))
            if currentTranscript is not None:
                try:
                    get_best(positions, indexer, currentTranscript, cargs)
                except Exception as err:
                    queue_logger.exception(err)
                    log_queue_listener.enqueue_sentinel()
                    handler.close()
                    log_queue_listener.stop()
                    args.queue_handler.close()
                    raise
                
            currentTranscript=transcript(row)
        elif row.is_exon is True:
            try:
                currentTranscript.addExon(row)
            except Exception as err:
                queue_logger.exception(err)
                #In case of error, signal the threads to exit
                args.queue.put("EXIT")
                args.queue.all_tasks_done = True
                args.refmap_queue.put("EXIT")
                break
        else:
            continue

    if currentTranscript is not None:
        try:
            get_best(positions, indexer, currentTranscript, cargs)
        except Exception as err:
            queue_logger.exception(err)
            log_queue_listener.enqueue_sentinel()
            handler.close()
            log_queue_listener.stop()
            args.queue_handler.close()
            raise
 
    queue_logger.info("Finished parsing")
    args.queue.join()
    args.queue.put("EXIT")
    args.queue.all_tasks_done = True

    queue_logger.debug("Sent EXIT signal")
    pp.join()
    args.refmap_queue.put("EXIT")
    refmap_proc.join()
    queue_logger.debug("Printer process alive: {0}".format(pp.is_alive()))
    queue_logger.debug("Refmap process alive: {0}".format(refmap_proc.is_alive()))
    
    args.stats_queue.join()
    args.stats_queue.put("EXIT")
    stat_thread.join()
    
    queue_logger.info("Finished")
    log_queue_listener.enqueue_sentinel()
    handler.close()
    log_queue_listener.stop()
    args.queue_handler.close()
    return


if __name__=='__main__': main()