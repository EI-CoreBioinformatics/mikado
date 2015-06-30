import csv
from logging import handlers as log_handlers
import logging
import collections
import argparse
import bisect
import operator
import itertools
from collections import namedtuple

from shanghai_lib.scales.result_storer import result_storer
from shanghai_lib.loci_objects.transcript import transcript
import shanghai_lib.exceptions
from shanghai_lib.scales.accountant import accountant


class assigner:
    
    def __init__(self, genes:dict, positions:collections.defaultdict, args:argparse.Namespace, stat_calculator:accountant):
        
        self.args = args
        self.queue_handler = log_handlers.QueueHandler(self.args.log_queue)
        self.logger = logging.getLogger("assigner")
        self.logger.addHandler(self.queue_handler)
        if args.verbose:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.INFO)
        self.logger.propagate=False
        
        self.genes = genes
        self.positions = positions
        self.gene_matches = dict()
        for gid in genes:
            self.gene_matches[gid]=dict()
            for tid in genes[gid]:
                self.gene_matches[gid][tid.id]=[]

        self.indexer = collections.defaultdict(list).fromkeys(self.positions)
        
        for chrom in positions:
            self.indexer[chrom]=sorted(self.positions[chrom].keys())
        
        self.tmap_out = open("{0}.tmap".format(args.out),'wt' )
        self.tmap_rower=csv.DictWriter(self.tmap_out, result_storer.__slots__, delimiter="\t"  )
        self.tmap_rower.writeheader()
        self.done=0
        self.stat_calculator = stat_calculator
        
    def add_to_refmap(self,result:result_storer) -> None:
        
        if result is not None and result.ccode != ("u",):
            #This is necessary for fusion genes
            if len(result.RefGene)==1:
                self.gene_matches[ result.RefGene[0] ][ result.RefId[0] ].append(result)
                    
            else: #Fusion gene
                
                for index,(gid,tid) in enumerate(zip(result.RefGene, result.RefId)):
                    new_result = result_storer( gid, tid, ["f",result.ccode[index+1]], result.TID, result.GID,
                                            result.n_prec[index], result.n_recall[index], result.n_f1[index],
                                            result.j_prec[index], result.j_recall[index], result.j_f1[index],
                                              0
                                              )
                    self.gene_matches[gid][tid].append(new_result)
        return

    def get_best(self,tr:transcript):
        
        '''This function will get the best possible assignment for each transcript.
        Fusion genes are called when the following conditions are verified:
        - the prediction intersects (at least) two transcripts in (at least) two different loci
        - the suspected fusion transcript lies on the same strand of all candidate fused genes
        - each candidate transcript has at least one fusion or 10% of its nucleotides covered by the fusion transcript.
        The 10% threshold is hard-coded in the function. 
        
        '''

        self.logger.debug("Started with {0}".format(tr.id))
        
        tr.set_logger(self.logger)
        try:
            tr.finalize()
        except shanghai_lib.exceptions.InvalidTranscript:
            return None
    #         args.queue.put_nowait("mock")
            self.logger.warn("Invalid transcript: {0}.".format(tr.id))
            self.done+=1
            self.print_tmap(None)
            return None
        
        if self.args.protein_coding is True and tr.combined_cds_length == 0:
    #         args.queue.put_nowait("mock")
            self.logger.debug("No CDS for {0}. Ignoring.".format(tr.id))
            self.done+=1
            self.print_tmap(None)
            return None
    
        if tr.chrom in self.indexer:    
            keys = self.indexer[tr.chrom]
        else:
            keys=[]        
        indexed = bisect.bisect(keys, (tr.start,tr.end) )
        
        found = []
    
        search_right = True
        search_left = True
        
    
        left_index=max(0,min(indexed, len(keys)-1)) #Must be a valid list index
        if len(keys)==0:# or left_index==0:
            search_left=False
        
        if len(keys)==0:
            search_right=False
        else:
            right_index=min(len(keys)-1, left_index+1)
            if right_index==left_index:
                search_right=False
        
        while search_left is True:
            if  keys[left_index][1]+self.args.distance<tr.start:
                search_left=False
                continue
            found.append( keys[left_index] )
            left_index-=1
            if left_index<0:
                search_left=False
        
        while search_right is True:
            if keys[right_index][0]-self.args.distance>tr.end:
                search_right=False
                continue
            found.append(keys[right_index])
            right_index+=1
            if right_index>=len(keys):
                search_right=False
    
        distances = []
        for key in found:
            distances.append( (key, max( 0, max(tr.start,key[0] ) - min(tr.end, key[1])   )   ) )
    
        distances = sorted(distances, key = operator.itemgetter(1))
    
        #Unknown transcript
        if len(found)==0 or distances[0][1]>self.args.distance:
            ccode = "u"
            best_result = result_storer( "-", "-", ccode, tr.id, ",".join(tr.parent), *[0]*6+["-"] )
            self.stat_calculator.store(tr, best_result, None)
            results = [best_result] 
    
        #Polymerase run-on
        else:
            if distances[0][1]>0:
                match = self.positions[tr.chrom][distances[0][0]][0]
                results = []
                for other in match:
                    results.append( self.calc_compare(tr, other))
                
                best_result = sorted(results, key=operator.attrgetter("distance"))[0]
            else:
                matches=list(filter(lambda x: x[1]==0, distances))
                
                if len(matches)>1:
                    #Possible fusion
                    matches = [self.positions[tr.chrom][x[0]][0] for x in matches]
                    
                    strands = set(x.strand for x in matches)
                    if len(strands)>1 and tr.strand in strands:
                        matches = list(filter(lambda match: match.strand == tr.strand, matches))
                    if len(matches)==0:
                        raise ValueError("I filtered out all matches. This is wrong!")
                    
                    best_fusion_results = []
                    results = [] #Final results
                    dubious = [] #Necessary for a double check.
                    
                    def get_f1(result):
                        return result.j_f1[0], result.n_f1[0]
                    
                    for match in matches:
                        m_res = sorted([self.calc_compare(tr, tra) for tra in match], reverse=True, key=get_f1  )
                        #A fusion is called only if I have one of the following conditions:
                        #the transcript gets one of the junctions of the other transcript
                        #the exonic overlap is >=10% (n_recall)_
                        
                        if m_res[0].j_f1[0]==0 and m_res[0].n_recall[0]<10:
                            dubious.append(m_res)
                            continue
                        results.extend(m_res)
                        best_fusion_results.append(m_res[0])

                    def dubious_getter(dubious):
                        getter=operator.attrgetter( "j_f1", "n_f1" )
                        return getter(dubious[0])

                    if len(results)==0:
                        #I have filtered out all the results, because I only match partially the reference genes
                        dubious = sorted( dubious, key=dubious_getter   )
                        results=dubious[0]
                        best_fusion_results = [results[0]]
                        
                    values = []
                    for key in result_storer.__slots__:
                        if key in ["GID", "TID", "distance"]:
                            values.append(getattr(best_fusion_results[0],key))
                        elif key=="ccode":
                            if len(best_fusion_results)>1:
                                values.append(tuple( ["f"]+[ getattr(x,key)[0]  for x in best_fusion_results])   ) #Add the "f" ccode for fusions
                            else:
                                values.append(tuple( getattr(best_fusion_results[0],key))   )
                        else:
                            val = tuple([ getattr(x,key)[0]  for x in best_fusion_results])   
                            assert len(val)==len(best_fusion_results)
                            assert type(val[0]) is not tuple, val
                            values.append(val)
                        
                    best_result = result_storer(*values)
                else:
                    match =  self.positions[tr.chrom][matches[0][0]][0]
                    results = sorted([self.calc_compare(tr, tra) for tra in match], reverse=True, key=operator.attrgetter( "j_f1", "n_f1" )  )
                    self.logger.debug(results)
                    assert len(results) == len(match.transcripts)
                    best_result = results[0]
    #     args.queue.put_nowait(result)
    #     args.refmap_queue.put_nowait(result)
    #     args.stats_queue.put_nowait((tr,result))
    #     return result
           
        for result in results:
            self.add_to_refmap( result)
        self.logger.debug("Finished with {0}".format(tr.id))
        self.print_tmap(best_result)
        self.done+=1
        return best_result
    
    def finish(self):
        self.logger.info("Finished parsing, total: {0} transcripts.".format(self.done))
        self.refmap_printer()
        self.stat_calculator.print_stats(self.genes)
    
    
    def calc_compare(self,tr:transcript, other:transcript) ->  result_storer:
    
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
        - _    Complete match, for monoexonic transcripts (nucleotide F1>=95% - i.e. min(precision,recall)>=90.4%
        - m    Exon overlap between two monoexonic transcripts
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
    
        other_exon = None
    
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
                    if nucl_f1>=0.95 and other.strand==tr.strand:
                        other_exon=other.exons[0]
                        ccode="_"
                    elif nucl_precision==1:
                        ccode="c" #contained
                    else:
                        ccode="m" #just a generic exon overlap b/w two monoexonic transcripts
        
        if ccode in ("e","o","c") and tr.strand is not None and other.strand is not None and tr.strand!=other.strand:
            ccode="x"
    
        result = result_storer(other.id, ",".join(other.parent), ccode, tr.id, ",".join(tr.parent), 
                         round(nucl_precision*100,2), round(100*nucl_recall,2),round(100*nucl_f1,2),
                         round(junction_precision*100,2), round(100*junction_recall,2), round(100*junction_f1,2),
                         distance
                         )
        if ccode is None:
            raise ValueError("Ccode is null;\n{0}".format(  repr(result)))
     
        self.stat_calculator.store(tr, result, other_exon)
    
        return result
    
    
    def print_tmap(self, res:result_storer):
        if self.done % 10000 == 0 and self.done>0:
            self.logger.info("Done {0:,} transcripts".format(self.done))
        elif self.done % 1000 == 0 and self.done>0:
            self.logger.debug("Done {0:,} transcripts".format(self.done))
        if res is not None:
            if type(res) is not result_storer:
                self.logger.exception("Wrong type for res: {0}".format(type(res)))
                self.logger.exception(repr(res))
                raise ValueError
            else:
                self.tmap_rower.writerow(res._asdict())
    

    def refmap_printer(self) -> None:
    
        '''Function to print out the best match for each gene.'''
        self.logger.info("Starting printing RefMap")
        with open( "{0}.refmap".format(self.args.out), 'wt' ) as out:
            fields = ["RefId", "ccode", "TID", "GID",  "RefGene", "best_ccode", "best_TID", "best_GID"  ]
            out_tuple = namedtuple("refmap", fields)
            
            rower=csv.DictWriter(out, fields, delimiter="\t"  )
            rower.writeheader()
            
            for gid in sorted(self.gene_matches.keys()):
                
                rows=[]
                best_picks = []
                best_pick = None
                assert len(self.gene_matches[gid].keys())>0
                for tid in sorted(self.gene_matches[gid].keys()):
                    if len(self.gene_matches[gid][tid])==0:
                        row = tuple([ gid, tid, "NA", "NA", "NA" ])
                    else:
                        best = sorted(self.gene_matches[gid][tid], key=operator.attrgetter( "j_f1", "n_f1" ), reverse=True)[0]
                        best_picks.append(best)
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
                        assert row[2]!="NA"
                        row = out_tuple( row[0], row[2], row[3], row[4], row[1], ",".join(best_pick.ccode), best_pick.TID, best_pick.GID   )
                    else:
                        row = out_tuple( row[0], "NA", "NA", "NA", row[1], "NA", "NA", "NA"  )
                    rower.writerow(row._asdict())
        self.logger.info("Finished printing RefMap")
        return None