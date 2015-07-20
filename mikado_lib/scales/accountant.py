import argparse
from logging import handlers as log_handlers
import logging
# import itertools
from mikado_lib.loci_objects.transcript import transcript
from mikado_lib.scales.result_storer import result_storer  
# import collections
import operator
import collections

class accountant:
    
    '''This class stores the data necessary to calculate the final statistics - base and exon Sn/Sp/F1 etc.'''
    
    def __init__(self, genes:dict, args:argparse.Namespace):

        '''Class constructor. It requires:
        - The reference gene dictionary (containing the various transcript instances)
        - A namespace (like those provided by argparse) containing the parameters for the run.        
        '''

        if hasattr(args,"log_queue"):
            self.queue_handler = log_handlers.QueueHandler(args.log_queue)
        else:
            self.queue_handler = logging.NullHandler
        self.logger = logging.getLogger("stat_logger")
        self.logger.addHandler(self.queue_handler)
        if args.verbose:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.INFO)
        self.logger.propagate=False
        self.logger.debug("Started with stat printing")
        
        self.args = args
     
#         
#         self.pred_introns = dict()
#         self.pred_intron_chains = dict()
#         self.matching_intron_chains = collections.Counter()
#         self.found_genes = set()
#         self.found_transcripts = collections.Counter()
#         self.pred_transcripts = set()
#         self.pred_genes = set()
        
        self.introns = dict()
        #This is the easy part
        self.exons = dict()
        self.starts = dict()
        self.ends = dict()
        self.intron_chains = collections.Counter()
        self.ref_genes = dict()
        self.pred_genes = dict()
        
        for gene in genes:
            self.ref_genes[gene] = dict()
            for tr in genes[gene]:
                self.ref_genes[gene][tr.id] = 0b00
#                 self.ref_transcript_num+=1
                if tr.chrom not in self.introns:
                    self.exons[tr.chrom]=dict([ ("+", dict()), ("-", dict())  ]  )
                    self.starts[tr.chrom]=dict([ ("+", dict()), ("-", dict())  ]  )
                    self.ends[tr.chrom]=dict([ ("+", dict()), ("-", dict())  ]  )
                    
                    self.introns[tr.chrom]=dict([ ("+", dict()), ("-", dict())  ]  )
                    self.intron_chains[tr.chrom]=dict([ ("+", dict()), ("-", dict())  ]  )
                if tr.strand is None: s="+"
                else: s=tr.strand

                # 0b000001: in reference
                # 0b000010: in prediction
                # 0b000100: single
                # 0b001000: internal
                # 0b010000: border
                # 0b100000: single match lenient
                
                if tr.exon_num == 1:
                    self.exons[tr.chrom][s][tr.exons[0]] = 0b00000
                    self.exons[tr.chrom][s][tr.exons[0]] |= 0b1
                    self.exons[tr.chrom][s][tr.exons[0]] |= 0b100
                else:
                    for intron in tr.introns:
                        self.introns[tr.chrom][s][intron] = 0b01
                    if tuple(tr.introns) not in self.intron_chains[tr.chrom][s]: 
                        self.intron_chains[tr.chrom][s][tuple(tr.introns)] = [set(), set()]
                    self.intron_chains[tr.chrom][s][tuple(tr.introns)][0].add(tr.id)
                         
                    for index,exon in enumerate(tr.exons):
                        if exon not in self.exons[tr.chrom][s]: 
                            self.exons[tr.chrom][s][exon] = 0b00000
                        self.exons[tr.chrom][s][exon] |= 0b01
                        if index==0:
                            self.exons[tr.chrom][s][exon] |= 0b010000
                            self.starts[tr.chrom][s][exon[1]] = 0b01
                        elif index==tr.exon_num-1:
                            self.exons[tr.chrom][s][exon] |= 0b010000
                            self.ends[tr.chrom][s][exon[0]] = 0b01
                        else:
                            self.exons[tr.chrom][s][exon] |= 0b01000



    def store(self, tr:transcript, result:result_storer, other_exon:tuple):

        '''Add exons introns intron chains etc. to the storage.
        A transcript is considered a perfect match if it has junction_f1==100 and nucleotide_f1==100;
        a lenient match if it has junction_f1==100 and nucleotide_f1>95,
        i.e. min(nucleotide_precision, nucleotide_recall)>90.4 
        '''

        for parent in tr.parent:
            if parent not in self.pred_genes:
                self.pred_genes[parent]=dict()
            if tr.id not in self.pred_genes[parent]: 
                self.pred_genes[parent][tr.id]=0b00

        if tr.strand is None: s="+"
        else: s=tr.strand

        
        if result.ccode!=("u",):
            self.pred_genes[parent][tr.id]|=0b100
            for refid,refgene in zip(result.RefId, result.RefGene):
                self.ref_genes[refgene][refid] |= 0b100

            if len(result.j_f1)==1:
                for refid,refgene,junc_f1,nucl_f1 in zip(result.RefId, result.RefGene, result.j_f1, result.n_f1):
                    if junc_f1==100 and nucl_f1>=95:
                        for parent in tr.parent:
                            if nucl_f1==100:
                                self.pred_genes[parent][tr.id] |= 0b01
                            self.pred_genes[parent][tr.id] |= 0b10
                        self.ref_genes[refgene][refid] |= 0b10
                        if nucl_f1==100:
                            self.ref_genes[refgene][refid] |= 0b01

        if tr.chrom not in self.exons:
            self.exons[tr.chrom]=dict([ ("+", dict()), ("-", dict())  ]  )
            self.starts[tr.chrom]=dict([ ("+", dict()), ("-", dict())  ]  )
            self.ends[tr.chrom]=dict([ ("+", dict()), ("-", dict())  ]  )
            self.introns[tr.chrom]=dict([ ("+", dict()), ("-", dict())  ]  )
            self.intron_chains[tr.chrom]=dict([ ("+", dict()), ("-", dict())  ]  )
            
        if tr.exon_num>1:
            ic_key =tuple(tr.introns) 
            if result.ccode==("=",):
                assert ic_key in self.intron_chains[tr.chrom][s]
                assert result.RefId[0] in self.intron_chains[tr.chrom][s][ic_key][0]
            
            if ic_key not in self.intron_chains[tr.chrom][s]:
                self.intron_chains[tr.chrom][s][ic_key]=[set(),set()]
            self.intron_chains[tr.chrom][s][ic_key][1].add(tr.id)

            for intron in tr.introns:
                if intron not in self.introns[tr.chrom][s]:
                    self.introns[tr.chrom][s][intron]=0b0
                self.introns[tr.chrom][s][intron] |= 0b10
            
            for index,exon in enumerate(tr.exons):
                if exon not in self.exons[tr.chrom][s]:
                    self.exons[tr.chrom][s][exon]=0b0
                self.exons[tr.chrom][s][exon] |= 0b10 #set it as "in prediction" 
                if index == 0:
                    self.exons[tr.chrom][s][exon] |= 0b010000
                    if exon[1] not in self.starts[tr.chrom][s]:
                        self.starts[tr.chrom][s][exon[1]]=0b0
                    self.starts[tr.chrom][s][exon[1]] |= 0b10
                elif index == tr.exon_num-1:
                    self.exons[tr.chrom][s][exon] |= 0b010000
                    if exon[0] not in self.ends[tr.chrom][s]:
                        self.ends[tr.chrom][s][exon[0]] = 0b00
                    self.ends[tr.chrom][s][exon[0]] |= 0b10
                else:
                    self.exons[tr.chrom][s][exon] |= 0b01000
        else:
            exon=tr.exons[0]
            if exon not in self.exons[tr.chrom][s]:
                self.exons[tr.chrom][s][exon]=0b0
            self.exons[tr.chrom][s][exon] |= 0b10
            self.exons[tr.chrom][s][exon] |= 0b100
            if other_exon is not None:
                assert type(other_exon) is tuple
                for refgene, refid in zip(result.RefGene, result.RefId):
                    for exon in self.ref_genes[refgene][refid].exons:
                        assert exon in self.exons[tr.chrom], (refid, exon, self.ref_genes[refgene][refid].exons)
                
                assert other_exon in self.exons[tr.chrom][s], (tr.id, tr.exons, other_exon)
                if not 0b100 & self.exons[tr.chrom][s][other_exon]:
                    self.exons[tr.chrom][s][other_exon] |= 0b100 #Check the other exon is marked as single 
#                 self.exons[tr.chrom][s][exon] |= 0b100000
                self.exons[tr.chrom][s][other_exon] |= 0b100000

    def print_stats(self, genes):
 
        num_ref_transcripts = sum(len(self.ref_genes[gid]) for gid in self.ref_genes) 
        num_pred_transcripts = sum(len(self.pred_genes[gid]) for gid in self.pred_genes)
 
        bases_common = 0
        bases_reference = 0
        bases_prediction = 0
 
        exon_pred_stringent = 0
        exon_ref_stringent = 0
        exon_common_stringent = 0

        exon_ref_lenient = 0
        exon_pred_lenient = 0
        exon_common_lenient = 0
        missed_lenient = []
 
        for chrom in self.starts:
            for strand in self.starts[chrom]:
                for start in self.starts[chrom][strand]:
                    exon_common_lenient += (0b1 & self.starts[chrom][strand][start]) & ((0b10 & self.starts[chrom][strand][start])>>1) 
                    exon_ref_lenient +=  0b01 & self.starts[chrom][strand][start]
                    exon_pred_lenient +=  (0b10 & self.starts[chrom][strand][start])>>1
                    if not (0b10 & self.starts[chrom][strand][start])>>1:
                        missed_lenient.append(("{0}{1}".format(chrom,strand),"start",start))
                    
        starts_common = exon_common_lenient
        starts_ref = exon_ref_lenient
        starts_pred = exon_pred_lenient
        self.logger.debug("Starts {0}".format([starts_common, starts_ref, starts_pred]))
                    
        for chrom in self.ends:
            for strand in self.ends[chrom]:
                for end in self.ends[chrom][strand]:
                    exon_common_lenient += (0b01 & self.ends[chrom][strand][end]) & ((0b10 & self.ends[chrom][strand][end])>>1)
                    exon_ref_lenient +=  0b01 & self.ends[chrom][strand][end]
                    exon_pred_lenient +=  (0b10 & self.ends[chrom][strand][end])>>1
                    if not (0b10 & self.ends[chrom][strand][end])>>1:
                        missed_lenient.append(("{0}{1}".format(chrom,strand),"end",end))

        ends_common = exon_common_lenient-starts_common
        ends_ref = exon_ref_lenient-starts_ref
        ends_pred = exon_pred_lenient-starts_pred
                    
        self.logger.debug("Ends {0}".format([ends_common, ends_ref, ends_pred]))

        intron_common = 0
        intron_ref = 0
        intron_pred = 0

        for chrom in self.introns:
            for strand in self.introns[chrom]:
                for intron in self.introns[chrom][strand]:
                    intron_common += (0b01 & self.introns[chrom][strand][intron]) & ((0b10 & self.introns[chrom][strand][intron])>>1) 
                    intron_ref +=  0b01 & self.introns[chrom][strand][intron]
                    intron_pred +=  (0b10 & self.introns[chrom][strand][intron])>>1
                    
        if intron_ref>0:
            intron_recall = intron_common/intron_ref
        else:
            intron_recall = 0
        if intron_pred>0:
            intron_precision = intron_common / intron_pred
        else:
            intron_precision =0
        if max(intron_precision, intron_recall)>0:
            intron_f1 = 2*(intron_recall*intron_precision)/(intron_recall + intron_precision)
        else:
            intron_f1 = 0


        intron_chains_common_nonred = 0
        intron_chains_common_ref = set()
        intron_chains_common_pred = set()
        intron_chains_ref_nonred = 0
        intron_chains_pred_nonred = 0
        intron_chains_pred = 0
        intron_chains_ref = 0
        for chrom in self.intron_chains:
            for strand in self.intron_chains[chrom]:
                for _,intron_val in self.intron_chains[chrom][strand].items():
                    if len(intron_val[0])>0 and len(intron_val[1])>0:
                        intron_chains_common_nonred+=1
                        intron_chains_common_ref.update(intron_val[0])
                        intron_chains_common_pred.update(intron_val[1])
                    if len(intron_val[0])>0:
                        intron_chains_ref_nonred +=1
                        intron_chains_ref+=len(intron_val[0])
                    if len(intron_val[1])>0:
                        intron_chains_pred_nonred +=1
                        intron_chains_pred+=len(intron_val[1])
                    
        intron_chains_common_ref = len(intron_chains_common_ref)
        intron_chains_common_pred = len(intron_chains_common_pred)

        self.logger.debug("Intron chains:\n\treference\t{0}\t{1}\n\tprediction\t{2}\t{3}\n\tcommon\t{4} {5} {6}".format(
                                                                                                     intron_chains_ref,
                                                                                                     intron_chains_ref_nonred,
                                                                                                     intron_chains_pred,
                                                                                                     intron_chains_pred_nonred,
                                                                                                     intron_chains_common_nonred,
                                                                                                     intron_chains_common_ref,
                                                                                                     intron_chains_common_pred,
                                                                                                     ))

        if intron_chains_ref_nonred>0:
            intron_chains_recall = intron_chains_common_nonred/intron_chains_ref_nonred
        else:
            intron_chains_recall = 0
        if intron_chains_pred_nonred>0:
            intron_chains_precision = intron_chains_common_nonred/intron_chains_pred_nonred
        else:
            intron_chains_precision =0
        if max(intron_chains_precision, intron_chains_recall)>0:
            intron_chains_f1 = 2*(intron_chains_recall*intron_chains_precision)/(intron_chains_recall + intron_chains_precision)
        else:
            intron_chains_f1 = 0

        self.logger.debug("Intron chain: {0}\t{1}\t{2}".format(intron_chains_recall, intron_chains_precision, intron_chains_f1 ))

#         missed_lenient = []
        for chrom in self.exons:
            for strand in self.exons[chrom]:
                # 0b000001: in reference
                # 0b000010: in prediction
                # 0b000100: single
                # 0b001000: internal
                # 0b010000: border
                # 0b100000: single match lenient
                curr_pred_bases = set()
                curr_ref_bases = set()
                curr_exon = (-1,-1)

                for exon in sorted(self.exons[chrom][strand], key=operator.itemgetter(0)):
                    if exon[0]>curr_exon[1]:

                        bases_common+=len(set.intersection(curr_pred_bases,curr_ref_bases))
                        bases_reference += len(curr_ref_bases)
                        bases_prediction += len(curr_pred_bases)
                        curr_exon = exon
                        curr_pred_bases = set()
                        curr_ref_bases = set()
                    else:
                        curr_exon=(curr_exon[0],exon[1])
                    
                    if (0b01 & self.exons[chrom][strand][exon])==0b1: #In reference
                        curr_ref_bases.update( set(range(exon[0],exon[1])) )
                        exon_ref_stringent+=1
                        #Internal (first condition)
                        #OR
                        #Single exon
                        if ((0b001000 & self.exons[chrom][strand][exon]) == 0b001000):
                            exon_ref_lenient+=1
                        elif (0b100000 & self.exons[chrom][strand][exon]):
                            exon_ref_lenient+=1

                                                
#                         if (0b100100 & self.exons[chrom][strand][exon])==0b100100: #Single exon lenient match
#                             exon_common_lenient+=1
#                             if 0b10 ^ self.exons[chrom][strand][exon]:
#                                 exon_pred_lenient+=1
                        
                    if (0b10 & self.exons[chrom][strand][exon])==0b10: #In prediction
                        curr_pred_bases.update( set(range(exon[0],exon[1])) )
                        exon_pred_stringent+=1
                        #Either internal, or a single exon which has not a lenient single match with the reference annotation
                        if ((0b001000 & self.exons[chrom][strand][exon]) == 0b001000):
                            exon_pred_lenient+=1
                        elif (0b100000 & self.exons[chrom][strand][exon]):
                            exon_pred_lenient+=1
                        
                    if (0b11 & self.exons[chrom][strand][exon])==0b11: #In both
                        exon_common_stringent+=1
                        if ((0b001000 & self.exons[chrom][strand][exon]) == 0b001000):
                            exon_common_lenient+=1
                        elif (0b100000 & self.exons[chrom][strand][exon]):
                            exon_common_lenient+=1
                            
                bases_common+=len(set.intersection(curr_pred_bases,curr_ref_bases))
                bases_reference += len(curr_ref_bases)
                bases_prediction += len(curr_pred_bases)

        self.logger.debug("Missed lenient: {0}".format(missed_lenient))
        self.logger.debug("Exon stringent {0}".format([exon_common_stringent,exon_ref_stringent,exon_pred_stringent]))
#         assert exon_common_stringent<=min(exon_pred_stringent, exon_ref_stringent), (exon_common_stringent, exon_ref_stringent, exon_pred_stringent)
        self.logger.debug("Exon lenient {0}".format([exon_common_lenient, exon_ref_lenient, exon_pred_lenient]))
        #assert exon_common_lenient<=min(exon_pred_lenient, exon_ref_lenient), (exon_common_lenient,exon_ref_lenient, exon_pred_lenient ) 
        
        self.logger.debug("Bases:\n\treference\t{0}\n\tprediction\t{1}\n\tcommon\t{2}".format(bases_reference, bases_prediction, bases_common))
        if bases_reference>0:
            bases_recall = bases_common / bases_reference
        else:
            bases_recall = 0
        if bases_prediction>0:
            bases_prec = bases_common / bases_prediction
        else:
            bases_prec = 0
        if max(bases_prec, bases_recall)>0:
            bases_f1 = 2*(bases_prec*bases_recall)/(bases_prec+bases_recall)
        else:
            bases_f1 = 0
        
        if exon_ref_stringent>0:
            exon_stringent_recall = exon_common_stringent / exon_ref_stringent
        else:
            exon_stringent_recall = 0
        if exon_pred_stringent>0:
            exon_stringent_precision = exon_common_stringent / exon_pred_stringent
        else:
            exon_stringent_precision=0
        if max(exon_stringent_precision, exon_stringent_recall)>0:
            exon_stringent_f1 = 2*(exon_stringent_precision*exon_stringent_recall)/((exon_stringent_precision+exon_stringent_recall))
        else:
            exon_stringent_f1=0
        
        if exon_ref_lenient>0:
            exon_lenient_recall = exon_common_lenient / exon_ref_lenient
        else:
            exon_lenient_recall = 0
        if exon_pred_lenient>0:
            exon_lenient_precision = exon_common_lenient / exon_pred_lenient
        else:
            exon_lenient_precision=0
        if max(exon_lenient_precision, exon_lenient_recall)>0:
            exon_lenient_f1 = 2*(exon_lenient_precision*exon_lenient_recall)/((exon_lenient_precision+exon_lenient_recall))
        else:
            exon_lenient_f1=0
        
        
        found_ref_transcripts_stringent = 0
        found_ref_transcripts_lenient = 0
        found_ref_genes_stringent = 0
        found_ref_genes_lenient = 0
        
        missed_genes = 0
        missed_transcripts = 0
        
        ref_transcripts = 0
        for ref_gene in self.ref_genes:
            gene_match = 0b00
            gene_not_found = 0b0
            for _,val in self.ref_genes[ref_gene].items():
                ref_transcripts+=1
                found_ref_transcripts_lenient += (0b10 & val ) >> 1 
                found_ref_transcripts_stringent += 0b1 & val
                missed_transcripts +=  (0b100 & val ) >> 2 ^ 0b1 #Will evaluate to 1 if the transcript has at least one match
                gene_not_found ^= (0b100 & val ) >> 2 ^ 0b1
                gene_match |= val 
            found_ref_genes_stringent += gene_match & 0b1
            found_ref_genes_lenient += (gene_match & 0b10) >> 1
            missed_genes += gene_not_found
        
        self.logger.debug("Found ref transcripts:\n\tstringent\t{0}\n\tlenient\t{1}\n\ttotal\t{2}".format(found_ref_transcripts_stringent, found_ref_transcripts_lenient, ref_transcripts)   )
        
        found_pred_transcripts_stringent = 0
        found_pred_transcripts_lenient = 0
        found_pred_genes_stringent = 0
        found_pred_genes_lenient = 0
        pred_transcripts = 0
        
        novel_genes = 0
        novel_transcripts = 0
        
        for pred_gene in self.pred_genes:
            gene_match = 0b000
            gene_not_found = 0b0
            for _,val in self.pred_genes[pred_gene].items():
                pred_transcripts+=1
                found_pred_transcripts_lenient += (0b10 & val ) >> 1 
                found_pred_transcripts_stringent += 0b1 & val
                novel_transcripts += (0b100 & val ) >> 2 ^ 0b1 #Will evaluate to 1 if the transcript has at least one match
                gene_not_found ^= (0b100 & val ) >> 2 ^ 0b1
                gene_match |= val 
            found_pred_genes_stringent += gene_match & 0b1
            found_pred_genes_lenient += (gene_match & 0b10)>>1
            novel_genes += gene_not_found
        
        self.logger.debug("Found pred transcripts:\n\tstringent\t{0}\n\tlenient\t{1}\n\ttotal\t{2}".format(found_pred_transcripts_stringent, found_pred_transcripts_lenient, pred_transcripts)   )

        #Transcript level        
        if ref_transcripts>0:
            transcript_recall_stringent = found_ref_transcripts_stringent/ref_transcripts
            transcript_recall_lenient = found_ref_transcripts_lenient/ref_transcripts
        else:
            transcript_recall_stringent = transcript_recall_lenient = 0
        if pred_transcripts>0:
            transcript_precision_stringent = found_pred_transcripts_stringent/pred_transcripts
            transcript_precision_lenient = found_pred_transcripts_lenient/pred_transcripts
        else:
            transcript_precision_stringent = transcript_precision_lenient =0
            
        if max(transcript_precision_stringent,transcript_recall_stringent)>0:
            transcript_f1_stringent = 2*(transcript_precision_stringent*transcript_recall_stringent)/(transcript_precision_stringent+transcript_recall_stringent)
        else:
            transcript_f1_stringent = 0
            
        if max(transcript_precision_lenient,transcript_recall_lenient)>0:
            transcript_f1_lenient = 2*(transcript_precision_lenient*transcript_recall_lenient)/(transcript_precision_lenient+transcript_recall_lenient)
        else:
            transcript_f1_lenient = 0

        #Gene level
        ref_genes = len(self.ref_genes)
        pred_genes = len(self.pred_genes)            
        if ref_genes>0:
            gene_recall_stringent = found_ref_genes_stringent/ref_genes
            gene_recall_lenient = found_ref_genes_lenient/ref_genes
        else:
            gene_recall_stringent = gene_recall_lenient = 0
        if pred_genes>0:
            gene_precision_stringent = found_pred_genes_stringent/pred_genes
            gene_precision_lenient = found_pred_genes_lenient/pred_genes
        else:
            gene_precision_stringent = gene_precision_lenient =0
            
        if max(gene_precision_stringent,gene_recall_stringent)>0:
            gene_f1_stringent = 2*(gene_precision_stringent*gene_recall_stringent)/(gene_precision_stringent+gene_recall_stringent)
        else:
            gene_f1_stringent = 0
        
        if max(gene_precision_lenient,gene_recall_lenient)>0:
            gene_f1_lenient = 2*(gene_precision_lenient*gene_recall_lenient)/(gene_precision_lenient+gene_recall_lenient)
        else:
            gene_f1_lenient = 0

        
        with open("{0}.stats".format(self.args.out),'wt') as out:
            
            print("Command line:\n{0:>10}".format(self.args.commandline), file=out)
            print(num_ref_transcripts, "reference RNAs in", len(self.ref_genes), "genes", file=out )
            print( num_pred_transcripts , "predicted RNAs in ", len(self.pred_genes), "genes", file=out  )
            
            print("-"*30, "|   Sn |   Sp |   F1 |", file=out  )
            print("                     {0} {1:.2f}  {2:.2f}  {3:.2f}".format("Base level:",  bases_recall*100, bases_prec*100, bases_f1*100 ) , file=out   )
            print("         {0} {1:.2f}  {2:.2f}  {3:.2f}".format("Exon level (stringent):",  exon_stringent_recall*100, exon_stringent_precision*100, exon_stringent_f1*100 )  , file=out  )
            print("           {0} {1:.2f}  {2:.2f}  {3:.2f}".format("Exon level (lenient):",  exon_lenient_recall*100, exon_lenient_precision*100, exon_lenient_f1*100 )  , file=out  )
            print("                   {0} {1:.2f}  {2:.2f}  {3:.2f}".format("Intron level:",  intron_recall*100, intron_precision*100, intron_f1*100 )  , file=out  )
            print("             {0} {1:.2f}  {2:.2f}  {3:.2f}".format("Intron chain level:",  intron_chains_recall*100, 
                                                         intron_chains_precision*100, intron_chains_f1*100 )  , file=out  )
            print("   {0} {1:.2f}  {2:.2f}  {3:.2f}".format("Transcript level (stringent):",  transcript_recall_stringent*100, transcript_precision_stringent*100, transcript_f1_stringent*100 )  , file=out  )
            print("     {0} {1:.2f}  {2:.2f}  {3:.2f}".format("Transcript level (lenient):",  transcript_recall_lenient*100, transcript_precision_lenient*100, transcript_f1_lenient*100 )  , file=out  )
            print("         {0} {1:.2f}  {2:.2f}  {3:.2f}".format("Gene level (stringent):",  gene_recall_stringent*100, gene_precision_stringent*100, gene_f1_stringent*100 )  , file=out  )
            print("           {0} {1:.2f}  {2:.2f}  {3:.2f}".format("Gene level (lenient):",  gene_recall_lenient*100, gene_precision_lenient*100, gene_f1_lenient*100 )  , file=out  )
            print(file=out)
            print("         Matching intron chains: {0}".format(intron_chains_common_pred), file=out)
            print("          Matched intron chains: {0}".format(intron_chains_common_ref), file=out)
            
            print("", file = out)            
#             print(" Matched intron chains: {0}".format(len(self.matching_chains)), file=out)
            
            print("       {0} {1}/{2}  ({3:.2f}%)".format("Missed exons (stringent):",
                                                       exon_ref_stringent-exon_common_stringent,
                                                       exon_ref_stringent,
                                                       100*(exon_ref_stringent-exon_common_stringent)/(exon_ref_stringent) if exon_ref_stringent>0 else 0 ),
                  file=out)
            print("        {0} {1}/{2}  ({3:.2f}%)".format("Novel exons (stringent):",
                                                       exon_pred_stringent-exon_common_stringent,
                                                       exon_pred_stringent,
                                                       100*(exon_pred_stringent-exon_common_stringent)/(exon_pred_stringent) if exon_pred_stringent>0 else 0 ),
                  file=out)
            print("         {0} {1}/{2}  ({3:.2f}%)".format("Missed exons (lenient):",
                                                       exon_ref_lenient-exon_common_lenient,
                                                       exon_ref_lenient,
                                                       100*(exon_ref_lenient-exon_common_lenient)/(exon_ref_lenient) if exon_ref_lenient>0 else 0 ),
                  file=out)
            print("          {0} {1}/{2}  ({3:.2f}%)".format("Novel exons (lenient):",
                                                       exon_pred_lenient-exon_common_lenient,
                                                       exon_pred_lenient,
                                                       100*(exon_pred_lenient-exon_common_lenient)/(exon_pred_lenient) if exon_pred_lenient>0 else 0 ),
                  file=out)

            
            print("                 {0} {1}/{2}  ({3:.2f}%)".format("Missed introns:",
                                                       intron_ref-intron_common,
                                                       intron_ref,
                                                       100*(intron_ref-intron_common)/intron_ref if intron_ref>0 else 0 ),
                  file=out)

            
            print("                  {0} {1}/{2}  ({3:.2f}%)".format("Novel introns:",
                                                       intron_pred-intron_common,
                                                       intron_pred,
                                                       100*(intron_pred-intron_common)/intron_pred if intron_pred>0 else 0 ),
                  file=out)
            print("", file = out)
            print("             {0} {1}/{2} ({3:.2f}%)".format("Missed transcripts:",
                                                                missed_transcripts,
                                                                num_ref_transcripts,
                                                                100*missed_transcripts/num_ref_transcripts if num_ref_transcripts>0 else 0
                                                                ), file=out)
            print("              {0} {1}/{2} ({3:.2f}%)".format("Novel transcripts:",
                                                                novel_transcripts,
                                                                num_pred_transcripts,
                                                                100*novel_transcripts/num_pred_transcripts if num_pred_transcripts>0 else 0
                                                                ), file=out)

            print("                   {0} {1}/{2} ({3:.2f}%)".format("Missed genes:",
                                                                missed_genes,
                                                                len(self.ref_genes),
                                                                100*missed_genes/len(self.ref_genes) if len(self.ref_genes)>0 else 0 
                                                                ), file=out)
            print("                    {0} {1}/{2} ({3:.2f}%)".format("Novel genes:",
                                                                novel_genes,
                                                                len(self.pred_genes),
                                                                100*novel_genes/len(self.pred_genes) if len(self.pred_genes)>0 else 0
                                                                ), file=out)


            
            
        self.logger.removeHandler(self.queue_handler)
        self.queue_handler.close()
        return
