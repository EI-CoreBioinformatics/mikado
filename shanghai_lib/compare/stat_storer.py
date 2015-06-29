import argparse
from logging import handlers as log_handlers
import logging
import itertools

from shanghai_lib.loci_objects.transcript import transcript
from shanghai_lib.compare.result_storer import result_storer  
import collections

class stat_storer:
    
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
        self.singles = dict()
        self.intron_chains = dict()
        
        
        for gene in genes:
            for tr in genes[gene]:
                self.ref_transcript_num+=1
                if tr.chrom not in self.ref_introns:
                    self.exons[tr.chrom]=dict([ ("+", set()), ("-", set())  ]  )
                    self.introns[tr.chrom]=dict([ ("+", set()), ("-", set())  ]  )
                    self.intron_chains[tr.chrom]=dict([ ("+", set()), ("-", set())  ]  )
                if tr.strand is None: s="+"
                else: s=tr.strand


                #Bitye array of the form: 
                # 0b1: in reference
                # 0b10: in prediction
                # 0b100: internal
                # 0b1000: border
                # 0b10000: single
                 
                num_exons = len(tr.exons)-1

                for index,exon in enumerate(tr.exons):
                    if exon not in self.exons:
                        self.exons[tr.chrom][s][exon] = 1
                    if num_exons==0:
                        self.exons[tr.chrom][s][exon] |= 0b10000
                        self.singles[tr.chrom][s][exon] = 1
                    else:
                        if index==0:
                            self.starts[tr.chrom][s][exon[1]] = 1
                            self.exons[tr.chrom][s][exon] |= 0b1000
                        elif index == num_exons:
                            self.ends[tr.chrom][s][exon[0]] = 1
                            self.exons[tr.chrom][s][exon] |= 0b1000
                        else:
                            self.exons[tr.chrom][s][exon]
                 
                if len(tr.exons) == 1:
                    if tr.exons[0] not in self.exons[tr.chrom][s]: 
                        self.exons[tr.chrom][s][tr.exons[0]] = 0
                    else:
                        self.exons[tr.chrom][s][tr.exons[0]][8] = True 
                        
                else:
                    if tr.exons[0] not in self.exons[tr.chrom][s]: 
                        self.exons[tr.chrom][s][tr.exons[0]] = [ True, False, False, False,
                                                                True, False,
                                                                False, False,
                                                                False, False ]
                    else:
                        self.exons[tr.chrom][s][tr.exons[0]][4] = True
                    if tr.exons[-1] not in self.exons[tr.chrom][s]:
                        self.exons[tr.chrom][s][tr.exons[-1]] = [ True, False, True, False, False, False  ]
                    else:
                        self.exons[tr.chrom][s][tr.exons[-1]][0] = True
                        
                    if len(tr.exons)>2:
                        for exon in tr.exons[1:-1]:
                            if exon not in self.exons[tr.chrom][s]:
                                self.exons[tr.chrom][s][exon] = [ False, False, True, False, False, False  ]
                            else:
                                self.exons[tr.chrom][s][2] = True
                    # 0 in reference
                    # 1 in prediction        
                    for intron in tr.introns:
                        self.introns[tr.chrom][s][intron] = [True, False]

    def store(self, tr:transcript, other_exons=list, results:[result_storer]):

        if tr.strand is None: s="+"
        else: s=tr.strand

        if tr.chrom not in self.exons:
            self.exons[tr.chrom]=dict([ ("+", set()), ("-", set())  ]  )
            self.introns[tr.chrom]=dict([ ("+", set()), ("-", set())  ]  )
            self.intron_chains[tr.chrom]=dict([ ("+", set()), ("-", set())  ]  )
            
        if len(tr.exons)>1:
            if tr.exons[0] not in self.exons[tr.chrom][s]:
                self.exons[tr.chrom][s] = [True, False, False, False, True, False, False]
            else:
                self.exons[tr.chrom][s]

                



    def print_stats(self, genes):
 
        self.ref_exons = dict([("internal", dict()), ("starting", dict()), ("terminal", dict()), ("single", dict()) ])  
        self.ref_introns = dict()
        self.ref_intron_chains = dict()
         
        self.ref_transcript_num = 0

 
    
        found_exons = dict().fromkeys(self.ref_exons.keys())
        new_exons = dict().fromkeys(self.ref_exons.keys())
        missed_exons = dict().fromkeys(self.ref_exons.keys())
    
        found_introns = dict()
        new_introns = dict()
        missed_introns = dict()
        
        missed_intron_chains = dict()
        new_intron_chains = dict()
    
        for chrom in filter(lambda key: key not in self.pred_introns, self.ref_introns.keys()):
            for key in self.ref_exons:
                missed_exons[key][chrom]=self.ref_exons[key][chrom]
            missed_introns[chrom]=self.ref_introns[chrom]
            missed_intron_chains[chrom]=self.ref_intron_chains[chrom]
    
        for chrom in filter(lambda key: key not in self.ref_exons, self.pred_exons.keys()):
            for key in self.pred_exons:
                new_exons[key][chrom]=self.pred_exons[key][chrom]
            new_introns[chrom]=self.pred_introns[chrom]
            new_intron_chains[chrom] = self.pred_intron_chains[chrom]
    
        for chrom in filter(lambda key: key in self.ref_exons, self.pred_exons.keys()):
            
            for key in new_exons:
                found_exons[key][chrom]=dict([ ("+", set()), ("-", set())  ])
                new_exons[key][chrom]=dict([ ("+", set()), ("-", set())  ])
                missed_exons[key][chrom]=dict([ ("+", set()), ("-", set())  ])
            
            found_introns[chrom]=dict([ ("+", set()), ("-", set())  ])
            new_introns[chrom]=dict([ ("+", set()), ("-", set())  ])
            missed_introns[chrom]=dict([ ("+", set()), ("-", set())  ])
            
            for strand in ("+","-"):
                
                found_exons[]
                
                
                found_exons[chrom][strand] = set.intersection(self.pred_exons[chrom][strand], self.ref_exons[chrom][strand])
                new_exons[chrom][strand] = set.difference(self.pred_exons[chrom][strand], self.ref_exons[chrom][strand])
                missed_exons[chrom][strand] = set.difference(self.ref_exons[chrom][strand], self.pred_exons[chrom][strand])
                
                found_introns[chrom][strand] = set.intersection(self.pred_introns[chrom][strand], self.ref_introns[chrom][strand])
                self.logger.debug( "Comparing introns for {0}{1}; prediction: {2}; reference: {3}".format(
                                                                                                    chrom,
                                                                                                    strand,
                                                                                                    len(self.pred_introns[chrom][strand]),
                                                                                                    len(self.ref_introns[chrom][strand])
                                                                                                     ) )
                new_introns[chrom][strand] = set.difference(self.pred_introns[chrom][strand], self.ref_introns[chrom][strand])
                missed_introns[chrom][strand] = set.difference(self.ref_introns[chrom][strand], self.pred_introns[chrom][strand])
        
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
                self.logger.debug("{0}{1} found introns: {2}".format(chrom, strand, len(found_introns[chrom][strand])))
                prediction_introns_found_num+=len(found_introns[chrom][strand])
        for chrom in new_introns:
            for strand in new_introns[chrom]:
                self.logger.debug("{0}{1} new introns: {2}".format(chrom, strand, len(new_introns[chrom][strand])))
                prediction_introns_new_num+=len(new_introns[chrom][strand])
    
        for chrom in missed_introns:
            for strand in missed_introns[chrom]:
                self.logger.debug("{0}{1} missed introns: {2}".format(chrom, strand, len(missed_introns[chrom][strand])))
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
    
        ref_intron_chain_num = 0
        for chrom in self.ref_intron_chains:
            for strand in self.ref_intron_chains[chrom]:
                ref_intron_chain_num+=len(self.ref_intron_chains[chrom][strand])

        pred_intron_chain_num = 0
        for chrom in self.pred_intron_chains:
            for strand in self.pred_intron_chains[chrom]:
                pred_intron_chain_num+=len(self.pred_intron_chains[chrom][strand])


        if pred_intron_chain_num>0:        
            intron_chain_precision = len(self.matching_chains)/pred_intron_chain_num
        else:
            intron_chain_precision=0
        if ref_intron_chain_num>0:
            intron_chain_recall = len(self.matching_chains) /ref_intron_chain_num
        else:
            intron_chain_recall=0
        if min(intron_chain_recall,intron_chain_precision)>0:
            intron_chain_f1 = 2*(intron_chain_recall*intron_chain_precision)/(intron_chain_recall+intron_chain_precision)
        else:
            intron_chain_f1 =0
    
        found_bases_num = 0
        new_bases_num = 0
        missing_bases_num = 0
        
        
        #For each chromosome, generate a set of integer values for both reference and prediction
        #These are the bases annotated in ref and pred, and using set.intersection/difference
        #we are able to derive the precision/recall stats
        for chrom in set.union(set(self.ref_exons.keys()), set(self.pred_exons.keys())  ):
            if chrom not in self.pred_exons:
                for strand in self.ref_exons[chrom]:
                    missing_bases_num += len(set(itertools.chain(*[range(x[0], x[1]+1) for x in self.ref_exons[chrom][strand]]  )))
            elif chrom not in self.ref_exons:
                for strand in self.pred_exons[chrom]:
                    new_bases_num += len(set(itertools.chain(*[range(x[0], x[1]+1) for x in self.pred_exons[chrom][strand]]  )))
            else:
                for strand in self.pred_exons[chrom]:
                    ref_bases = set(itertools.chain(*[range(x[0], x[1]+1) for x in self.ref_exons[chrom][strand]]  ))
                    pred_bases = set(itertools.chain(*[range(x[0], x[1]+1) for x in self.pred_exons[chrom][strand]]  ))
                    found_bases_num += len(set.intersection(pred_bases,ref_bases )  )
                    new_bases_num +=  len(set.difference(pred_bases,ref_bases )  )
                    missing_bases_num +=  len(set.difference(ref_bases, pred_bases )  )
                    
        self.logger.debug("Matching bases: {0:,}".format(found_bases_num))
        self.logger.debug("New bases: {0:,}".format(new_bases_num))
        self.logger.debug("Missed bases: {0:,}".format(missing_bases_num))
        
        if found_bases_num+new_bases_num>0:
            bases_prec = found_bases_num/(found_bases_num+new_bases_num)
        else:
            bases_prec = 0
        bases_recall = found_bases_num/(found_bases_num+missing_bases_num)
        if max(bases_recall,bases_prec)>0:
            bases_f1 = 2*(bases_recall*bases_prec)/(bases_recall+bases_prec)
        else:
            bases_f1 = 0
        
        transcript_recall = len(self.found_ref_transcripts)/self.ref_transcript_num
        if len(self.pred_transcripts)>0:
            transcript_prec = len(self.found_pred_transcripts)/len(self.pred_transcripts)
        else:
            transcript_prec = 0
        if max(transcript_prec,transcript_recall)>0:
            transcript_f1 = 2*(transcript_recall*transcript_prec)/(transcript_recall+transcript_prec)
        else:
            transcript_f1 = 0
        
        gene_recall = len(self.found_ref_genes)/self.ref_gene_num
        if len(self.pred_genes)>0:
            gene_prec = len(self.found_pred_genes)/len(self.pred_genes)
        else:
            gene_prec = 0
        if max(gene_prec, gene_recall)>0:
            gene_f1 = 2*(gene_recall*gene_prec)/(gene_recall+gene_prec)
        else:
            gene_f1 = 0
        
        ref_introns_num=0
        for chrom in self.ref_introns:
            for strand in self.ref_introns[chrom]:
                ref_introns_num+=len(self.ref_introns[chrom][strand])
        
        
        with open("{0}.stats".format(self.args.out),'wt') as out:
            
            print("Command line:\n{0:>10}".format(self.args.commandline), file=out)
            print(self.ref_transcript_num, "reference RNAs in", self.ref_gene_num, "genes", file=out )
            print( len(self.pred_transcripts) , "predicted RNAs in ", len(self.pred_genes), "genes", file=out  )
            
            print("-"*20, "|   Sn |   Sp |   F1 |", file=out  )
            print("           {0} {1:.2f}  {2:.2f}  {3:.2f}".format("Base level:",  bases_recall*100, bases_prec*100, bases_f1*100 ) , file=out   )
            print("           {0} {1:.2f}  {2:.2f}  {3:.2f}".format("Exon level:",  exon_recall*100, exon_prec*100, exon_f1*100 )  , file=out  )
            print("         {0} {1:.2f}  {2:.2f}  {3:.2f}".format("Intron level:",  intron_recall*100, intron_prec*100, intron_f1*100 )  , file=out  )
            print("   {0} {1:.2f}  {2:.2f}  {3:.2f}".format("Intron chain level:",  intron_chain_recall*100, 
                                                         intron_chain_precision*100, intron_chain_f1*100 )  , file=out  )
            print("     {0} {1:.2f}  {2:.2f}  {3:.2f}".format("Transcript level:",  transcript_recall*100, transcript_prec*100, transcript_f1*100 )  , file=out  )
            print("           {0} {1:.2f}  {2:.2f}  {3:.2f}".format("Gene level:",  gene_recall*100, gene_prec*100, gene_f1*100 )  , file=out  )
            print(file=out)
            print(" Matching intron chains: {0}".format(sum(self.matching_chains.values())), file=out)
            print(" Matched intron chains: {0}".format(len(self.matching_chains)), file=out)
            
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
            
            
        self.logger.removeHandler(self.queue_handler)
        self.queue_handler.close()
        return
