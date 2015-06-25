import argparse
from logging import handlers as log_handlers
import logging
import itertools

from shanghai_lib.loci_objects.transcript import transcript
from shanghai_lib.compare.result_storer import result_storer  
import collections

class stat_storer:
    
    def __init__(self, genes:dict, args:argparse.Namespace):

        self.queue_handler = log_handlers.QueueHandler(args.log_queue)
        self.logger = logging.getLogger("stat_logger")
        self.logger.addHandler(self.queue_handler)
        if args.verbose:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.INFO)
        self.logger.propagate=False
        self.logger.debug("Started with stat printing")
        
        self.ref_exons = dict([("internal", dict()), ("starting", dict()), ("terminal", dict()), ("single", dict()) ])  
        self.ref_introns = dict()
        self.ref_intron_chains = dict()
         
        self.ref_transcript_num = 0
        for gene in genes:
            for tr in genes[gene]:
                self.ref_transcript_num+=1
                if tr.chrom not in self.ref_exons:
                    for key in self.ref_exons:
                        self.ref_exons[key][tr.chrom]=dict([ ("+", set()), ("-", set())  ]  )
                    self.ref_introns[tr.chrom]=dict([ ("+", set()), ("-", set())  ]  )
                    self.ref_intron_chains[tr.chrom]=dict([ ("+", set()), ("-", set())  ]  )
                if tr.strand is None: s="+"
                else: s=tr.strand

                if tr.monoexonic is False:
                    exons=sorted(tr.exons)
                    starting = exons[0]
                    ending = exons[-1]
                    if len(exons)>2:
                        internal=set(exons[1:-1])
                    else:
                        internal=set()
                    
                    self.ref_exons["internal"][tr.chrom][s].update(internal)
                    self.ref_exons["starting"][tr.chrom][s].add(starting)
                    self.ref_exons["terminal"][tr.chrom][s].add(ending)
                                                       
                    self.ref_introns[tr.chrom][s].update(tr.introns)
                    self.ref_intron_chains[tr.chrom][s].add( tuple(sorted(tr.introns)) )
                else:
                    self.ref_exons["single"][tr.chrom][s].add(tr.exons[0])
     
        self.pred_exons = dict([("internal", dict()), ("starting", dict()), ("terminal", dict()), ("single", dict()) ])
        self.pred_introns = dict()
        self.pred_transcripts = set()
        self.pred_genes = set()
        self.found_ref_transcripts = set()
        self.found_pred_transcripts = set()
        self.matching_chains=collections.Counter()
        self.new_transcripts = set()
        self.new_genes = set()
        self.found_ref_genes = set()
        self.found_pred_genes = set()
        self.ref_gene_num = len(genes)

        self.found_exons = dict()
        self.found_border_exons=dict()
        self.new_exons = dict()
        self.missed_exons = dict()
    
        self.found_introns = dict()
        self.new_introns = dict()
        self.missed_introns = dict()



        self.pred_intron_chains = dict()

    def store(self, tr:transcript, result:result_storer, genes:dict):

        if result is None:
            return
        if tr.strand is None: s="+"
        else: s=tr.strand


        for store in (self.found_exons, self.new_exons, self.missed_exons):
            if tr.chrom not in store:
                store[tr.chrom]=dict( [ ("+", set()), ("-", set()) ] )
                
        for store in (self.found_introns, self.new_introns, self.missed_introns):
            if tr.chrom not in store:
                store[tr.chrom]=dict( [ ("+", set()), ("-", set()) ] )

        if tr.monoexonic is False:
            self.found_introns[tr.chrom][s].update( set.intersection( set(tr.introns), self.ref_introns[tr.chrom][s]  )  )
            self.new_introns[tr.chrom][s].update(set.difference( set(tr.introns), self.ref_introns[tr.chrom][s]  )  )
            
            exons = sorted(tr.exons)
            if len(exons)>2:
                internal = exons[1:-1] 
                self.found_exons[tr.chrom][s].update( set.intersection(set(internal),  self.ref_exons["internal"][tr.chrom][s]  )  )
            
            





    def print_stats(self, args:argparse.Namespace):
    
        found_exons = dict()
        new_exons = dict()
        missed_exons = dict()
    
        found_introns = dict()
        new_introns = dict()
        missed_introns = dict()
    
    
        for chrom in filter(lambda key: key not in self.pred_exons, self.ref_exons.keys()):
            missed_exons[chrom]=self.ref_exons[chrom]
            missed_introns[chrom]=self.ref_introns[chrom]
    
        for chrom in filter(lambda key: key not in self.ref_exons, self.pred_exons.keys()):
            new_exons[chrom]=self.pred_exons[chrom]
            new_introns[chrom]=self.pred_introns[chrom]
    
        for chrom in filter(lambda key: key in self.ref_exons, self.pred_exons.keys()):
            found_exons[chrom]=dict([ ("+", set()), ("-", set())  ])
            new_exons[chrom]=dict([ ("+", set()), ("-", set())  ])
            missed_exons[chrom]=dict([ ("+", set()), ("-", set())  ])
            
            found_introns[chrom]=dict([ ("+", set()), ("-", set())  ])
            new_introns[chrom]=dict([ ("+", set()), ("-", set())  ])
            missed_introns[chrom]=dict([ ("+", set()), ("-", set())  ])
            
            for strand in ("+","-"):
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
        
        
        with open("{0}.stats".format(args.out),'wt') as out:
            
            print("Command line:\n{0:>10}".format(args.commandline), file=out)
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
