#!/usr/bin/env python3

import re
import multiprocessing
import csv
import loci_objects
import os


class Creator:
    
    def __init__(self, json_conf):
        
        self.json_conf = json_conf
        
        self.input_file = self.json_conf["input"]
        self.input = self.to_input(self.input_file) #Serialize the input

        
        
    def to_input(self, string):
        '''Function to select the correct serializer for the input file (either GTF or GFF)'''
        
        assert os.path.exists(string)
        if string.endswith(".gtf"):
            parser = loci_objects.GTF.GTF
        else:
            parser = loci_objects.GFF.GFF3 
        
        checked=False
        for row in parser(string):
            if row.header: continue
            else:
                checked=True
                break
        if checked is False:
            raise TypeError("Unrecognized input file format.")
        return parser(string)
        
    def printer(self):
        
        '''Listener process that will print out the loci recovered by the analyse_locus function.'''
        
        sub_metrics=csv.DictWriter(open(args.sub_metrics,'a'), superlocus.available_metrics, delimiter="\t")
        sub_scores=csv.DictWriter(open(args.sub_scores,'a'), args.score_keys, delimiter="\t")
        locus_metrics=csv.DictWriter(open(args.locus_metrics,'a'), superlocus.available_metrics, delimiter="\t")
        locus_scores=csv.DictWriter(open(args.locus_scores,'a'), args.score_keys, delimiter="\t")
        sub_out=open(args.sub_out,'a')
        mono_out=open(args.mono_out,'a')
        locus_out=open(args.locus_out,'a') 
        with open(args.sub_metrics, "w") as out_file:
            csv_out=csv.DictWriter( out_file, superlocus.available_metrics, delimiter="\t" )
            csv_out.writeheader()
        with open(args.sub_scores,'w') as out_file:
            csv_out=csv.DictWriter( out_file, args.score_keys, delimiter="\t" )
            csv_out.writeheader()
        with open(args.locus_metrics, "w") as out_file:
            csv_out=csv.DictWriter( out_file, superlocus.available_metrics, delimiter="\t" )
            csv_out.writeheader()
        with open(args.locus_scores,'w') as out_file:
            csv_out=csv.DictWriter( out_file, args.score_keys, delimiter="\t" )
            csv_out.writeheader()
        
        
        print('##gff-version 3', file=sub_out)
        print('##gff-version 3', file=mono_out)
        print('##gff-version 3', file=locus_out)
        
        while True:
            stranded_locus=self.printer_queue.get()
            if stranded_locus is None:
                break
            sub_lines = stranded_locus.__str__(level="subloci", print_cds=not args.no_cds )
            sub_metrics_rows = [x for x in stranded_locus.print_subloci_metrics()]
            sub_scores_rows = [x for x in stranded_locus.print_subloci_scores()]
            mono_lines = stranded_locus.__str__(level="monosubloci", print_cds=not args.no_cds)
            locus_metrics_rows=[x for x in stranded_locus.print_monoholder_metrics()]
            locus_scores_rows = [x for x in stranded_locus.print_monoholder_scores()]
            locus_lines = stranded_locus.__str__(print_cds=not args.no_cds)
            for row in sub_metrics_rows: sub_metrics.writerow(row)
            for row in sub_scores_rows: sub_scores.writerow(row)
            for row in locus_metrics_rows: locus_metrics.writerow(row)
            for row in locus_scores_rows: locus_scores.writerow(row)
            print(sub_lines, file=sub_out)
            if mono_lines!='':
                print(mono_lines, file=mono_out)
            if locus_lines!='':
                print(locus_lines, file=locus_out)    
        return
    
    def analyse_locus(self, slocus, args, queue, lock=None ):

        '''This function takes as input a "superlocus" instance and the pipeline configuration.
        It also accepts as optional keywords a dictionary with the CDS information (derived from a bed12Parser)
        and a "lock" used for avoiding writing collisions during multithreading.
        The function splits the superlocus into its strand components and calls the relevant methods
        to define the loci. It also prints out the results to the requested output files.
        '''

        #Load the CDS information
        slocus.load_transcript_data()
        #Split the superlocus in the stranded components
        stranded_loci = sorted(list(slocus.split_strands()))
        #Define the loci
        for stranded_locus in stranded_loci:
            stranded_locus.define_loci()
    
        #Remove overlapping fragments.
        #This part should be rewritten in order to make it more flexible and powerful.
        for stranded_locus in stranded_loci:
            if args.remove_overlapping_fragments is True and len(stranded_loci)>1:
                for final_locus in stranded_locus.loci:
                    for other_superlocus in filter(lambda x: x!=stranded_locus, stranded_loci):
                        for other_final_locus in other_superlocus.loci:
                            if other_final_locus.other_is_fragment( final_locus ) is True and final_locus in stranded_locus.loci: 
                                stranded_locus.loci.remove(final_locus)
    #                             except ValueError as err:
    #                                 if "not in list" in err: pass
    #                                 else:
    #                                     raise ValueError(err)
    #                             finally:
    #                                 break
            self.printer_queue.put(stranded_locus)
        return
    
    
    def __call__(self):
        
        '''Method to call the pipeline.'''
        
        self.context=multiprocessing.get_context("fork") #@UndefinedVariable
        self.manager=ctx.Manager() # @UndefinedVariable
        self.printer_queue=self.manager.Queue() 
        self.lock=self.manager.RLock()
        self.pool=ctx.Pool(processes=self.json_conf["threads"]) # @UndefinedVariable
        self.printer_process=multiprocessing.context.Process(target=self.printer)
        self.printer_process.start()        
    