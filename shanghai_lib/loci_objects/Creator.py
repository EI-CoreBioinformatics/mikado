#!/usr/bin/env python3

#import sys
#from logging import Logger
import re
import multiprocessing
import shanghai_lib.loci_objects
import shanghai_lib.parsers
import csv
import os
import logging
import logging.handlers

class Creator:
    
    def __init__(self, json_conf):
        
        if type(json_conf) is str:
            assert os.path.exists(json_conf)
            json_conf = shanghai_lib.json_utils.to_json(json_conf) 
        else:
            assert type(json_conf) is dict
        
        self.json_conf = json_conf
        self.threads = self.json_conf["run_options"]["threads"]
        self.input_file = self.json_conf["input"]
        _ = self.define_input() #Check the input file
        self.sub_out = self.json_conf["subloci_out"]
        self.monolocus_out = self.json_conf["monoloci_out"]
        self.locus_out = self.json_conf["loci_out"]
        if self.locus_out is None:
            raise shanghai_lib.exceptions.InvalidJson("No output prefix specified for the final loci. Key: \"loci_out\"")
        
    def define_input(self):
        '''Function to check that the input file exists and is valid. It returns the parser.'''
        
        if self.input_file.endswith(".gtf"):
            parser=shanghai_lib.parsers.GTF.GTF
        else:
            parser=shanghai_lib.parsers.GFF.GFF3
            
        verified=False
        for row in parser(self.input_file):
            if row.header is False:
                verified=True
                break
        if verified is False:
            raise shanghai_lib.exceptions.InvalidJson("Invalid input file: {0}".format(self.input_file))

        return parser(self.input_file)

    def logging_utility(self):
        self.formatter = logging.Formatter("{asctime} - {levelname} - {lineno} - {funcName} - {name} - {message}",
                                           style="{"
                                           )
        self.logger=logging.getLogger("main_logger")
        if self.json_conf["log_settings"]["log"] is None  or self.json_conf["log_settings"]["log"]=="stream":
            self.log_handler = logging.StreamHandler()
        else:
            self.log_handler=logging.FileHandler(self.json_conf["log_settings"]["log"], 'w')
        self.log_level = self.json_conf["log_settings"]["log_level"] 
                
        self.log_handler.setFormatter(self.formatter)
        self.logger.setLevel(self.log_level)
        self.logger.addHandler(self.log_handler)
        self.logger.info("Begun analysis of {0}".format(self.input_file)  )
        while True:
            record = self.logging_queue.get()
            if record is None:
                break
            self.logger.handle(record)
            
        self.logger.info("Finished analysis of {0}".format(self.input_file)  )
        return

        
    def printer(self):

        '''Listener process that will print out the loci recovered by the analyse_locus function.'''
        

        score_keys = ["tid","parent","score"] + sorted(list(self.json_conf["scoring"].keys()))
        #Define mandatory output files        
        self.locus_metrics_file = re.sub("$",".metrics.tsv",  re.sub(".gff3$", "", self.locus_out  ))
        self.locus_scores_file = re.sub("$",".scores.tsv",  re.sub(".gff3$", "", self.locus_out  ))
        locus_metrics=csv.DictWriter(open(self.locus_metrics_file,'w'), shanghai_lib.loci_objects.superlocus.superlocus.available_metrics, delimiter="\t")
        locus_metrics.writeheader()
        locus_scores=csv.DictWriter(open(self.locus_scores_file,'a'), score_keys, delimiter="\t")
        locus_scores.writeheader()
        locus_out=open(self.locus_out,'w')
        print('##gff-version 3', file=locus_out)

        if self.sub_out is not None:
            self.sub_metrics_file=re.sub("$",".metrics.tsv",  re.sub(".gff3$", "", self.sub_out  ))
            self.sub_scores_file=re.sub("$",".scores.tsv",  re.sub(".gff3$", "", self.sub_out  ))
            sub_metrics=csv.DictWriter(open(self.sub_metrics_file,'w'), shanghai_lib.loci_objects.superlocus.superlocus.available_metrics, delimiter="\t")
            sub_metrics.writeheader()
            sub_scores=csv.DictWriter(open(self.sub_scores_file,'w'), score_keys, delimiter="\t")
            sub_scores.writeheader()
            sub_out=open(self.sub_out,'w')
            print('##gff-version 3', file=sub_out)
        else:
            self.sub_metrics=self.sub_scores=None
        
        if self.monolocus_out is not None:
            mono_out=open(self.monolocus_out,'w')
            print('##gff-version 3', file=mono_out)
        
        while True:
            stranded_locus=self.printer_queue.get()
            if stranded_locus=="EXIT":
                break #Poison pill - once we receive a "EXIT" signal, we exit
            if self.sub_out is not None: #Skip this section if no sub_out is defined
                sub_lines = stranded_locus.__str__(level="subloci", print_cds=not self.json_conf["run_options"]["exclude_cds"] )
                sub_metrics_rows = [x for x in stranded_locus.print_subloci_metrics()]
                sub_scores_rows = [x for x in stranded_locus.print_subloci_scores()]
                for row in sub_metrics_rows: sub_metrics.writerow(row)
                for row in sub_scores_rows: sub_scores.writerow(row)
                print(sub_lines, file=sub_out)
            if self.monolocus_out is not None:
                mono_lines = stranded_locus.__str__(level="monosubloci", print_cds=not self.json_conf["run_options"]["exclude_cds"])
                if mono_lines!='':
                    print(mono_lines, file=mono_out)
            locus_metrics_rows=[x for x in stranded_locus.print_monoholder_metrics()]
            locus_scores_rows = [x for x in stranded_locus.print_monoholder_scores()]
            locus_lines = stranded_locus.__str__(print_cds=not self.json_conf["run_options"]["exclude_cds"])
            for row in locus_metrics_rows: locus_metrics.writerow(row)
            for row in locus_scores_rows: locus_scores.writerow(row)
            
            if locus_lines!='':
                print(locus_lines, file=locus_out)
                    
        return
        
    def analyse_locus(self, slocus ):

        '''This function takes as input a "superlocus" instance and the pipeline configuration.
        It also accepts as optional keywords a dictionary with the CDS information (derived from a bed12Parser)
        and a "lock" used for avoiding writing collisions during multithreading.
        The function splits the superlocus into its strand components and calls the relevant methods
        to define the loci. It also prints out the results to the requested output files.
        '''
    
        #Define the logger
        handler = logging.handlers.QueueHandler(self.logging_queue)
        logger = logging.getLogger( "{chr}:{start}-{end}".format(chr=slocus.chrom, start=slocus.start, end=slocus.end) )
        logger.addHandler(handler)
        logger.setLevel(self.json_conf["log_settings"]["log_level"]) #We need to set this to the lowest possible level, otherwise we overwrite the global configuration
    
        #Load the CDS information
        logger.info("Logging transcript data")
        slocus.load_transcript_data()
        #Split the superlocus in the stranded components
        logger.info("Splitting by strand")
        stranded_loci = sorted(list(slocus.split_strands()))
        #Define the loci        
        logger.debug("Divided into {0} loci".format(len(stranded_loci)))
        
        logger.info("Defining loci")
        for stranded_locus in stranded_loci:
            try:
                stranded_locus.define_loci()
                logger.info("Defined loci for {0}:{1}-{2}, strand: {3}".format(stranded_locus.chrom,
                                                                               stranded_locus.start,
                                                                               stranded_locus.err,
                                                                               stranded_locus.strand))
            except Exception as err:
                logger.exception("Error in defining loci for {0}:{1}-{2}, strand: {3}".format(stranded_locus.chrom,
                                                                               stranded_locus.start,
                                                                               stranded_locus.err,
                                                                               stranded_locus.strand))
                logger.exception("Exception: {0}".format(err))
                stranded_loci.remove(stranded_locus)
                
        logger.info("Defined loci")
        #Remove overlapping fragments.
        #This part should be rewritten in order to make it more flexible and powerful.
        for stranded_locus in stranded_loci:
            if self.json_conf["run_options"]["remove_overlapping_fragments"] is True and len(stranded_loci)>1:
                for final_locus in stranded_locus.loci:
                    for other_superlocus in filter(lambda x: x!=stranded_locus, stranded_loci):
                        for other_final_locus in other_superlocus.loci:
                            if other_final_locus.other_is_fragment( final_locus ) is True and final_locus in stranded_locus.loci: 
                                stranded_locus.loci.remove(final_locus)
            if len(stranded_locus.transcripts)==0:
                continue
            else:
                try:
                    self.printer_queue.put(stranded_locus)
                    logger.info("Finished")
                except Exception as err:
                    logger.exception(err)
        
        #close up shop
        logger.removeHandler(handler)
        handler.close()
            
        return

    def __getstate__(self):
        state = self.__dict__.copy()
        if "pool" in state:
            del state['pool']
        if 'manager' in state:
            del state['manager']
        if 'log_process' in state:
            del state['log_process']
            
        if 'printer_process' in state:
            del state['printer_process']
            
        if 'log_process' in state:
            del state['log_process']
            
        return state
        

    def __call__(self):
        
        '''This method will activate the class and start the analysis of the input file.'''
        
        
        #NOTE: Pool, Process and Manager must NOT become instance attributes!
        #Otherwise it will raise all sorts of mistakes
        
        self.ctx = multiprocessing.get_context() # @UndefinedVariable 
        self.manager=self.ctx.Manager()
        self.printer_queue=self.manager.Queue()
        self.logging_queue = self.manager.Queue() #queue for logging
        
        self.lock=self.manager.RLock()
        self.pool=self.ctx.Pool(processes=self.threads) # @UndefinedVariable
        self.log_process = self.ctx.Process(target = self.logging_utility)
        self.log_process.start()
        
        self.printer_process=self.ctx.Process(target=self.printer) # @UndefinedVariable
        self.printer_process.start()
        
        currentLocus = None
        currentTranscript = None
        
        self.logger_queue_handler = logging.handlers.QueueHandler(self.logging_queue)
        self.queue_logger = logging.getLogger( "parser")
        self.queue_logger.addHandler(self.logger_queue_handler)
        self.queue_logger.setLevel(self.json_conf["log_settings"]["log_level"]) #We need to set this to the lowest possible level, otherwise we overwrite the global configuration
        
        jobs=[]
        for row in self.define_input():
            if row.is_exon is True:
                currentTranscript.addExon(row)
            elif row.is_transcript is True:
                if currentTranscript is not None:
                    if shanghai_lib.loci_objects.superlocus.superlocus.in_locus(currentLocus, currentTranscript) is True:
                        currentLocus.add_transcript_to_locus(currentTranscript, check_in_locus=False)
                        assert currentTranscript.id in currentLocus.transcripts
                    else:
                        if currentLocus is not None:
                            jobs.append(self.pool.apply_async(self.analyse_locus, args=(currentLocus,)))
                        currentLocus=shanghai_lib.loci_objects.superlocus.superlocus(currentTranscript, stranded=False, json_dict=self.json_conf)
                currentTranscript=shanghai_lib.loci_objects.transcript.transcript(row, source=self.json_conf["source"])
            else:
                continue
        
        if currentTranscript is not None:
            if shanghai_lib.loci_objects.superlocus.superlocus.in_locus(currentLocus, currentTranscript) is True:
                currentLocus.add_transcript_to_locus(currentTranscript)
            else:
                jobs.append(self.pool.apply_async(self.analyse_locus, args=(currentLocus,)))
                currentLocus=shanghai_lib.loci_objects.superlocus.superlocus(currentTranscript, stranded=False, json_dict=self.json_conf)
                
        if currentLocus is not None:
            jobs.append(self.pool.apply_async(self.analyse_locus, args=(currentLocus,)))
            
        list(map(lambda job: job.get(), jobs)) #Finish up all jobs
        self.pool.close()
        self.pool.join()
        
        self.printer_queue.put("EXIT")
        self.printer_process.join()
        self.logging_queue.put(None)
        self.log_process.join()
        return 0