#!/usr/bin/env python3

import re
import csv
import os
import logging
from logging import handlers as logging_handlers
import time

#SQLAlchemy/DB imports
from sqlalchemy.engine import create_engine
from sqlalchemy.orm.session import sessionmaker
import sqlalchemy
import sqlite3

#Shanghai imports
import mikado_lib.loci_objects
import mikado_lib.parsers
import mikado_lib.serializers.blast_utils
from mikado_lib.loci_objects.superlocus import superlocus
import multiprocessing
from multiprocessing.context import Process


#For profiling
#from memory_profiler import profile
# if "line_profiler" not in dir(): #@UndefinedVariable
#     def profile(function):
#         def inner(*args, **kwargs):
#             return function(*args, **kwargs)
#         return inner

class Creator:

    #@profile
    def __init__(self, json_conf: dict, commandline = ""):
        
        '''Constructor. It takes a single argument as input - the JSON/YAML configuration,
        prepared by the json_utils functions.
        Optional keyword: commandline (i.e. the commandline used to start the program)'''
        
        
        if type(json_conf) is str:
            assert os.path.exists(json_conf)
            json_conf = mikado_lib.json_utils.to_json(json_conf) 
        else:
            assert type(json_conf) is dict
        
        self.commandline = commandline
        self.json_conf = json_conf
        self.threads = self.json_conf["run_options"]["threads"]
        self.input_file = self.json_conf["input"]
        _ = self.define_input() #Check the input file
        self.sub_out = self.json_conf["subloci_out"]
        self.monolocus_out = self.json_conf["monoloci_out"]
        self.locus_out = self.json_conf["loci_out"]
        self.context = multiprocessing.get_context() #@UndefinedVariable
        self.manager = self.context.Manager() #@UndefinedVariable
        self.printer_queue = self.manager.Queue(-1)
        self.logging_queue = self.manager.Queue(-1) #queue for logging
      
        self.setup_logger()
        self.logger_queue_handler = logging_handlers.QueueHandler(self.logging_queue) # @UndefinedVariable
        self.queue_logger = logging.getLogger( "parser")
        self.queue_logger.addHandler(self.logger_queue_handler)
        self.queue_logger.setLevel(self.json_conf["log_settings"]["log_level"]) #We need to set this to the lowest possible level, otherwise we overwrite the global configuration
        self.queue_logger.propagate = False
        
        if self.locus_out is None:
            raise mikado_lib.exceptions.InvalidJson("No output prefix specified for the final loci. Key: \"loci_out\"")

    def define_input(self):
        '''Function to check that the input file exists and is valid. It returns the parser.'''
        
        if self.input_file.endswith(".gtf"):
            parser=mikado_lib.parsers.GTF.GTF
        else:
            parser=mikado_lib.parsers.GFF.GFF3
            
        verified=False
        for row in parser(self.input_file):
            if row.header is False:
                verified=True
                break
        if verified is False:
            raise mikado_lib.exceptions.InvalidJson("Invalid input file: {0}".format(self.input_file))

        return parser(self.input_file)

    def setup_logger(self):
        
        '''This function sets up the logger for the class. It creates the instance attribute "log_writer", which
        is itself a logging.handlers.QueueListener instance listening on the logging_queue instance attribute
        (which is a normal mp.Manager.Queue instance).'''
        
        self.formatter = logging.Formatter("{asctime} - {levelname} - {module}:{lineno} - {funcName} - {name} - {message}",
                                           style="{"
                                           )
        self.main_logger=logging.getLogger("main_logger")
        self.logger=logging.getLogger("listener")
        self.logger.propagate=False
        if self.json_conf["log_settings"]["log"] is None  or self.json_conf["log_settings"]["log"]=="stream":
            self.log_handler = logging.StreamHandler()
        else:
            self.log_handler=logging.FileHandler(self.json_conf["log_settings"]["log"], 'w')
        self.log_level = self.json_conf["log_settings"]["log_level"] # For the main logger I want to keep it at the "INFO" level
                
        self.log_handler.setFormatter(self.formatter)
        self.logger.setLevel(self.log_level)
        self.logger.addHandler(self.log_handler)
        
        self.main_logger.setLevel(logging.INFO)
        self.main_logger.addHandler(self.log_handler)
        
        self.main_logger.info("Begun analysis of {0}".format(self.input_file)  )
        if self.commandline!='':
            self.main_logger.info("Command line: {0}".format(self.commandline))
        else:
            self.main_logger.info("Analysis launched directly, without using the launch script.")
        
        if self.json_conf["chimera_split"]["blast_check"] is True:
            engine=create_engine( "{0}://".format(self.json_conf["dbtype"]), creator=self.db_connection)
            Session = sessionmaker()
            Session.configure(bind=engine)
            session=Session()

            evalue=self.json_conf["chimera_split"]["blast_params"]["evalue"]
            queries_with_hits = session.query(mikado_lib.serializers.blast_utils.Hit.query_id ).filter(
                                                                                                 mikado_lib.serializers.blast_utils.Hit.evalue<=evalue,
                                                                                            ).distinct().count()
            total_queries = session.query(mikado_lib.serializers.blast_utils.Query).count()
            self.main_logger.info("Queries with at least one hit at evalue<={0}: {1} out of {2} ({3}%)".format(
                                                                                                               evalue,
                                                                                                               queries_with_hits,
                                                                                                               total_queries,
                                                                                                               round(100*queries_with_hits/total_queries,2)
                                                                                                               ))
            session.close()
        
        self.log_writer = logging_handlers.QueueListener(self.logging_queue, self.logger)
        self.log_writer.start()
        
        return


    def printer(self):

        '''Listener process that will print out the loci recovered by the analyse_locus function.'''
        

        handler = logging_handlers.QueueHandler(self.logging_queue) # @UndefinedVariable
        logger = logging.getLogger( "queue_listener")
        logger.propagate=False
        logger.addHandler(handler)
        logger.setLevel(self.json_conf["log_settings"]["log_level"])

        score_keys = ["tid","parent","score"] + sorted(list(self.json_conf["scoring"].keys()))
        #Define mandatory output files        
        self.locus_metrics_file = re.sub("$",".metrics.tsv",  re.sub(".gff.?$", "", self.locus_out  ))
        self.locus_scores_file = re.sub("$",".scores.tsv",  re.sub(".gff.?$", "", self.locus_out  ))
        locus_metrics=csv.DictWriter(open(self.locus_metrics_file,'w'), mikado_lib.loci_objects.superlocus.superlocus.available_metrics, delimiter="\t")
        locus_metrics.writeheader()
        locus_scores=csv.DictWriter(open(self.locus_scores_file,'w'), score_keys, delimiter="\t")
        locus_scores.writeheader()
        locus_out=open(self.locus_out,'w')
        print('##gff-version 3', file=locus_out)

        if self.sub_out is not None:
            self.sub_metrics_file=re.sub("$",".metrics.tsv",  re.sub(".gff.?$", "", self.sub_out  ))
            self.sub_scores_file=re.sub("$",".scores.tsv",  re.sub(".gff.?$", "", self.sub_out  ))
            sub_metrics=csv.DictWriter(open(self.sub_metrics_file,'w'), mikado_lib.loci_objects.superlocus.superlocus.available_metrics, delimiter="\t")
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
#             except multiprocessing.managers.RemoteError as err:
#                 logger.exception(err)
#                 continue
            
            if stranded_locus=="EXIT":
                return #Poison pill - once we receive a "EXIT" signal, we exit
            logger.debug("Received {0}".format(stranded_locus.id))
            stranded_locus.set_logger(logger)
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
            self.printer_queue.task_done()
        return
    
    def db_connection(self):
        '''Creator function for the database connection. It necessitates the following information from
        the json_conf dictionary:
        
        - dbtype (one of sqlite, mysql, postgresql)
        - db (name of the database file, for sqlite, otherwise name of the database)
        
        If the database is MySQL/PostGreSQL, the method also requires:
        
        - dbuser
        - dbhost
        - dbpasswd
        - dbport
        
        These are controlled and added automatically by the json_utils functions.
        
        '''
        
        if self.json_conf["dbtype"] == "sqlite":
            return sqlite3.connect(database= self.json_conf["db"]) #@UndefinedVariable
        elif self.json_conf["dbtype"] == "mysql":
            import MySQLdb
            return MySQLdb.connect(host = self.json_conf["dbhost"],
                                   user = self.json_conf["dbuser"],
                                   passwd = self.json_conf["dbpasswd"],
                                   db = self.json_conf["db"],
                                   port = self.json_conf["dbport"]
                                    )
        elif self.json_conf["dbtype"] == "postgresql":
            import psycopg2
            return psycopg2.connect(
                                    host = self.json_conf["dbhost"],
                                   user = self.json_conf["dbuser"],
                                   password = self.json_conf["dbpasswd"],
                                   database = self.json_conf["db"],
                                   port = self.json_conf["dbport"]
                                   ) 
        
#             import postgresql
#             return postgresql.open("pq://{user}:{passwd}@{host}/{db}".format(
#                                                                                      host = self.json_conf["dbhost"],
#                                                                                      user = self.json_conf["dbuser"],
#                                                                                      passwd = self.json_conf["dbpasswd"],
#                                                                                      db = self.json_conf["db"]                         
#                                                                                      )  )
            
    
    
    #@profile
    def analyse_locus(self, slocus: superlocus ) -> [superlocus]:

        '''This function takes as input a "superlocus" instance and the pipeline configuration.
        It also accepts as optional keywords a dictionary with the CDS information (derived from a bed12Parser)
        and a "lock" used for avoiding writing collisions during multithreading.
        The function splits the superlocus into its strand components and calls the relevant methods
        to define the loci.
        When it is finished, it transmits the superloci to the printer function.
        '''
    
        #Define the logger
        if slocus is None: return
        handler = logging_handlers.QueueHandler(self.logging_queue) # @UndefinedVariable
        logger = logging.getLogger( "{chr}:{start}-{end}".format(chr=slocus.chrom, start=slocus.start, end=slocus.end) )
        logger.addHandler(handler)
        logger.setLevel(self.json_conf["log_settings"]["log_level"]) #We need to set this to the lowest possible level, otherwise we overwrite the global configuration
        logger.propagate = False
    
        #Load the CDS information
        logger.info("Started with {0}".format(slocus.id))
        logger.debug("Loading transcript data")
        slocus.set_logger(logger)
        slocus.load_all_transcript_data(pool = self.connection_pool)
        #Split the superlocus in the stranded components
        logger.debug("Splitting by strand")
        stranded_loci = sorted(list(slocus.split_strands()))
        #Define the loci        
        logger.debug("Divided into {0} loci".format(len(stranded_loci)))
        
        logger.info("Defining loci")
        for stranded_locus in stranded_loci:
            try:
                stranded_locus.define_loci()
                logger.debug("Defined loci for {0}:{1}-{2}, strand: {3}".format(stranded_locus.chrom,
                                                                               stranded_locus.start,
                                                                               stranded_locus.end,
                                                                               stranded_locus.strand))
            except Exception as err:
                logger.exception("Error in defining loci for {0}:{1}-{2}, strand: {3}".format(stranded_locus.chrom,
                                                                               stranded_locus.start,
                                                                               stranded_locus.end,
                                                                               stranded_locus.strand))
                logger.exception("Exception: {0}".format(err))
                stranded_loci.remove(stranded_locus)
                
        logger.debug("Defined loci")
        
        #Remove overlapping fragments.
        loci_to_check = {True: set(), False:set()}
        for stranded_locus in stranded_loci:
            for _,locus_instance in stranded_locus.loci.items():
                loci_to_check[locus_instance.monoexonic].add(locus_instance)
        
        
        for stranded_locus in stranded_loci:
            for locus_id,locus_instance in stranded_locus.loci.items():
                if locus_instance in loci_to_check[True]:
                    logger.debug("Checking if {0} is a fragment".format( _ ))
                    for other_locus in loci_to_check[False]:
                        if other_locus.other_is_fragment( locus_instance, 
                                                          minimal_cds_length=self.json_conf["run_options"]["fragments_maximal_cds"] ) is True:
                            if self.json_conf["run_options"]["remove_overlapping_fragments"] is False:
                                stranded_locus.loci[locus_id].is_fragment=True
                            else:
                                del stranded_locus.loci[locus_id]
                        break
                    
            putter_counter = 0
            while True:
                try:
                    self.printer_queue.put(stranded_locus)
                    logger.debug("Finished for {0}:{1}-{2}, strand: {3}".format(stranded_locus.chrom,
                                                                               stranded_locus.start,
                                                                               stranded_locus.end,
                                                                               stranded_locus.strand))
                    break
                except Exception as err:
                    if putter_counter<10:
                        putter_counter+=1
                        time.sleep(0.0001)
                    else:
                        logger.exception("Error in reporting for {0}:{1}-{2}, strand: {3}".format(stranded_locus.chrom,
                                                                                                  stranded_locus.start,
                                                                                                  stranded_locus.end,
                                                                                                  stranded_locus.strand))
                        logger.exception(err)
                        break
        
        #close up shop
        logger.info("Finished with {0}".format(slocus.id))
        logger.removeHandler(handler)
        handler.close()
        return


    def __getstate__(self):
        self.not_pickable = ["queue_logger", "manager", "printer_process", "log_process", "pool", "main_logger", "log_handler", "log_writer", "logger"]
        state = self.__dict__.copy()
        for not_pickable in self.not_pickable:
            if not_pickable in state:
                del state[not_pickable]            
             
        return state
  
    #@profile
    def __call__(self):
        
        '''This method will activate the class and start the analysis of the input file.'''
        
        
        #NOTE: Pool, Process and Manager must NOT become instance attributes!
        #Otherwise it will raise all sorts of mistakes
        
        self.connection_pool = sqlalchemy.pool.QueuePool( self.db_connection, pool_size=self.threads, max_overflow=self.threads*2 )
        self.printer_process=Process(target=self.printer) # @UndefinedVariable
        self.printer_process.start()
        
        currentLocus = None
        currentTranscript = None
        
        jobs = []
        
        self.logger.debug("Source: {0}".format(self.json_conf["source"]))        
        for row in self.define_input():
            if row.is_exon is True:
                currentTranscript.addExon(row)
            elif row.is_transcript is True:
                if currentTranscript is not None:
                    if mikado_lib.loci_objects.superlocus.superlocus.in_locus(currentLocus, currentTranscript) is True:
                        currentLocus.add_transcript_to_locus(currentTranscript, check_in_locus=False)
                        assert currentTranscript.id in currentLocus.transcripts
                    else:
#                         self.analyse_locus(currentLocus)
                        while len(jobs)>=self.threads:
                            for job in jobs:
                                if job.is_alive() is False:
                                    jobs.remove(job)
                            time.sleep(0.1)

                        job = Process(target=self.analyse_locus, args=(currentLocus,))
                        job.start()
                        jobs.append(job)

                        currentLocus=mikado_lib.loci_objects.superlocus.superlocus(currentTranscript, stranded=False, json_dict=self.json_conf)
                currentTranscript=mikado_lib.loci_objects.transcript.transcript(row, source=self.json_conf["source"])
            else:
                continue
        
        if currentTranscript is not None:
            if mikado_lib.loci_objects.superlocus.superlocus.in_locus(currentLocus, currentTranscript) is True:
                currentLocus.add_transcript_to_locus(currentTranscript)
            else:
#                 self.analyse_locus(currentLocus)
                while len(jobs)>=self.threads+1:
                    for job in jobs:
                        if job.is_alive() is False:
                            jobs.remove(job)
                    time.sleep(0.1)

                job = Process(target=self.analyse_locus, args=(currentLocus,))
                job.start()
                jobs.append(job)

                

                currentLocus=mikado_lib.loci_objects.superlocus.superlocus(currentTranscript, stranded=False, json_dict=self.json_conf)
                
        if currentLocus is not None:
#             self.analyse_locus(currentLocus)
            while len(jobs)>=self.threads+1:
                for job in jobs:
                    if job.is_alive() is False:
                        jobs.remove(job)
                time.sleep(0.1)

                job = Process(target=self.analyse_locus, args=(currentLocus,))
                job.start()
                jobs.append(job)
                
            job = Process(target=self.analyse_locus, args=(currentLocus,))
            job.start()
            jobs.append(job)
            
        for job in jobs:
            if job.is_alive() is True:
                job.join()

        self.printer_queue.join()
        self.printer_queue.put("EXIT")
        #The printing process must be started AFTER we have put the stopping signal  into the queue
        self.printer_process.join()
        self.main_logger.info("Finished analysis of {0}".format(self.input_file)  )
        self.log_writer.stop()
        return 0