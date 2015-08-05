#!/usr/bin/env python3
# coding=utf-8

"""
This module defines the Creator class, which is the main workhorse for Mikado pick.
"""

import re
import csv
import os
import logging
from logging import handlers as logging_handlers
import collections
import functools

# SQLAlchemy/DB imports
from sqlalchemy.engine import create_engine
from sqlalchemy.orm.session import sessionmaker
import sqlalchemy.pool
import sqlalchemy
import sqlite3

# Shanghai imports
import mikado_lib.loci_objects
import mikado_lib.parsers
import mikado_lib.serializers.blast_utils
from mikado_lib.loci_objects.superlocus import Superlocus
import multiprocessing
import multiprocessing.managers
from multiprocessing.context import Process

# For profiling
# from memory_profiler import profile
# if "line_profiler" not in dir():  # @UndefinedVariable
#     def profile(function):
#         """
#         Mock wrapper to make the program function also without memory_profile/cProfile
#         enabled.
#         :param function: the function to be wrapped
#         """
#         def inner(*args, **kwargs):
#             """Inner function of the wrapper
#             :param args: positional arguments
#             :param kwargs: keyword arguments
#             """
#             return function(*args, **kwargs)
#         return inner
#


# @profile
def connector(json_conf):
    """Creator function for the database connection. It necessitates the following information from
    the json_conf dictionary:

    - dbtype (one of sqlite, mysql, postgresql)
    - db (name of the database file, for sqlite, otherwise name of the database)

    If the database is MySQL/PostGreSQL, the method also requires:

    - dbuser
    - dbhost
    - dbpasswd
    - dbport

    These are controlled and added automatically by the json_utils functions.

    :param json_conf: configuration dictionary
    :type json_conf: dict

    :rtype : MySQLdb.connect | sqlite3.connect | psycopg2.connect

    """

    if json_conf["dbtype"] == "sqlite":
        return sqlite3.connect(database=json_conf["db"])  # @UndefinedVariable
    elif json_conf["dbtype"] == "mysql":
        import MySQLdb
        return MySQLdb.connect(host=json_conf["dbhost"],
                               user=json_conf["dbuser"],
                               passwd=json_conf["dbpasswd"],
                               db=json_conf["db"],
                               port=json_conf["dbport"]
                               )
    elif json_conf["dbtype"] == "postgresql":
        import psycopg2
        return psycopg2.connect(
            host=json_conf["dbhost"],
            user=json_conf["dbuser"],
            password=json_conf["dbpasswd"],
            database=json_conf["db"],
            port=json_conf["dbport"]
        )


# @profile
def analyse_locus(slocus: Superlocus,
                  json_conf: dict,
                  printer_queue: multiprocessing.managers.AutoProxy,
                  logging_queue: multiprocessing.managers.AutoProxy,
                  connection_pool,
                  data_dict
                  ) -> [Superlocus]:

    """
    :param slocus: a superlocus instance
    :type slocus: mikado_lib.loci_objects.superlocus.Superlocus

    :param json_conf: the configuration dictionary
    :type json_conf: dict

    :param logging_queue: the logging queue
    :type logging_queue: multiprocessing.managers.AutoProxy

    :param printer_queue: the printing queue
    :type printer_queue: multiprocessing.managers.AutoProxy

    :param connection_pool: the connection pool

    :param data_dict: a dictionary of preloaded data
    :type data_dict: dict

    This function takes as input a "superlocus" instance and the pipeline configuration.
    It also accepts as optional keywords a dictionary with the CDS information (derived from a Bed12Parser)
    and a "lock" used for avoiding writing collisions during multithreading.
    The function splits the superlocus into its strand components and calls the relevant methods
    to define the loci.
    When it is finished, it transmits the superloci to the printer function.
    """

    # Define the logger
    if slocus is None:
        return

    if connection_pool is None:
        # Create the connection pool from scratch
        db_connection = functools.partial(connector, json_conf)
        connection_pool = sqlalchemy.pool.QueuePool(db_connection, pool_size=1, max_overflow=2)

    handler = logging_handlers.QueueHandler(logging_queue)  # @UndefinedVariable
    logger = logging.getLogger("{chr}:{start}-{end}".format(chr=slocus.chrom, start=slocus.start, end=slocus.end))
    logger.addHandler(handler)
    # We need to set this to the lowest possible level, otherwise we overwrite the global configuration
    logger.setLevel(json_conf["log_settings"]["log_level"])
    logger.propagate = False

    # Load the CDS information
    logger.info("Started with {0}".format(slocus.id))
    logger.debug("Loading transcript data")
    slocus.logger = logger
    slocus.load_all_transcript_data(pool=connection_pool, data_dict=data_dict)
    # Split the superlocus in the stranded components
    logger.debug("Splitting by strand")
    stranded_loci = sorted(list(slocus.split_strands()))
    # Define the loci
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

    # Remove overlapping fragments.
    loci_to_check = {True: set(), False: set()}
    for stranded_locus in stranded_loci:
        for _, locus_instance in stranded_locus.loci.items():
            locus_instance.logger = logger
            loci_to_check[locus_instance.monoexonic].add(locus_instance)

    for stranded_locus in stranded_loci:
        for locus_id, locus_instance in stranded_locus.loci.items():
            if locus_instance in loci_to_check[True]:
                logger.debug("Checking if {0} is a fragment".format(locus_instance.id))
                for other_locus in loci_to_check[False]:
                    if other_locus.other_is_fragment(locus_instance,
                                                     minimal_cds_length=json_conf["run_options"][
                                                         "fragments_maximal_cds"]) is True:
                        if json_conf["run_options"]["remove_overlapping_fragments"] is False:
                            stranded_locus.loci[locus_id].is_fragment = True
                        else:
                            del stranded_locus.loci[locus_id]
                    break

        # putter_counter = 0
        printer_queue.put(stranded_locus)
        logger.debug("Finished for {0}:{1}-{2}, strand: {3}".format(stranded_locus.chrom,
                                                                    stranded_locus.start,
                                                                    stranded_locus.end,
                                                                    stranded_locus.strand))
        # while True:
        #     try:
        #         self.printer_queue.put(stranded_locus)
        #         logger.debug("Finished for {0}:{1}-{2}, strand: {3}".format(stranded_locus.chrom,
        #                                                                     stranded_locus.start,
        #                                                                     stranded_locus.end,
        #                                                                     stranded_locus.strand))
        #         break
        #     except Exception as err:
        #         if putter_counter < 10:
        #             putter_counter += 1
        #             time.sleep(0.0001)
        #         else:
        #             message = "Error in reporting for {0}:{1}-{2}".format(stranded_locus.chrom,
        #                                                                   stranded_locus.start,
        #                                                                   stranded_locus.end)
        #             message += ", strand {0}".format(stranded_locus.strand)
        #             logger.exception(message)
        #             logger.exception(err)
        #             break

    # close up shop
    logger.info("Finished with {0}".format(slocus.id))
    logger.removeHandler(handler)
    handler.close()
    return


class Creator:

    """
    This class is used to launch the main Mikado pipeline. Its purpose is to parse
    an input sorted annotation file, locate the loci, and perform the selection analysis
    using the parameters provided in the input configuration file.
    """

    # @profile
    def __init__(self, json_conf, commandline=""):

        """Constructor. It takes a single argument as input - the JSON/YAML configuration,
        prepared by the json_utils functions.

        :param json_conf: Either a configuration dictionary or the configuration file.
        :type json_conf: str,dict

        :param commandline: optional, the commandline used to start the program
        :type commandline: str
        """

        # Mock variables
        self.formatter = self.main_logger = self.log_writer = self.log_handler = self.logger = None
        self.log_level = "WARN"

        # Now we start the real work
        if type(json_conf) is str:
            assert os.path.exists(json_conf)
            json_conf = mikado_lib.json_utils.to_json(json_conf)
        else:
            assert type(json_conf) is dict

        self.commandline = commandline
        self.json_conf = json_conf
        self.threads = self.json_conf["run_options"]["threads"]
        self.input_file = self.json_conf["input"]
        _ = self.define_input()  # Check the input file
        self.sub_out = self.json_conf["subloci_out"]
        self.monolocus_out = self.json_conf["monoloci_out"]
        self.locus_out = self.json_conf["loci_out"]
        self.context = multiprocessing.get_context()
        self.manager = self.context.Manager()
        self.printer_queue = self.manager.Queue(-1)
        self.logging_queue = self.manager.Queue(-1)

        # self.printer_queue = self.manager.Queue(-1)
        # self.logging_queue = self.manager.Queue(-1)  # queue for logging

        self.db_connection = functools.partial(connector, self.json_conf)
        self.setup_logger()
        self.logger_queue_handler = logging_handlers.QueueHandler(self.logging_queue)
        self.queue_logger = logging.getLogger("parser")
        self.queue_logger.addHandler(self.logger_queue_handler)
        # We need to set this to the lowest possible level, otherwise we overwrite the global configuration
        self.queue_logger.setLevel(self.json_conf["log_settings"]["log_level"])
        self.queue_logger.propagate = False
        if self.json_conf["single_thread"] is True:
            # Reset threads to 1
            self.logger.warning("Reset number of threads to 1 as requested")
            self.threads = self.json_conf["run_options"]["threads"] = 1

        if self.locus_out is None:
            raise mikado_lib.exceptions.InvalidJson("No output prefix specified for the final loci. Key: \"loci_out\"")

    def define_input(self):
        """Function to check that the input file exists and is valid. It returns the parser."""

        if self.input_file.endswith(".gtf"):
            parser = mikado_lib.parsers.GTF.GTF
        else:
            parser = mikado_lib.parsers.GFF.GFF3

        verified = False
        for row in parser(self.input_file):
            if row.header is False:
                verified = True
                break
        if verified is False:
            raise mikado_lib.exceptions.InvalidJson("Invalid input file: {0}".format(self.input_file))

        return parser(self.input_file)

    def setup_logger(self):

        """This function sets up the logger for the class. It creates the instance attribute "log_writer", which
        is itself a logging.handlers.QueueListener instance listening on the logging_queue instance attribute
        (which is a normal mp.Manager.Queue instance)."""

        self.formatter = logging.Formatter(
            "{asctime} - {levelname} - {module}:{lineno} - {funcName} - {name} - {message}",
            style="{"
            )
        self.main_logger = logging.getLogger("main_logger")
        self.logger = logging.getLogger("listener")
        self.logger.propagate = False
        if self.json_conf["log_settings"]["log"] is None or self.json_conf["log_settings"]["log"] == "stream":
            self.log_handler = logging.StreamHandler()
        else:
            self.log_handler = logging.FileHandler(self.json_conf["log_settings"]["log"], 'w')
        self.log_level = self.json_conf["log_settings"][
            "log_level"]  # For the main logger I want to keep it at the "INFO" level

        self.log_handler.setFormatter(self.formatter)
        self.logger.setLevel(self.log_level)
        self.logger.addHandler(self.log_handler)

        if self.log_level == "DEBUG":
            self.main_logger.setLevel(logging.DEBUG)
        else:
            self.main_logger.setLevel(logging.INFO)
        self.main_logger.addHandler(self.log_handler)

        self.main_logger.info("Begun analysis of {0}".format(self.input_file))
        if self.commandline != '':
            self.main_logger.info("Command line: {0}".format(self.commandline))
        else:
            self.main_logger.info("Analysis launched directly, without using the launch script.")

        if self.json_conf["chimera_split"]["blast_check"] is True:
            engine = create_engine("{0}://".format(self.json_conf["dbtype"]),
                                   creator=self.db_connection)
            smaker = sessionmaker()
            smaker.configure(bind=engine)
            session = smaker()

            evalue = self.json_conf["chimera_split"]["blast_params"]["evalue"]
            queries_with_hits = session.query(mikado_lib.serializers.blast_utils.Hit.query_id).filter(
                mikado_lib.serializers.blast_utils.Hit.evalue <= evalue,
            ).distinct().count()
            total_queries = session.query(mikado_lib.serializers.blast_utils.Query).count()
            self.main_logger.debug("Queries with at least one hit at evalue<={0}: {1} out of {2} ({3}%)".format(
                evalue,
                queries_with_hits,
                total_queries,
                round(100 * queries_with_hits / total_queries, 2)
            ))
            session.close()

        self.log_writer = logging_handlers.QueueListener(self.logging_queue, self.logger)
        self.log_writer.start()

        return

    def printer(self):

        """Listener process that will print out the loci recovered by the analyse_locus function."""

        handler = logging_handlers.QueueHandler(self.logging_queue)  # @UndefinedVariable
        logger = logging.getLogger("queue_listener")
        logger.propagate = False
        logger.addHandler(handler)
        logger.setLevel(self.json_conf["log_settings"]["log_level"])

        score_keys = ["tid", "parent", "score"] + sorted(list(self.json_conf["scoring"].keys()))
        # Define mandatory output files
        locus_metrics_file = re.sub("$", ".metrics.tsv", re.sub(".gff.?$", "", self.locus_out))
        locus_scores_file = re.sub("$", ".scores.tsv", re.sub(".gff.?$", "", self.locus_out))
        locus_metrics = csv.DictWriter(open(locus_metrics_file, 'w'),
                                       mikado_lib.loci_objects.superlocus.Superlocus.available_metrics, delimiter="\t")
        locus_metrics.writeheader()
        locus_scores = csv.DictWriter(open(locus_scores_file, 'w'), score_keys, delimiter="\t")
        locus_scores.writeheader()
        locus_out = open(self.locus_out, 'w')
        print('##gff-version 3', file=locus_out)

        if self.sub_out is not None:
            sub_metrics_file = re.sub("$", ".metrics.tsv", re.sub(".gff.?$", "", self.sub_out))
            sub_scores_file = re.sub("$", ".scores.tsv", re.sub(".gff.?$", "", self.sub_out))
            sub_metrics = csv.DictWriter(open(sub_metrics_file, 'w'),
                                         mikado_lib.loci_objects.superlocus.Superlocus.available_metrics,
                                         delimiter="\t")
            sub_metrics.writeheader()
            sub_scores = csv.DictWriter(open(sub_scores_file, 'w'), score_keys, delimiter="\t")
            sub_scores.writeheader()
            sub_out = open(self.sub_out, 'w')
            print('##gff-version 3', file=sub_out)
        else:
            sub_metrics = sub_scores = sub_out = None

        if self.monolocus_out is not None:
            mono_out = open(self.monolocus_out, 'w')
            print('##gff-version 3', file=mono_out)
        else:
            mono_out = None

        while True:
            stranded_locus = self.printer_queue.get()

            if stranded_locus == "EXIT":
                return  # Poison pill - once we receive a "EXIT" signal, we exit
            logger.debug("Received {0}".format(stranded_locus.id))
            stranded_locus.logger = logger
            if self.sub_out is not None:  # Skip this section if no sub_out is defined
                sub_lines = stranded_locus.__str__(level="subloci",
                                                   print_cds=not self.json_conf["run_options"]["exclude_cds"])
                sub_metrics_rows = [x for x in stranded_locus.print_subloci_metrics()]
                sub_scores_rows = [x for x in stranded_locus.print_subloci_scores()]
                for row in sub_metrics_rows:
                    sub_metrics.writerow(row)
                for row in sub_scores_rows:
                    sub_scores.writerow(row)
                print(sub_lines, file=sub_out)
            if self.monolocus_out is not None:
                mono_lines = stranded_locus.__str__(level="monosubloci",
                                                    print_cds=not self.json_conf["run_options"]["exclude_cds"])
                if mono_lines != '':
                    print(mono_lines, file=mono_out)
            locus_metrics_rows = [x for x in stranded_locus.print_monoholder_metrics()]
            locus_scores_rows = [x for x in stranded_locus.print_monoholder_scores()]
            locus_lines = stranded_locus.__str__(print_cds=not self.json_conf["run_options"]["exclude_cds"])
            for row in locus_metrics_rows:
                locus_metrics.writerow(row)
            for row in locus_scores_rows:
                locus_scores.writerow(row)

            if locus_lines != '':
                print(locus_lines, file=locus_out)
            self.printer_queue.task_done()
        return

    def __getstate__(self):
        self.not_pickable = ["queue_logger", "manager", "printer_process", "log_process", "pool", "main_logger",
                             "log_handler", "log_writer", "logger"]
        state = self.__dict__.copy()
        for not_pickable in self.not_pickable:
            if not_pickable in state:
                del state[not_pickable]

        return state

    # @profile
    def __call__(self):

        """This method will activate the class and start the analysis of the input file."""

        # NOTE: Pool, Process and Manager must NOT become instance attributes!
        # Otherwise it will raise all sorts of mistakes

        self.connection_pool = sqlalchemy.pool.QueuePool(self.db_connection, pool_size=self.threads,
                                                         max_overflow=self.threads * 2)
        self.printer_process = Process(target=self.printer)
        self.printer_process.start()

        current_locus = None
        current_transcript = None

        jobs = []

        data_dict = None
        if self.json_conf["run_options"]["preload"] is True:
            self.main_logger.info("Starting to preload the database into memory")

            data_dict = dict()
            engine = create_engine("{0}://".format(self.json_conf["dbtype"]),
                                   creator=self.db_connection)
            session = sqlalchemy.orm.sessionmaker(bind=engine)()

            data_dict["junctions"] = self.manager.dict()
            for x in session.query(mikado_lib.serializers.junction.Junction):
                data_dict["junctions"][(x.chrom, x.junctionStart, x.junctionEnd, x.strand)] = None

            self.main_logger.info("{0} junctions loaded".format(len(data_dict["junctions"])))
            queries = dict((x.query_id, x) for x in engine.execute("select * from query"))

            # Then load ORFs
            orfs = collections.defaultdict(list)

            for x in engine.execute("select * from orf"):
                query_name = queries[x.query_id].query_name
                orfs[query_name].append(mikado_lib.serializers.orf.Orf.as_bed12_static(x, query_name))

            data_dict['orf'] = self.manager.dict()
            for key in orfs:
                data_dict["orf"][key] = orfs[key]
            del orfs

            self.main_logger.info("{0} ORFs loaded".format(len(data_dict["orf"])))

            # Finally load BLAST

            data_dict["hit"] = self.manager.dict()

            if self.json_conf["chimera_split"]["execute"] is True and \
                    self.json_conf["chimera_split"]["blast_check"] is True:
                hsps = dict()
                for hsp in engine.execute("select * from hsp where hsp_evalue <= {0}".format(
                    self.json_conf["chimera_split"]["blast_params"]["hsp_evalue"]
                )).fetchall():
                    if hsp.query_id not in hsps:
                        hsps[hsp.query_id] = collections.defaultdict(list)
                    hsps[hsp.query_id][hsp.target_id].append(hsp)

                self.main_logger.info("{0} HSPs prepared".format(len(hsps)))

                targets = dict((x.target_id, x) for x in engine.execute("select * from target"))

                data_dict["hit"] = self.manager.dict()

                hit_counter = 0
                hits = engine.execute("select * from hit where evalue <= {0} order by query_id,evalue;".format(
                    self.json_conf["chimera_split"]["blast_params"]["evalue"]
                ))

                # self.main_logger.info("{0} BLAST hits to analyse".format(hits))
                current_counter = 0
                current_hit = None
                hits_dict = collections.defaultdict(list)

                for hit in hits:
                    if current_hit != hit.query_id:
                        current_hit = hit.query_id
                        current_counter = 0

                    current_counter += 1
                    if current_counter > self.json_conf["chimera_split"]["blast_params"]["max_target_seqs"]:
                        continue
                    my_query = queries[hit.query_id]
                    my_target = targets[hit.target_id]

                    # We HAVE to use the += approach because extend/append
                    # leave the original list empty
                    hits_dict[my_query.query_name].append(
                        mikado_lib.serializers.blast_utils.Hit.as_full_dict_static(
                            hit,
                            hsps[hit.query_id][hit.target_id],
                            my_query,
                            my_target
                        )
                    )
                    hit_counter += 1
                    if hit_counter >= 2*10**4 and hit_counter % (2*10**4) == 0:
                        self.main_logger.info("Loaded {0} BLAST hits in database".format(hit_counter))

                for key in hits_dict:
                    data_dict["hit"][key] = hits_dict[key]

                del hsps
                del hits_dict
                assert len(data_dict["hit"]) <= len(queries)
                self.main_logger.info("{0} BLAST hits loaded for {1} queries".format(
                    hit_counter,
                    len(data_dict["hit"])
                ))
                self.main_logger.debug("{0}".format(", ".join([str(x) for x in list(data_dict["hit"].keys())[:10]])))
            else:
                self.main_logger.info("Skipping BLAST loading")

        self.main_logger.info("Finished to preload the database into memory")

        pool = multiprocessing.Pool(processes=self.threads)

        self.logger.debug("Source: {0}".format(self.json_conf["source"]))
        for row in self.define_input():
            if row.is_exon is True:
                current_transcript.add_exon(row)
            elif row.is_transcript is True:
                if current_transcript is not None:
                    if mikado_lib.loci_objects.superlocus.Superlocus.in_locus(current_locus, current_transcript) \
                            is True:
                        current_locus.add_transcript_to_locus(current_transcript, check_in_locus=False)
                        assert current_transcript.id in current_locus.transcripts
                    else:
                        if self.json_conf["single_thread"] is True:
                            analyse_locus(current_locus,
                                          self.json_conf,
                                          self.printer_queue,
                                          self.logging_queue,
                                          None,
                                          data_dict
                                          # self.connection_pool
                                          )
                        else:
                            # while len(jobs) >= self.threads:
                            #     for job in jobs:
                            #         if job.is_alive() is False:
                            #             jobs.remove(job)
                            #
                            # job = Process(target=analyse_locus, args=(current_locus,
                            #                                           self.json_conf,
                            #                                           self.printer_queue,
                            #                                           self.logging_queue,
                            #                                           self.connection_pool
                            #                                           ))
                            # job.start()
                            # jobs.append(job)
                            jobs.append(pool.apply_async(analyse_locus, args=(current_locus,
                                                                              self.json_conf,
                                                                              self.printer_queue,
                                                                              self.logging_queue,
                                                                              None,
                                                                              data_dict
                                                                              # self.connection_pool
                                                                              )))

                        current_locus = mikado_lib.loci_objects.superlocus.Superlocus(current_transcript,
                                                                                      stranded=False,
                                                                                      json_dict=self.json_conf)
                current_transcript = mikado_lib.loci_objects.transcript.Transcript(row, source=self.json_conf["source"])
            else:
                continue

        if current_transcript is not None:
            if mikado_lib.loci_objects.superlocus.Superlocus.in_locus(current_locus, current_transcript) is True:
                current_locus.add_transcript_to_locus(current_transcript)
            else:
                if self.json_conf["single_thread"] is True:
                    analyse_locus(current_locus,
                                  self.json_conf,
                                  self.printer_queue,
                                  self.logging_queue,
                                  None,
                                  data_dict)
                else:
                    # while len(jobs) >= self.threads:
                    #     for job in jobs:
                    #         if job.is_alive() is False:
                    #             jobs.remove(job)
                    #
                    # job = Process(target=analyse_locus, args=(current_locus,
                    #                                           self.json_conf,
                    #                                           self.printer_queue,
                    #                                           self.logging_queue,
                    #                                           self.connection_pool
                    #                                           ))
                    # job.start()
                    # jobs.append(job)
                    jobs.append(pool.apply_async(analyse_locus, args=(current_locus,
                                                                      self.json_conf,
                                                                      self.printer_queue,
                                                                      self.logging_queue,
                                                                      None,
                                                                      data_dict
                                                                      )))

                current_locus = mikado_lib.loci_objects.superlocus.Superlocus(current_transcript, stranded=False,
                                                                              json_dict=self.json_conf)

        if current_locus is not None:
            # while len(multiprocessing.active_children()) >= self.threads + 2:
            #         continue
            if self.json_conf["single_thread"] is True:
                analyse_locus(current_locus,
                              self.json_conf,
                              self.printer_queue,
                              self.logging_queue,
                              None,
                              data_dict
                              # self.connection_pool
                              )
            else:
                # while len(jobs) >= self.threads:
                #     for job in jobs:
                #         if job.is_alive() is False:
                #             jobs.remove(job)
                #
                # job = Process(target=self.analyse_locus, args=(current_locus,
                #                                                self.json_conf,
                #                                                self.printer_queue,
                #                                                self.logging_queue,
                #                                                self.connection_pool
                #                                                ))
                # job.start()
                # jobs.append(job)
                jobs.append(pool.apply_async(analyse_locus, args=(current_locus,
                                                                  self.json_conf,
                                                                  self.printer_queue,
                                                                  self.logging_queue,
                                                                  None,
                                                                  data_dict
                                                                  )))
        for job in jobs:
            job.get()

        # while len(jobs) > 0:
        #     for job in jobs:
        #         if job.is_alive() is False:
        #             jobs.remove(job)
        #
        pool.close()
        pool.join()

        self.printer_queue.join()
        self.printer_queue.put("EXIT")
        # The printing process must be started AFTER we have put the stopping signal  into the queue
        self.printer_process.join()
        self.main_logger.info("Finished analysis of {0}".format(self.input_file))
        self.log_writer.stop()
        return 0
