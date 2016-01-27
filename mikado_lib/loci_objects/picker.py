#!/usr/bin/env python3
# coding=utf-8

"""
This module defines the Picker class, which is the main workhorse for mikado pick.
"""

import sys
import re
import csv
import os
import shutil
import tempfile
import logging
from logging import handlers as logging_handlers
import collections
import functools
import multiprocessing
from sqlalchemy.engine import create_engine  # SQLAlchemy/DB imports
from sqlalchemy.orm.session import sessionmaker
from sqlalchemy.pool import QueuePool as SqlPool
import sqlalchemy
from ..utilities import path_join
from ..parsers.GTF import GTF
from ..parsers.GFF import GFF3
from ..serializers.blast_serializer import Hit, Query
from ..serializers.junction import Junction, Chrom
from ..serializers.orf import Orf
from .superlocus import Superlocus, Transcript
# I have to use an absolute import to make nosetests import function ... strange
from mikado_lib.configuration.configurator import to_json
# from .. import configuration
from ..utilities import dbutils
from ..exceptions import UnsortedInput, InvalidJson
# import mikado_lib.exceptions
# from concurrent.futures import ProcessPoolExecutor

# pylint: disable=no-name-in-module
from multiprocessing import Process  # , Pool
# pylint: enable=no-name-in-module
import multiprocessing.managers


def remove_fragments(stranded_loci, json_conf, logger):

    """This method checks which loci are possible fragments, according to the
    parameters provided in the configuration file, and tags/remove them according
    to the configuration specification.

    :param stranded_loci: a list of the loci to consider for fragment removal
    :type stranded_loci: list[Superlocus]

    :param json_conf: the configuration dictionary
    :type json_conf: dict

    :param logger: the logger
    :type logger: logging.Logger

    """

    loci_to_check = {True: set(), False: set()}
    for stranded_locus in stranded_loci:
        for _, locus_instance in stranded_locus.loci.items():
            locus_instance.logger = logger
            loci_to_check[locus_instance.monoexonic].add(locus_instance)

    mcdl = json_conf["pick"]["run_options"]["fragments_maximal_cds"]
    bool_remove_fragments = json_conf["pick"]["run_options"]["remove_overlapping_fragments"]
    for stranded_locus in stranded_loci:
        to_remove = set()
        for locus_id, locus_instance in stranded_locus.loci.items():
            if locus_instance in loci_to_check[True]:
                logger.debug("Checking if %s is a fragment", locus_instance.id)

                for other_locus in loci_to_check[False]:
                    if other_locus.other_is_fragment(locus_instance,
                                                     minimal_cds_length=mcdl) is True:
                        if bool_remove_fragments is False:
                            # Just mark it as a fragment
                            stranded_locus.loci[locus_id].is_fragment = True
                        else:
                            to_remove.add(locus_id)
                            # del stranded_locus.loci[locus_id]
                        break
        for locus_id in to_remove:
            del stranded_locus.loci[locus_id]
        yield stranded_locus


def analyse_locus(slocus: Superlocus,
                  counter: int,
                  json_conf: dict,
                  printer_queue: [multiprocessing.managers.AutoProxy, None],
                  logging_queue: multiprocessing.managers.AutoProxy) -> [Superlocus]:

    """
    :param slocus: a superlocus instance
    :type slocus: mikado_lib.loci_objects.superlocus.Superlocus

    :param counter: an integer which is used to create the proper name for the locus.
    :type counter: int

    :param json_conf: the configuration dictionary
    :type json_conf: dict

    :param logging_queue: the logging queue
    :type logging_queue: multiprocessing.managers.AutoProxy

    :param printer_queue: the printing queue
    :type printer_queue: multiprocessing.managers.AutoProxy

    # :param data_dict: a dictionary of preloaded data
    # :type data_dict: dict

    This function takes as input a "superlocus" instance and the pipeline configuration.
    It also accepts as optional keywords a dictionary with the CDS information
    (derived from a Bed12Parser) and a "lock" used for avoiding writing collisions
    during multithreading.
    The function splits the superlocus into its strand components and calls the relevant methods
    to define the loci.
    When it is finished, it transmits the superloci to the printer function.
    """

    # Define the logger
    if slocus is None:
        # printer_dict[counter] = []
        if printer_queue:
            printer_queue.put_nowait(([], counter))
            return
        else:
            return []

    handler = logging_handlers.QueueHandler(logging_queue)
    logger = logging.getLogger("{0}:{1}-{2}".format(
        slocus.chrom, slocus.start, slocus.end))
    logger.addHandler(handler)

    # We need to set this to the lowest possible level,
    # otherwise we overwrite the global configuration
    logger.setLevel(json_conf["log_settings"]["log_level"])
    logger.propagate = False
    logger.debug("Started with %s, counter %d",
                 slocus.id, counter)
    if slocus.stranded is True:
        logger.warn("%s is stranded already! Resetting",
                    slocus.id)
        slocus.stranded = False

    slocus.logger = logger

    # Load the CDS information if necessary
    if slocus.initialized is False:
        # This happens when all transcripts have been removed from the locus,
        # due to errors that have been hopefully logged
        logger.warning(
            "%s had all transcripts failing checks, ignoring it",
            slocus.id)
        # printer_dict[counter] = []
        if printer_queue:
            printer_queue.put_nowait(([], counter))
            return
        else:
            return []

    # Split the superlocus in the stranded components
    logger.debug("Splitting by strand")
    stranded_loci = sorted(list(slocus.split_strands()))
    # Define the loci
    logger.debug("Divided into %d loci", len(stranded_loci))

    for stranded_locus in stranded_loci:
        stranded_locus.define_loci()
        logger.debug("Defined loci for %s:%f-%f, strand: %s",
                     stranded_locus.chrom,
                     stranded_locus.start,
                     stranded_locus.end,
                     stranded_locus.strand)

    # Remove overlapping fragments.
    loci_to_check = {True: set(), False: set()}
    for stranded_locus in stranded_loci:
        for _, locus_instance in stranded_locus.loci.items():
            locus_instance.logger = logger
            loci_to_check[locus_instance.monoexonic].add(locus_instance)

    # Check if any locus is a fragment, if so, tag/remove it
    stranded_loci = sorted(list(remove_fragments(stranded_loci, json_conf, logger)))
    try:
        logger.debug("Size of the loci to send: {0}, for {1} loci".format(
            sys.getsizeof(stranded_loci),
            len(stranded_loci)))
    except Exception as err:
        logger.error(err)
        pass
    # printer_dict[counter] = stranded_loci
    if printer_queue:
        printer_queue.put_nowait((stranded_loci, counter))
        logger.debug("Finished with %s, counter %d", slocus.id, counter)
        logger.removeHandler(handler)
        handler.close()
        return
    else:
        logger.debug("Finished with %s, counter %d", slocus.id, counter)
        logger.removeHandler(handler)
        handler.close()
        return stranded_loci


class LociProcesser(Process):

    def __init__(self, json_conf, data_dict, locus_queue, printer_queue, logging_queue):

        self.logging_queue = logging_queue
        self.json_conf = json_conf
        self.handler = logging_handlers.QueueHandler(self.logging_queue)
        self.logger = logging.getLogger("LociProcesser")
        self.logger.addHandler(self.handler)
        self.logger.setLevel(self.json_conf["log_settings"]["log_level"])
        self.logger.debug("Starting Process")

        self.data_dict = data_dict
        self.locus_queue = locus_queue
        self.printer_queue = printer_queue
        super(LociProcesser, self).__init__()

        self.logger.debug("Starting the pool for {0}".format(self.name))
        try:
            if self.json_conf["pick"]["run_options"]["preload"] is False:
                db_connection = functools.partial(dbutils.create_connector,
                                                  self.json_conf,
                                                  self.logger)
                self.connection_pool = sqlalchemy.pool.QueuePool(db_connection,
                                                                 pool_size=1,
                                                                 max_overflow=2)

            else:
                self.connection_pool = None
        except KeyboardInterrupt:
            raise
        except EOFError:
            raise
        except Exception as exc:
            self.logger.exception(exc)
            return

    def run(self):

        self.logger.debug("Starting to parse data for {0}".format(self.name))
        while True:
            slocus, counter = self.locus_queue.get()
            if slocus == "EXIT":
                self.locus_queue.put((slocus, counter))
                if self.connection_pool is not None:
                    self.connection_pool.dispose()
                return

            if slocus is not None:
                try:
                    slocus.load_all_transcript_data(pool=self.connection_pool,
                                                    data_dict=self.data_dict)
                except KeyboardInterrupt:
                    raise
                except Exception as exc:
                    self.logger.error("Error while loading data for %s",
                                      slocus.id)
                    self.logger.exception(exc)
                self.logger.debug("Loading transcript data for %s", slocus.id)
            analyse_locus(slocus, counter,
                          self.json_conf,
                          self.printer_queue,
                          self.logging_queue)


# pylint: disable=too-many-instance-attributes
class Picker:

    """
    This class is used to launch the main mikado_lib pipeline. Its purpose is to parse
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

        # Things that have to be deleted upon serialisation
        self.not_pickable = ["queue_logger", "manager", "printer_process",
                             "log_process", "pool", "main_logger",
                             "log_handler", "log_writer", "logger"]

        # Now we start the real work
        if isinstance(json_conf, str):
            assert os.path.exists(json_conf)
            json_conf = to_json(json_conf)
        else:
            assert isinstance(json_conf, dict)

        self.commandline = commandline
        self.json_conf = json_conf

        self.threads = self.json_conf["pick"]["run_options"]["threads"]
        self.input_file = self.json_conf["pick"]["files"]["input"]
        _ = self.define_input()  # Check the input file

        if self.json_conf["pick"]["files"]["subloci_out"]:
            self.sub_out = path_join(
                self.json_conf["pick"]["files"]["output_dir"],
                self.json_conf["pick"]["files"]["subloci_out"]
            )
        else:
            self.sub_out = ""
        if self.json_conf["pick"]["files"]["monoloci_out"]:
            self.monolocus_out = path_join(
                self.json_conf["pick"]["files"]["output_dir"],
                self.json_conf["pick"]["files"]["monoloci_out"]
            )
        else:
            self.monolocus_out = ""
        self.locus_out = path_join(
            self.json_conf["pick"]["files"]["output_dir"],
            self.json_conf["pick"]["files"]["loci_out"])
        # pylint: disable=no-member
        self.context = multiprocessing.get_context()
        # pylint: enable=no-member
        self.manager = self.context.Manager()
        # self.result_dict = self.manager.dict()
        self.printer_queue = self.manager.Queue(-1)
        self.logging_queue = self.manager.Queue(-1)

        # self.printer_queue = self.manager.Queue(-1)
        # self.logging_queue = self.manager.Queue(-1)  # queue for logging

        self.setup_logger()
        self.db_connection = functools.partial(
            dbutils.create_connector,
            self.json_conf,
            self.logger)
        self.logger_queue_handler = logging_handlers.QueueHandler(self.logging_queue)
        self.queue_logger = logging.getLogger("parser")
        self.queue_logger.addHandler(self.logger_queue_handler)

        # Configure SQL logging
        sqllogger = logging.getLogger("sqlalchemy.engine")
        # if json_conf["log_settings"]["log_level"] == "DEBUG":
        #     sqllogger.setLevel("DEBUG")
        # else:
        sqllogger.setLevel(json_conf["log_settings"]["sql_level"])
        sqllogger.addHandler(self.logger_queue_handler)

        # We need to set this to the lowest possible level,
        # otherwise we overwrite the global configuration
        self.queue_logger.setLevel(self.json_conf["log_settings"]["log_level"])
        self.queue_logger.propagate = False
        if self.json_conf["pick"]["run_options"]["single_thread"] is True:
            # Reset threads to 1
            if self.json_conf["pick"]["run_options"]["threads"] > 1:
                self.main_logger.warning("Reset number of threads to 1 as requested")
                self.threads = 1
        elif self.json_conf["pick"]["run_options"]["threads"] == 1:
            self.json_conf["pick"]["run_options"]["single_thread"] = True

        if self.json_conf["pick"]["run_options"]["preload"] is True and self.threads > 1:
            self.logger.warning("Preloading using multiple threads can be extremely \
                                memory intensive, proceed with caution!")
        #     self.json_conf["pick"]["run_options"]["threads"] = 1
        #     self.json_conf["pick"]["run_options"]["single_thread"] = True

        if self.locus_out is None:
            raise InvalidJson(
                "No output prefix specified for the final loci. Key: \"loci_out\"")

        self.printer_process = Process(target=self.printer)
        self.queue_pool = None

    def define_input(self):
        """Function to check that the input file exists and is valid. It returns the parser."""

        if self.input_file.endswith(".gtf"):
            parser = GTF
        else:
            parser = GFF3

        verified = False
        for row in parser(self.input_file):
            if row.header is False:
                verified = True
                break
        if verified is False:
            raise InvalidJson(
                "Invalid input file: {0}".format(self.input_file))

        return parser(self.input_file)

    def setup_shm_db(self):
        """
        This method will copy the SQLite input DB into memory.
        """

        self.main_logger.info("Copy into a SHM db: %s",
                              self.json_conf["pick"]["run_options"]["shm"])
        if self.json_conf["pick"]["run_options"]["shm"] is True:
            self.json_conf["pick"]["run_options"]["shm_shared"] = False
            self.main_logger.info("Copying the DB into memory")
            assert self.json_conf["db_settings"]["dbtype"] == "sqlite"
            self.json_conf["pick"]["run_options"]["preload"] = False
            if self.json_conf["pick"]["run_options"]["shm_db"] is not None:
                self.json_conf["pick"]["run_options"]["shm_db"] = os.path.join(
                    "/dev/shm/",
                    self.json_conf["pick"]["run_options"]["shm_db"])
                self.json_conf["pick"]["run_options"]["shm_shared"] = True
            else:
                # Create temporary file
                temp = tempfile.mktemp(suffix=".db",
                                       prefix="/dev/shm/")
                if os.path.exists(temp):
                    os.remove(temp)
                self.json_conf["pick"]["run_options"]["shm_db"] = temp
            if self.json_conf["pick"]["run_options"]["shm"]:
                if not os.path.exists(self.json_conf["pick"]["run_options"]["shm_db"]):
                    self.main_logger.info("Copying {0} into {1}".format(
                        self.json_conf["db_settings"]["db"],
                        self.json_conf["pick"]["run_options"]["shm_db"]))
                    try:
                        shutil.copy2(self.json_conf["db_settings"]["db"],
                                     self.json_conf["pick"]["run_options"]["shm_db"])
                    except PermissionError:
                        self.main_logger.warn(
                            """Permission to write on /dev/shm denied.
                            Back to using the DB on disk.""")
                        self.json_conf["pick"]["run_options"]["shm"] = False
                else:
                    self.main_logger.info("%s exists already. Doing nothing.",
                                          self.json_conf["pick"]["run_options"]["shm_db"])
            self.main_logger.info("DB copied into memory")

    def setup_logger(self):

        """This function sets up the logger for the class.
        It creates the instance attribute "log_writer", which is itself a
        logging.handlers.QueueListener instance listening on the logging_queue
        instance attribute (which is a normal mp.Manager.Queue instance)."""

        self.formatter = logging.Formatter(
            "{asctime} - {levelname} - {module}:{lineno} - {funcName} - {name} - {message}",
            style="{"
            )
        self.main_logger = logging.getLogger("main_logger")
        if not os.path.exists(self.json_conf["pick"]["files"]["output_dir"]):
            try:
                os.makedirs(self.json_conf["pick"]["files"]["output_dir"])
            except (OSError,PermissionError) as exc:
                self.logger.error("Failed to create the output directory!")
                self.logger.exception(exc)
                raise
        elif not os.path.isdir(self.json_conf["pick"]["files"]["output_dir"]):
            self.logger.error(
                "The specified output directory %s exists and is not a file; aborting",
                self.json_conf["pick"]["files"]["output_dir"])
            raise OSError("The specified output directory %s exists and is not a file; aborting" %
                          self.json_conf["pick"]["files"]["output_dir"])

        self.logger = logging.getLogger("listener")
        self.logger.propagate = False
        if (self.json_conf["pick"]["files"]["log"] is None or
                self.json_conf["pick"]["files"]["log"] == "stream"):
            self.log_handler = logging.StreamHandler()
        else:
            self.log_handler = logging.FileHandler(
                path_join(
                    self.json_conf["pick"]["files"]["output_dir"],
                    self.json_conf["pick"]["files"]["log"]), 'w')
        # For the main logger I want to keep it at the "INFO" level
        self.log_level = self.json_conf["log_settings"]["log_level"]

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
            self.main_logger.info(
                "Analysis launched directly, without using the launch script.")

        # Create the shared DB if necessary
        self.setup_shm_db()

        if self.json_conf["pick"]["chimera_split"]["blast_check"] is True and \
                self.json_conf["log_settings"]["log_level"] == "DEBUG":
            engine = dbutils.connect(self.json_conf, self.main_logger)
            smaker = sessionmaker()
            smaker.configure(bind=engine)
            session = smaker()

            evalue = self.json_conf["pick"]["chimera_split"]["blast_params"]["evalue"]
            queries_with_hits = session.query(
                Hit.query_id).filter(Hit.evalue <= evalue).distinct().count()
            total_queries = session.query(Query).count()
            self.main_logger.debug(
                "Queries with at least one hit at evalue<=%f: %d out of %d (%f%%)",
                evalue,
                queries_with_hits,
                total_queries,
                0 if total_queries == 0 else round(100 * queries_with_hits / total_queries, 2)
            )
            session.close()

        self.log_writer = logging_handlers.QueueListener(
            self.logging_queue, self.logger)
        self.log_writer.start()

        return

    def __print_gff_headers(self, locus_out, score_keys):
        """Private method to print the GFF headers of the output files.
        Moreover, it will determine whether to start output files for
        subloci and/or monoloci.
        Returns: [subloci_metrics_handle, subloci_scores_handle, subloci_handle],
         monosubloci_outfile
        """

        print('##gff-version 3', file=locus_out)
        engine = create_engine("{0}://".format(self.json_conf["db_settings"]["dbtype"]),
                               creator=self.db_connection)
        session = sqlalchemy.orm.sessionmaker(bind=engine)()
        for chrom in session.query(Chrom).order_by(Chrom.name.desc()):
            print("##sequence-region {0} 1 {1}".format(chrom.name, chrom.length),
                  file=locus_out)

        if self.sub_out != '':
            assert isinstance(self.sub_out, str)
            sub_metrics_file = re.sub("$", ".metrics.tsv",
                                      re.sub(".gff.?$", "", self.sub_out))
            sub_scores_file = re.sub("$", ".scores.tsv",
                                     re.sub(".gff.?$", "", self.sub_out))
            sub_metrics = csv.DictWriter(
                open(sub_metrics_file, 'w'),
                Superlocus.available_metrics,
                delimiter="\t")
            sub_metrics.writeheader()
            sub_scores = csv.DictWriter(
                open(sub_scores_file, 'w'), score_keys, delimiter="\t")
            sub_scores.writeheader()
            sub_out = open(self.sub_out, 'w')
            print('##gff-version 3', file=sub_out)
            for chrom in session.query(Chrom).order_by(Chrom.name.desc()):
                print("##sequence-region {0} 1 {1}".format(chrom.name, chrom.length),
                      file=sub_out)
            sub_files = [sub_metrics, sub_scores, sub_out]
        else:
            sub_files = [None, None, None]

        if self.monolocus_out != '':
            mono_out = open(self.monolocus_out, 'w')
            print('##gff-version 3', file=mono_out)
            for chrom in session.query(Chrom).order_by(Chrom.name.desc()):
                print("##sequence-region {0} 1 {1}".format(chrom.name, chrom.length),
                      file=mono_out)
        else:
            mono_out = None

        session.close()
        engine.dispose()
        return sub_files, mono_out

    def __get_output_files(self):

        """
        Private method used by printer to prepare all the output files.
        :return: ((locus_metrics, locus_scores, locus_out),
                (sub_metrics, sub_scores, sub_out),
                mono_out)
                all of these are file handles
        """

        score_keys = ["tid", "parent", "score"]
        score_keys += sorted(list(self.json_conf["scoring"].keys()))
        # Define mandatory output files
        locus_metrics_file = re.sub("$", ".metrics.tsv", re.sub(
            ".gff.?$", "", self.locus_out))
        locus_scores_file = re.sub("$", ".scores.tsv", re.sub(
            ".gff.?$", "", self.locus_out))
        locus_metrics = csv.DictWriter(
            open(locus_metrics_file, 'w'),
            Superlocus.available_metrics,
            delimiter="\t")

        locus_metrics.writeheader()
        locus_scores = csv.DictWriter(open(locus_scores_file, 'w'), score_keys, delimiter="\t")
        locus_scores.writeheader()
        locus_out = open(self.locus_out, 'w')
        sub_files, mono_out = self.__print_gff_headers(locus_out, score_keys)

        return ((locus_metrics, locus_scores, locus_out),
                sub_files, mono_out)

    # This method has many local variables, but most (9!) are
    # actually file handlers. I cannot trim them down for now.
    # pylint: disable=too-many-locals

    def _print_locus(self, stranded_locus, gene_counter, logger=None, handles=()):

        """
        Private method that handles a single superlocus for printing.
        It also detects and flags/discard fragmentary loci.
        :param stranded_locus: the stranded locus to analyse
        :param gene_counter: A counter used to rename the genes/transcripts progressively
        :param logger: logger instance
        :param handles: the handles to print to
        :return:
        """

        locus_metrics, locus_scores, locus_out = handles[0]
        sub_metrics, sub_scores, sub_out = handles[1]
        mono_out = handles[2]

        stranded_locus.logger = logger
        if self.sub_out != '':  # Skip this section if no sub_out is defined
            sub_lines = stranded_locus.__str__(
                level="subloci",
                print_cds=not self.json_conf["pick"]["run_options"]["exclude_cds"])
            if sub_lines != '':
                print(sub_lines, file=sub_out)
            sub_metrics_rows = [x for x in stranded_locus.print_subloci_metrics()
                                if x != {} and "tid" in x]
            sub_scores_rows = [x for x in stranded_locus.print_subloci_scores()
                               if x != {} and "tid" in x]
            for row in sub_metrics_rows:
                sub_metrics.writerow(row)
            for row in sub_scores_rows:
                sub_scores.writerow(row)
        if self.monolocus_out != '':
            mono_lines = stranded_locus.__str__(
                level="monosubloci",
                print_cds=not self.json_conf["pick"]["run_options"]["exclude_cds"])
            if mono_lines != '':
                print(mono_lines, file=mono_out)
        locus_metrics_rows = [x for x in stranded_locus.print_monoholder_metrics()
                              if x != {} and "tid" in x]
        locus_scores_rows = [x for x in stranded_locus.print_monoholder_scores()
                             if x != {} and "tid" in x]

        for locus in stranded_locus.loci:
            gene_counter += 1
            fragment_test = (
                self.json_conf["pick"]["run_options"]["remove_overlapping_fragments"]
                is True and stranded_locus.loci[locus].is_fragment is True)

            if fragment_test is True:
                continue
            new_id = "{0}.{1}G{2}".format(
                self.json_conf["pick"]["output_format"]["id_prefix"],
                stranded_locus.chrom, gene_counter)
            stranded_locus.loci[locus].id = new_id

        locus_lines = stranded_locus.__str__(
            print_cds=not self.json_conf["pick"]["run_options"]["exclude_cds"])
        for row in locus_metrics_rows:
            locus_metrics.writerow(row)
        for row in locus_scores_rows:
            locus_scores.writerow(row)

        if locus_lines != '':
            print(locus_lines, file=locus_out)
        return gene_counter

    def printer(self):

        """Listener process that will print out the loci recovered by
        the analyse_locus function."""

        handler = logging_handlers.QueueHandler(self.logging_queue)
        logger = logging.getLogger("queue_listener")
        logger.propagate = False
        logger.addHandler(handler)
        logger.setLevel(self.json_conf["log_settings"]["log_level"])

        handles = self.__get_output_files()

        locus_printer = functools.partial(self._print_locus,
                                          logger=logger,
                                          handles=handles)

        last_printed = -1
        cache = dict()
        curr_chrom = None
        gene_counter = 0
        logger.debug("Starting to wait for loci to print")

        while True:
            current = self.printer_queue.get()
            stranded_loci, counter = current

            cache[counter] = stranded_loci
            if stranded_loci != "EXIT":
                logger.debug("Received one locus, counter: %d, total %d. names %s",
                               counter,
                               len(stranded_loci),
                               ", ".join([_.id for _ in stranded_loci]))
            if counter == last_printed + 1 or stranded_loci == "EXIT":
                cache[counter] = stranded_loci
                for num in sorted(cache.keys()):
                    if num > last_printed + 1:
                        break
                    else:
                        if num % 1000 == 0 and num > 0:
                            logger.info("Printed %d superloci", num)
                        for stranded_locus in cache[num]:
                            if stranded_locus.chrom != curr_chrom:
                                curr_chrom = stranded_locus.chrom
                                gene_counter = 0
                            gene_counter = locus_printer(stranded_locus, gene_counter)
                        last_printed += 1
                        del cache[num]
                if stranded_loci == "EXIT":
                    self.printer_queue.task_done()
                    logger.info("Final number of superloci: %d", last_printed)
                    break
            self.printer_queue.task_done()

        return
    # pylint: enable=too-many-locals

    def __getstate__(self):

        state = self.__dict__.copy()
        for not_pickable in self.not_pickable:
            if not_pickable in state:
                del state[not_pickable]

        return state

    # pylint: disable=too-many-locals
    def __preload_blast(self, engine, queries):

        """Private method to load all the blast information
        into a specific dictionary.

        :param engine: the connection engine from the preloading thread
        :param queries: a dictionary containing the name=>ID relationship
         for queries
        :returns hits_dict: a dictionary with the loaded BLAST data
        """

        hits_dict = collections.defaultdict(list)
        hsps = dict()
        for hsp in engine.execute("select * from hsp"):
            if hsp.query_id not in hsps:
                hsps[hsp.query_id] = collections.defaultdict(list)
            hsps[hsp.query_id][hsp.target_id].append(hsp)

        self.main_logger.info("{0} HSPs prepared".format(len(hsps)))

        targets = dict((x.target_id, x) for x in engine.execute("select * from target"))

        hit_counter = 0
        hits = engine.execute(
            " ".join(["select * from hit where evalue <= {0} and hit_number <= {1}",
                      "order by query_id, evalue asc;"]).format(
                    self.json_conf["pick"]["chimera_split"]["blast_params"]["evalue"],
                    self.json_conf["pick"]["chimera_split"]["blast_params"]["max_target_seqs"]))

        # self.main_logger.info("{0} BLAST hits to analyse".format(hits))
        current_counter = 0
        current_hit = None
        previous_evalue = -1

        max_targets = self.json_conf["pick"]["chimera_split"]["blast_params"]["max_target_seqs"]
        for hit in hits:
            if current_hit != hit.query_id:
                current_hit = hit.query_id
                current_counter = 0
                previous_evalue = -1

            if current_counter > max_targets and previous_evalue < hit.evalue:
                continue
            elif previous_evalue < hit.evalue:
                previous_evalue = hit.evalue

            current_counter += 1

            my_query = queries[hit.query_id]
            my_target = targets[hit.target_id]

            # We HAVE to use the += approach because extend/append
            # leaves the original list empty
            hits_dict[my_query.query_name].append(
                Hit.as_full_dict_static(
                    hit,
                    hsps[hit.query_id][hit.target_id],
                    my_query,
                    my_target
                )
            )
            hit_counter += 1
            if hit_counter >= 2*10**4 and hit_counter % (2*10**4) == 0:
                self.main_logger.debug("Loaded %d BLAST hits in database",
                                       hit_counter)

        del hsps
        assert len(hits_dict) <= len(queries)
        self.main_logger.info("%d BLAST hits loaded for %d queries",
                              hit_counter,
                              len(hits_dict))
        self.main_logger.debug("%s",
                               ", ".join(
                                   [str(x) for x in list(hits_dict.keys())[:10]]))
        return hits_dict
    # pylint: enable=too-many-locals

    def preload(self):
        """
        This method preloads the data from the DB into a dictionary ("data_dict").
        The information on what to extract and how to connect to the
        DB is retrieved from the json_conf dictionary.
        :return: data_dict
        :rtype: dict
        """

        self.main_logger.info("Starting to preload the database into memory")

        # data_dict = self.manager.dict(lock=False)
        data_dict = dict()
        engine = create_engine("{0}://".format(self.json_conf["db_settings"]["dbtype"]),
                               creator=self.db_connection)
        session = sqlalchemy.orm.sessionmaker(bind=engine)()

        junc_dict = dict()
        for junc in session.query(Junction):
            key = (junc.chrom, junc.junction_start, junc.junction_end)
            assert key not in junc_dict
            junc_dict[key] = junc.strand
        data_dict["junctions"] = junc_dict

        # data_dict["junctions"] = self.manager.dict(data_dict["junctions"], lock=False)

        self.main_logger.info("%d junctions loaded",
                              len(data_dict["junctions"]))
        self.main_logger.debug("Example junctions:\n{0}".format(
            "\n".join(str(junc) for junc in list(
                data_dict["junctions"])[:min(10, len(data_dict["junctions"]))])))

        queries = dict((que.query_id, que) for que in engine.execute("select * from query"))

        # Then load ORFs
        orf_dict = collections.defaultdict(list)

        for orf in engine.execute("select * from orf"):

            query_name = queries[orf.query_id].query_name
            orf_dict[query_name].append(
                Orf.as_bed12_static(orf, query_name)
            )

        data_dict["orfs"] = orf_dict
        assert len(data_dict["orfs"]) == engine.execute(
            "select count(distinct(query_id)) from orf").fetchone()[0]

        # data_dict['orf'] = self.manager.dict(orfs, lock=False)

        self.main_logger.info("%d ORFs loaded",
                              len(data_dict["orfs"]))
        self.main_logger.debug(",".join(
            list(data_dict["orfs"].keys())[:10]
        ))

        # Finally load BLAST

        # if self.json_conf["pick"]["chimera_split"]["execute"] is True and \
        #         self.json_conf["pick"]["chimera_split"]["blast_check"] is True:
        data_dict["hits"] = self.__preload_blast(engine, queries)
        # else:
        #     data_dict["hits"] = dict()
        #     self.main_logger.info("Skipping BLAST loading")

        self.main_logger.info("Finished to preload the database into memory")
        return data_dict

    def _submit_locus(self, slocus, counter, data_dict=None, pool=None):
        """
        Private method to submit / start the analysis of a superlocus in input.
        :param slocus: the locus to analyse.
        :param data_dict: the preloaded data in memory
        :param pool: thread execution pool
        :return: job object / None
        """

        if slocus is not None:
            self.main_logger.debug("Submitting %s", slocus.id)
            
        if self.json_conf["pick"]["run_options"]["preload"] is True and slocus is not None:
            self.logger.debug("Loading data from dict for %s", slocus.id)
            slocus.logger = self.logger
            slocus.load_all_transcript_data(pool=None,
                                            data_dict=data_dict)
            # slocus_id = slocus.id
            if slocus.initialized is False:
                # This happens when we have removed all transcripts from the locus
                # due to errors which should have been caught and logged
                self.logger.warning(
                    "%s had all transcripts failing checks, ignoring it",
                    slocus.id)
                # Exit
                return None

        if self.json_conf["pick"]["run_options"]["single_thread"] is True:
            return analyse_locus(slocus,
                                 counter,
                                 self.json_conf,
                                 None,
                                 self.logging_queue)
        else:
            job = pool.apply_async(analyse_locus,
                                   args=(
                                      slocus,
                                      counter,
                                      self.json_conf,
                                      self.printer_queue,
                                      self.logging_queue)
                                        )
            return job

    def __unsorted_interrupt(self, row, current_transcript, counter):
        """
        Private method that brings the program to a screeching halt
         if the GTF/GFF is not properly sorted.
        :param row:
        :param current_transcript:
        :return:
        """

        # self.result_dict[counter] = [_]
        self.printer_queue.put_nowait(("EXIT", float("inf")))
        current = "\t".join([str(x) for x in [row.chrom,
                                              row.start,
                                              row.end,
                                              row.strand]])
        previous = "\t".join([str(x) for x in [current_transcript.chrom,
                                               current_transcript.start,
                                               current_transcript.end,
                                               current_transcript.strand]])

        error_msg = """Unsorted input file, the results will not be correct.
            Please provide a properly sorted input. Error:
            {0}
            {1}""".format(current,
                          previous)
        self.logger.critical(error_msg)
        error_msg = "CRITICAL - {0}".format(" ".join(
            [l.strip() for l in error_msg.split("\n")]))
        raise UnsortedInput(error_msg)

    def __test_sortedness(self, row, current_transcript, counter):
        """
        Private method to test whether a row and the current transcript are actually in the expected
        sorted order.
        :param row:
        :param current_transcript:
        :return:
        """
        test = True
        if current_transcript.chrom > row.chrom:
            test = False
        elif (current_transcript.chrom == row.chrom and
              current_transcript.start > row.start):
            test = False
        elif (current_transcript.chrom == row.chrom and
              current_transcript.start == row.start and
              current_transcript.end > row.end):
            test = False

        if test is False:
            self.__unsorted_interrupt(row, current_transcript, counter)

    def __submit_multi_threading(self, data_dict):

        """
        Method to execute mikado pick in multi threaded mode.

        :param data_dict: The data dictionary
        :return:
        """

        intron_range = self.json_conf["soft_requirements"]["intron_range"]
        self.logger.info("Intron range: %s", intron_range)

        current_locus = None
        current_transcript = None

        self.printer_process.start()

        locus_queue = self.manager.Queue(-1)

        # with Pool(processes=self.threads, maxtasksperchild=None) as pool:
        working_processes = [LociProcesser(self.json_conf,
                                           data_dict,
                                           locus_queue,
                                           self.printer_queue, self.logging_queue)
                             for _ in range(self.threads)]
        # Start all processes
        [_.start() for _ in working_processes]
        # No sense in keeping this data available on the main thread now
        del data_dict

        counter = -1
        for row in self.define_input():
            if row.is_exon is True:
                current_transcript.add_exon(row)
            elif row.is_transcript is True:
                if current_transcript is not None:
                    self.__test_sortedness(row, current_transcript, counter)
                    if Superlocus.in_locus(
                            current_locus, current_transcript) is True:
                        current_locus.add_transcript_to_locus(current_transcript,
                                                              check_in_locus=False)
                    else:
                        counter += 1
                        self.logger.debug("Submitting locus # %d", counter)
                        locus_queue.put((current_locus, counter))
                        # submit_locus(current_locus, counter)
                        # if job is not None:
                        #     jobs.append(job)
                        current_locus = Superlocus(
                            current_transcript,
                            stranded=False,
                            json_conf=self.json_conf)

                if current_transcript is None or row.chrom != current_transcript.chrom:
                    if current_transcript is not None:
                        self.logger.info("Finished chromosome %s",
                                         current_transcript.chrom)
                    self.logger.info("Starting chromosome %s", row.chrom)

                current_transcript = Transcript(
                    row,
                    source=self.json_conf["pick"]["output_format"]["source"],
                    intron_range=intron_range)

        if current_transcript is not None:
            if Superlocus.in_locus(
                    current_locus, current_transcript) is True:
                current_locus.add_transcript_to_locus(
                    current_transcript, check_in_locus=False)
            else:
                counter += 1
                self.logger.debug("Submitting locus # %d", counter)
                locus_queue.put((current_locus, counter))

                # jobs.append(submit_locus(current_locus, counter))

                current_locus = Superlocus(
                    current_transcript,
                    stranded=False,
                    json_conf=self.json_conf)
                self.logger.debug("Created last locus %s",
                                  current_locus)
        self.logger.info("Finished chromosome %s", current_locus.chrom)

        counter += 1
        locus_queue.put((current_locus, counter))
        # submit_locus(current_locus, counter)
        self.logger.debug("Submitting locus %s, counter %d",
                          current_locus.id, counter)
        locus_queue.put(("EXIT", float("inf")))
        [_.join() for _ in working_processes]
        # pool.close()
        # pool.join()

        self.printer_queue.join()
        self.printer_queue.put_nowait(("EXIT", float("inf")))

        self.printer_process.join()

    def __submit_single_threaded(self, data_dict):

        """
        Method to execute mikado pick in single threaded mode.
        :param data_dict:
        :return:
        """

        current_locus = None
        current_transcript = None
        submit_locus = functools.partial(self._submit_locus, **{"data_dict": data_dict,
                                                                "pool": None})

        handler = logging_handlers.QueueHandler(self.logging_queue)
        logger = logging.getLogger("queue_listener")
        logger.propagate = False
        logger.addHandler(handler)
        logger.setLevel(self.json_conf["log_settings"]["log_level"])
        logger.debug("Begun single-threaded run")

        intron_range = self.json_conf["soft_requirements"]["intron_range"]
        logger.info("Intron range: %s", intron_range)

        handles = self.__get_output_files()

        locus_printer = functools.partial(self._print_locus,
                                          logger=logger,
                                          handles=handles)

        # last_printed = -1
        curr_chrom = None
        gene_counter = 0

        if self.json_conf["pick"]["run_options"]["preload"] is False:
            db_connection = functools.partial(dbutils.create_connector,
                                              self.json_conf,
                                              self.logger)
            self.connection_pool = sqlalchemy.pool.QueuePool(db_connection,
                                                             pool_size=1,
                                                             max_overflow=2)
        else:
            self.connection_pool = None

        counter = -1
        for row in self.define_input():
            if row.is_exon is True:
                current_transcript.add_exon(row)
            elif row.is_transcript is True:
                if current_transcript is not None:
                    self.__test_sortedness(row, current_transcript, counter)
                    if Superlocus.in_locus(
                            current_locus, current_transcript) is True:
                        current_locus.add_transcript_to_locus(current_transcript,
                                                              check_in_locus=False)
                    else:
                        counter += 1
                        self.logger.debug("Analysing locus # %d", counter)
                        current_locus.load_all_transcript_data(pool=self.connection_pool,
                                                               data_dict=data_dict)
                        for stranded_locus in submit_locus(current_locus, counter):
                            if stranded_locus.chrom != curr_chrom:
                                curr_chrom = stranded_locus.chrom
                                gene_counter = 0
                            gene_counter = locus_printer(stranded_locus, gene_counter)
                        current_locus = Superlocus(
                            current_transcript,
                            stranded=False,
                            json_conf=self.json_conf)

                if current_transcript is None or row.chrom != current_transcript.chrom:
                    if current_transcript is not None:
                        self.logger.info("Finished chromosome %s",
                                         current_transcript.chrom)
                    self.logger.info("Starting chromosome %s", row.chrom)

                current_transcript = Transcript(
                    row,
                    source=self.json_conf["pick"]["output_format"]["source"],
                    intron_range=intron_range)

        if current_transcript is not None:
            if Superlocus.in_locus(
                    current_locus, current_transcript) is True:
                current_locus.add_transcript_to_locus(
                    current_transcript, check_in_locus=False)
            else:
                counter += 1
                self.logger.debug("Analysing locus # %d", counter)
                current_locus.load_all_transcript_data(pool=self.connection_pool,
                                                       data_dict=data_dict)
                for stranded_locus in submit_locus(current_locus, counter):
                                if stranded_locus.chrom != curr_chrom:
                                    curr_chrom = stranded_locus.chrom
                                    gene_counter = 0
                                gene_counter = locus_printer(stranded_locus, gene_counter)

                current_locus = Superlocus(
                    current_transcript,
                    stranded=False,
                    json_conf=self.json_conf)
                self.logger.debug("Created last locus %s",
                                  current_locus)
        self.logger.info("Finished chromosome %s", current_locus.chrom)

        counter += 1
        self.logger.debug("Analysing locus # %d", counter)
        current_locus.load_all_transcript_data(pool=self.connection_pool,
                                               data_dict=data_dict)
        for stranded_locus in submit_locus(current_locus, counter):
            if stranded_locus.chrom != curr_chrom:
                curr_chrom = stranded_locus.chrom
                gene_counter = 0
            gene_counter = locus_printer(stranded_locus, gene_counter)
        # submit_locus(current_locus, counter)
        logger.info("Final number of superloci: %d", counter)

    def _parse_and_submit_input(self, data_dict):

        """
        This method does the parsing of the input and submission of the loci to the
        _submit_locus method.
        :param data_dict: The cached data from the database
        :return: jobs (the list of all jobs already submitted)
        """

        if self.json_conf["pick"]["run_options"]["single_thread"] is False:
            self.__submit_multi_threading(data_dict)
        else:
            self.__submit_single_threaded(data_dict)
        return

    def __call__(self):

        """This method will activate the class and start the analysis of the input file."""

        # NOTE: Pool, Process and Manager must NOT become instance attributes!
        # Otherwise it will raise all sorts of mistakes

        data_dict = None

        if self.json_conf["pick"]["run_options"]["preload"] is True:
            # Use the preload function to create the data dictionary
            data_dict = self.preload()
        # pylint: disable=no-member
        # pylint: enable=no-member

        self.logger.debug("Source: %s",
                          self.json_conf["pick"]["output_format"]["source"])
        if self.json_conf["db_settings"]["dbtype"] == "sqlite" and data_dict is not None:
            self.queue_pool = sqlalchemy.pool.QueuePool(
                self.db_connection,
                pool_size=self.threads,
                max_overflow=0)

        try:
            self._parse_and_submit_input(data_dict)
        except UnsortedInput as _:
            self.logger.error(
                "The input files were not properly sorted! Please run prepare and retry.")
            sys.exit(1)

        # list(map(job.get() for job in jobs if job is not None))
        # for job in iter(x for x in jobs if x is not None):
        #     job.get()

        self.log_writer.stop()
        if self.queue_pool is not None:
            self.queue_pool.dispose()

        # Clean up the DB copied to SHM
        if (self.json_conf["pick"]["run_options"]["shm"] is True and
                self.json_conf["pick"]["run_options"]["shm_shared"] is False):
            self.main_logger.info("Removing shared memory DB %s",
                                  self.json_conf["pick"]["run_options"]["shm_db"])
            os.remove(self.json_conf["pick"]["run_options"]["shm_db"])

        self.main_logger.info("Finished analysis of %s", self.input_file)
        return 0
# pylint: enable=too-many-instance-attributes
