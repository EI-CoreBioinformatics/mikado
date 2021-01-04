#!/usr/bin/env python3
# coding=utf-8

"""
This module defines the Picker class, which is the main workhorse for Mikado pick.
"""

import sys
import re
import csv
import os
import shutil
import tempfile
import random
import logging
from logging import handlers as logging_handlers
import functools
import multiprocessing
from sqlalchemy.engine import create_engine  # SQLAlchemy/DB imports
from sqlalchemy.orm.session import sessionmaker
from sqlalchemy.pool import QueuePool as SqlPool
import sqlalchemy
import sqlalchemy.exc
from ..utilities import path_join
from ..utilities.log_utils import formatter
from ..parsers.GTF import GTF, GtfLine
from ..parsers.GFF import GFF3
from ..parsers.bed12 import Bed12Parser
from ..parsers import Parser
from ..serializers.junction import Chrom
from ..serializers.external import ExternalSource
from ..transcripts import Transcript
from ..loci.superlocus import Superlocus
from ..configuration.configurator import to_json, check_json  # Necessary for nosetests
from ..utilities import dbutils
from ..exceptions import UnsortedInput, InvalidJson, InvalidTranscript
from .loci_processer import analyse_locus, LociProcesser, merge_loci
from ._locus_single_printer import print_locus
import multiprocessing.managers
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
import pickle
import warnings
import pyfaidx
import sqlite3
import msgpack
from numpy import percentile
logging.captureWarnings(True)
warnings.simplefilter("always")
try:
    import rapidjson as json
except ImportError:
    import json


# pylint: disable=too-many-instance-attributes
class Picker:

    """
    This class is used to launch the main Mikado pipeline. Its purpose is to parse
    an input sorted annotation file, locate the loci, and perform the selection analysis
    using the parameters provided in the input configuration file.
    """

    # @profile
    def __init__(self, json_conf, commandline="", regions=None):

        """Constructor. It takes a single argument as input - the JSON/YAML configuration,
        prepared by the json_utils functions.

        :param json_conf: Either a configuration dictionary or the configuration file.

        :param commandline: optional, the commandline used to start the program
        :type commandline: str
        """

        # Mock variables
        self.formatter = self.main_logger = self.log_writer = self.log_handler = self.logger = None
        self.log_level = "WARN"

        # # Things that have to be deleted upon serialisation
        # self.not_pickable = ["queue_logger", "manager", "printer_process",
        #                      "log_process", "pool", "main_logger",
        #                      "log_handler", "log_writer", "logger", "engine"]
        #
        # Now we start the real work
        self.commandline = commandline
        self.json_conf = json_conf

        self.__load_configuration()
        self.__regions = regions
        self.regressor = None

        self.procs = self.json_conf["threads"]

        # Check the input file
        with self.define_input() as _:
            pass

        self.__create_output_handles()
        # pylint: disable=no-member
        multiprocessing.set_start_method(self.json_conf["multiprocessing_method"],
                                         force=True)

        # self.setup_logger()
        self.logger.info("Random seed: %s", self.json_conf["seed"])
        if self.json_conf["seed"] is not None:
            # numpy.random.seed((self.json_conf["seed"]) % (2 ** 32 - 1))
            random.seed((self.json_conf["seed"]) % (2 ** 32 - 1))
        else:
            # numpy.random.seed(None)
            random.seed(None)
        self.logger.debug("Multiprocessing method: %s", self.json_conf["multiprocessing_method"])

        # pylint: enable=no-member
        self.manager = self.context.Manager()

        self.db_connection = functools.partial(
            dbutils.create_connector,
            self.json_conf,
            self.logger)

        # Update the database if necessary
        engine = create_engine("{0}://".format(self.json_conf["db_settings"]["dbtype"]), creator=self.db_connection)
        dbutils.DBBASE.metadata.create_all(engine)
        engine.dispose()

        if self.json_conf["pick"]["run_options"]["single_thread"] is True:
            # Reset threads to 1
            if self.json_conf["threads"] > 1:
                self.main_logger.warning("Reset number of threads to 1 as requested")
                self.procs = 1
        # elif self.json_conf["threads"] == 1:
        #     self.json_conf["pick"]["run_options"]["single_thread"] = True

        if self.locus_out is None:
            raise InvalidJson(
                "No output prefix specified for the final loci. Key: \"loci_out\"")

        # self.printer_process = Process(target=self.printer)
        self.queue_pool = None

    def define_input(self, multithreading=False):
        """Function to check that the input file exists and is valid. It returns the parser."""

        if ".gtf" in self.input_file:
            parser = GTF
            fformat = "gtf"
        elif ".gff" in self.input_file:
            parser = GFF3
            fformat = "gff3"
        elif ".bed" in self.input_file:
            parser = Bed12Parser
            fformat = "bed12"
        else:
            raise InvalidJson("Invalid input file: {0}".format(self.input_file))

        verified = False
        with parser(self.input_file) as testing:
            for row in testing:
                if row.header is False:
                    verified = True
                    break
        if verified is False:
            raise InvalidJson("Invalid input file: {0}".format(self.input_file))

        if multithreading:
            parser = Parser(self.input_file)
            parser.format = fformat
            return parser
        else:
            return parser(self.input_file)

    def __load_configuration(self):

        """Private method to load the configuration"""

        if isinstance(self.json_conf, str):
            assert os.path.exists(self.json_conf)
            self.json_conf = to_json(self.json_conf, logger=self.logger)
            # pylint: disable=no-member
            multiprocessing.set_start_method(self.json_conf["multiprocessing_method"],
                                             force=True)
            self.input_file = self.json_conf["pick"]["files"]["input"]
            self.setup_logger()
        elif isinstance(self.json_conf, dict):
            # pylint: disable=no-member
            self.input_file = self.json_conf["pick"]["files"]["input"]
            multiprocessing.set_start_method(self.json_conf["multiprocessing_method"],
                                             force=True)
            self.setup_logger()
            self.logger.debug("Checking the configuration dictionary")
            try:
                self.json_conf = check_json(self.json_conf, logger=self.logger)
                self.logger.debug("Configuration dictionary passes checks")
            except Exception as exc:
                self.logger.critical("Something went wrong with the configuration, critical error, aborting.")
                self.logger.critical(exc)
                sys.exit(1)
        else:
            raise TypeError(type(self.json_conf))
        assert isinstance(self.json_conf, dict)

        for key in ("remove_overlapping_fragments", "flank", "purge"):
            if key in self.json_conf["pick"]["run_options"]:
                # Put warnings in place for the deprecation of some options.

                if key == "remove_overlapping_fragments":
                    self.json_conf["pick"]["fragments"]["remove"] = self.json_conf["pick"]["run_options"].pop(key)
                    new_home = "fragments/remove"
                else:
                    self.json_conf["pick"]["clustering"][key] = self.json_conf["pick"]["run_options"].pop(key)
                    new_home = "clustering/{}".format(key)
                warns = PendingDeprecationWarning(
                    """The \"{}\" property has now been moved to pick/{}. \
Please update your configuration files in the future.""".format(
                        key, new_home))
                self.logger.warning(warns)

        if self.json_conf.get("pick", {}).get("alternative_splicing", {}).get("pad", False) is True:
            # Check that, when asks for padding, the reference genome is present
            self.logger.debug("Checking for the presence of the reference genome")
            try:
                _ = pyfaidx.Fasta(self.json_conf["reference"]["genome"])
            except (pyfaidx.FastaIndexingError, FileNotFoundError, pyfaidx.FastaNotFoundError):
                self.logger.error("Transcript padding cannot be executed without a valid genome file.\
                 Please, either disable the padding or provide a valid genome sequence.")
                sys.exit(1)
            self.logger.debug("Valid reference genome found")
        else:
            pass

        self.context = multiprocessing.get_context()
        if self.json_conf["pick"]["scoring_file"].endswith((".pickle", ".model")):
            with open(self.json_conf["pick"]["scoring_file"], "rb") as forest:
                self.regressor = pickle.load(forest)
            if not isinstance(self.regressor["scoring"], (RandomForestRegressor, RandomForestClassifier)):
                exc = TypeError("Invalid regressor provided, type: %s", type(self.regressor["scoring"]))
                self.logger.critical(exc)
                return
        else:
            self.regressor = None

        self.logger.debug("Configuration loaded successfully")

    def __create_output_handles(self):

        """Create all the output-related variables."""

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

        assert self.locus_out != ''
        assert self.locus_out != self.sub_out and self.locus_out != self.monolocus_out
        assert (not self.sub_out and not self.monolocus_out) or (self.sub_out != self.monolocus_out)

    def setup_shm_db(self):
        """
        This method will copy the SQLite input DB into memory.
        """

        if self.json_conf["pick"]["run_options"]["shm"] is True:
            self.main_logger.info("Copying Mikado database into a SHM db")
            assert self.json_conf["db_settings"]["dbtype"] == "sqlite"
            # Create temporary file
            temp = tempfile.mktemp(suffix=".db",
                                   prefix="/dev/shm/")
            if os.path.exists(temp):
                os.remove(temp)
            self.main_logger.debug("Copying {0} into {1}".format(
                self.json_conf["db_settings"]["db"],
                temp))
            try:
                shutil.copy2(self.json_conf["db_settings"]["db"],
                             temp)
                self.json_conf["db_settings"]["db"] = temp
            except PermissionError:
                self.main_logger.warning(
                    """Permission to write on /dev/shm denied.
                    Back to using the DB on disk.""")
                self.json_conf["pick"]["run_options"]["shm"] = False
            self.main_logger.info("DB copied into memory")

    def setup_logger(self):

        """This function sets up the logger for the class.
        It creates the instance attribute "log_writer", which is itself a
        logging.handlers.QueueListener instance listening on the logging_queue
        instance attribute (which is a normal mp.Manager.Queue instance)."""

        self.logging_queue = multiprocessing.Queue(-1)
        self.printer_queue = multiprocessing.Queue(-1)
        self.formatter = formatter
        self.main_logger = logging.getLogger("main_logger")
        if not os.path.exists(self.json_conf["pick"]["files"]["output_dir"]):
            try:
                os.makedirs(self.json_conf["pick"]["files"]["output_dir"])
            except (OSError, PermissionError) as exc:
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
                self.json_conf["pick"]["files"]["log"] in ("stream", "")):
            self.log_handler = logging.StreamHandler()
        else:
            if os.path.basename(self.json_conf["pick"]["files"]["log"]) == self.json_conf["pick"]["files"]["log"]:
                fname = path_join(self.json_conf["pick"]["files"]["output_dir"],
                                  self.json_conf["pick"]["files"]["log"])
            else:
                fname = self.json_conf["pick"]["files"]["log"]

            self.log_handler = logging.FileHandler(filename=fname, mode='w')
            assert os.path.exists(fname)

        # For the main logger I want to keep it at the "INFO" level
        self.log_level = self.json_conf["log_settings"]["log_level"]
        self.log_handler.setFormatter(self.formatter)
        self.logger.setLevel(self.log_level)
        self.logger.addHandler(self.log_handler)

        if self.log_level == "DEBUG" and self.json_conf["threads"] > 1:
            self.main_logger.setLevel(logging.DEBUG)
            self.main_logger.warning(
                    "Due to a Python design bug, we have to force Mikado to go in single-threaded mode when debugging.")
            self.procs = self.json_conf["threads"] = 1
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

        self.log_writer = logging_handlers.QueueListener(
            self.logging_queue, self.logger)
        self.log_writer.start()

        self.logger_queue_handler = logging_handlers.QueueHandler(self.logging_queue)
        self.queue_logger = logging.getLogger("parser")
        self.queue_logger.addHandler(self.logger_queue_handler)

        self.queue_logger.setLevel(logging.getLevelName(self.json_conf["log_settings"]["log_level"]))
        self.logger.warning("Current level for queue: %s", logging.getLevelName(self.queue_logger.level))

        self.queue_logger.propagate = False

        # Configure SQL logging
        sqllogger = logging.getLogger("sqlalchemy.engine")
        sqllogger.setLevel(self.json_conf["log_settings"]["sql_level"])
        sqllogger.addHandler(self.logger_queue_handler)

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
        try:
            for chrom in session.query(Chrom).order_by(Chrom.name.asc()):
                print("##sequence-region {0} 1 {1}".format(chrom.name, chrom.length),
                      file=locus_out)
                locus_out.flush()
        except sqlalchemy.exc.OperationalError as _:
            self.logger.error("Empty database! Creating a mock one")
            self.json_conf["db_settings"]["dbtype"] = "sqlite"
            self.json_conf["db_settings"]["db"] = tempfile.mktemp()
            engine = create_engine("{0}://".format(self.json_conf["db_settings"]["dbtype"]),
                                   creator=self.db_connection)
            session = sqlalchemy.orm.sessionmaker(bind=engine)()
            dbutils.DBBASE.metadata.create_all(engine)

        metrics = Superlocus.available_metrics[5:]
        metrics.extend(["external.{}".format(_.source) for _ in session.query(ExternalSource.source).all()])
        metrics = Superlocus.available_metrics[:5] + sorted(metrics)

        if self.sub_out != '':
            assert isinstance(self.sub_out, str)
            sub_metrics_file = open(re.sub("$", ".metrics.tsv",
                                    re.sub(".gff.?$", "", self.sub_out)), "w")
            sub_scores_file = open(re.sub("$", ".scores.tsv",
                                   re.sub(".gff.?$", "", self.sub_out)), "w")
            sub_metrics = csv.DictWriter(
                sub_metrics_file,
                metrics,
                extrasaction="ignore",
                delimiter="\t")
            sub_metrics.writeheader()
            sub_scores = csv.DictWriter(
                sub_scores_file, score_keys, delimiter="\t", extrasaction="ignore",)
            sub_scores.writeheader()
            # Not a very clean wya to do things ... attaching the handles as properties
            sub_scores.handle = sub_scores_file
            sub_scores.flush = sub_scores.handle.flush
            sub_scores.flush()
            sub_scores.close = sub_scores.handle.close
            sub_scores.name = sub_scores.handle.name
            sub_metrics.handle = sub_metrics_file
            sub_metrics.flush = sub_metrics.handle.flush
            sub_metrics.flush()
            sub_metrics.close = sub_metrics.handle.close
            sub_metrics.name = sub_metrics.handle.name
            sub_out = open(self.sub_out, 'w')
            print('##gff-version 3', file=sub_out)
            for chrom in session.query(Chrom).order_by(Chrom.name.desc()):
                print("##sequence-region {0} 1 {1}".format(chrom.name, chrom.length),
                      file=sub_out)
            sub_out.flush()
            sub_files = [sub_metrics,
                         sub_scores,
                         sub_out]
        else:
            sub_files = [None, None, None]

        if self.monolocus_out != '':
            mono_metrics_file = open(re.sub("$", ".metrics.tsv",
                                     re.sub(".gff.?$", "", self.monolocus_out)), "w")
            mono_scores_file = open(re.sub("$", ".scores.tsv",
                                    re.sub(".gff.?$", "", self.monolocus_out)), "w")          
            mono_metrics = csv.DictWriter(
                mono_metrics_file,
                metrics,
                extrasaction="ignore",
                delimiter="\t")
            mono_metrics.writeheader()
            mono_scores = csv.DictWriter(
                mono_scores_file, score_keys, delimiter="\t",
                extrasaction="ignore")
            mono_scores.writeheader()
            # Not a very clean wya to do things ... attaching the handles as properties
            mono_scores.handle = mono_scores_file
            mono_scores.flush = mono_scores.handle.flush
            mono_scores.flush()
            mono_scores.close = mono_scores.handle.close
            mono_scores.name = mono_scores.handle.name
            mono_metrics.handle = mono_metrics_file
            mono_metrics.flush = mono_metrics.handle.flush
            mono_metrics.flush()
            mono_metrics.close = mono_metrics.handle.close
            mono_metrics.name = mono_metrics.handle.name
            mono_out = open(self.monolocus_out, 'w')
            print('##gff-version 3', file=mono_out)
            for chrom in session.query(Chrom).order_by(Chrom.name.desc()):
                print("##sequence-region {0} 1 {1}".format(chrom.name, chrom.length),
                      file=mono_out)

            mono_out.flush()
            mono_files = [mono_metrics,
                          mono_scores,
                          mono_out]
        else:
            mono_files = [None] * 3

        session.close()
        engine.dispose()
        return sub_files, mono_files

    def __get_output_files(self):

        """
        Private method used by printer to prepare all the output files.
        :return: ((locus_metrics, locus_scores, locus_out),
                (sub_metrics, sub_scores, sub_out),
                mono_out)
                all of these are file handles
        """

        engine = create_engine("{0}://".format(self.json_conf["db_settings"]["dbtype"]),
                               creator=self.db_connection)
        session = sqlalchemy.orm.sessionmaker(bind=engine)()

        external_metrics = ["external.{}".format(_.source) for _ in session.query(ExternalSource.source).all()]

        score_keys = ["source_score"]
        if self.regressor is None:
            __scores = sorted(list(self.json_conf["scoring"].keys()))
            # Check that the external scores are all present. If they are not, raise a warning.
            __externals = set([_ for _ in __scores if _.startswith("external.")])
            if __externals - set(external_metrics):
                self.logger.error(
                    ("The following external metrics, found in the scoring file, are not present in the database. " +
                     "Please check their existence:\n" + "\n".join(
                                ["    - {metric}".format(metric=metric) for metric in sorted(
                                    __externals - set(external_metrics))]
                                ))
                )
                sys.exit(1)
                # __scores = sorted(set(__scores) - (__externals - set(external_metrics)))

            score_keys += __scores
        else:
            score_keys += self.regressor["scoring"].metrics

        score_keys = ["tid", "alias", "parent", "score"] + sorted(score_keys)
        # Define mandatory output files
        locus_metrics_file = open(re.sub("$", ".metrics.tsv", re.sub(
            ".gff.?$", "", self.locus_out)), "w")
        locus_scores_file = open(re.sub("$", ".scores.tsv", re.sub(
            ".gff.?$", "", self.locus_out)), "w")

        metrics = Superlocus.available_metrics[5:]
        metrics.extend(external_metrics)
        metrics = Superlocus.available_metrics[:5] + sorted(metrics)
        session.close()
        engine.dispose()

        locus_metrics = csv.DictWriter(
            locus_metrics_file,
            metrics,
            extrasaction="ignore",
            delimiter="\t")
        locus_metrics.writeheader()
        locus_scores = csv.DictWriter(locus_scores_file, score_keys, delimiter="\t", extrasaction="ignore",)
        locus_scores.writeheader()

        locus_metrics.handle = locus_metrics_file
        locus_metrics.flush = locus_metrics.handle.flush
        locus_metrics.close = locus_metrics.handle.close
        locus_metrics.name = locus_metrics.handle.name
        locus_scores.handle = locus_scores_file
        locus_scores.flush = locus_scores.handle.flush
        locus_scores.close = locus_scores.handle.close
        locus_scores.name = locus_scores.handle.name

        locus_out = open(self.locus_out, 'w')
        sub_files, mono_files = self.__print_gff_headers(locus_out, score_keys)

        return ((locus_metrics,
                 locus_scores,
                 locus_out),
                sub_files, mono_files)

    # This method has many local variables, but most (9!) are
    # actually file handlers. I cannot trim them down for now.
    # pylint: disable=too-many-locals

    def _submit_locus(self, slocus, counter, data_dict=None, engine=None):
        """
        Private method to submit / start the analysis of a superlocus in input.
        :param slocus: the locus to analyse.
        :param data_dict: the preloaded data in memory
        :param engine: connection engine
        :return: job object / None
        """

        if slocus is not None:
            self.main_logger.debug("Submitting %s", slocus.id)
        else:
            return []
            
        self.logger.debug("Loading data for %s", slocus.id)
        slocus.logger = self.logger
        return analyse_locus(slocus=slocus,
                             counter=counter,
                             json_conf=self.json_conf,
                             logging_queue=self.logging_queue,
                             data_dict=data_dict,
                             engine=engine)

    def __unsorted_interrupt(self, row, current_transcript):
        """
        Private method that brings the program to a screeching halt
         if the GTF/GFF is not properly sorted.
        :param row:
        :param current_transcript:
        :return:
        """

        # self.result_dict[counter] = [_]
        self.printer_queue.put_nowait(("EXIT", None))
        # self.printer_queue.put(("EXIT", counter + 1))
        current = "\t".join([str(x) for x in row])
        previous = "\t".join([str(x) for x in [current_transcript.chrom,
                                               current_transcript.start,
                                               current_transcript.end]])

        error_msg = """Unsorted input file, the results will not be correct.
            Please provide a properly sorted input. Error:
            {0}
            {1}""".format(current,
                          previous)
        self.logger.critical(error_msg)
        error_msg = "CRITICAL - {0}".format(" ".join(
            [l.strip() for l in error_msg.split("\n")]))
        raise UnsortedInput(error_msg)

    def __test_sortedness(self, row_coords, coords):
        """
        Private method to test whether a row and the current transcript are actually in the expected
        sorted order.
        :param row_coords:
        :param coords:
        :return:
        """
        test = True
        if coords is None:
            return test

        if not hasattr(coords, "chrom"):
            chrom, start, end = coords
        else:
            chrom, start, end = coords.chrom, coords.start, coords.end

        if not hasattr(row_coords, "chrom"):
            row_chrom, row_start, row_end = coords
        else:
            row_chrom, row_start, row_end = row_coords.chrom, row_coords.start, row_coords.end

        if any(_ is None for _ in (chrom, start, end)):
            return True
        if chrom > row_chrom:
            test = False
        elif (chrom == row_chrom and
              start > row_start):
            test = False
        elif (chrom == row_chrom and
              start == row_start and
              end > row_end):
            test = False

        if test is False:
            self.__unsorted_interrupt(row_coords, coords)

    @staticmethod
    def add_to_index(transcripts: dict,
                     counter: int,
                     locus_queue,
                     mapper: dict,):

        """Method to create the simple indexed database for features."""

        chroms = set([_["chrom"] for _ in transcripts.values()])
        if len(chroms) > 1:
            raise AssertionError(chroms)

        transcripts = msgpack.dumps([_ for _ in transcripts.values()])
        # cursor.execute("INSERT INTO transcripts VALUES (?, ?)", (counter, transcripts))
        # max_submit = 1000
        chrom = chroms.pop()
        if "done" not in mapper:
            mapper["done"] = set()
        if "submit" not in mapper:
            mapper["submit"] = set()

        if chrom not in mapper:
            mapper[chrom] = dict()
            mapper[chrom]["submit"] = set()
            mapper[chrom]["done"] = set()

        mapper[chrom]["submit"].add(counter)
        mapper["submit"].add(counter)
        mapper[counter] = chrom
        locus_queue.put((counter, transcripts))

        return mapper

    def __submit_multi_threading(self):

        """
        Method to execute Mikado pick in multi threaded mode.
        :return:
        """

        intron_range = self.json_conf["pick"]["run_options"]["intron_range"]
        self.logger.debug("Intron range: %s", intron_range)

        locus_queue = self.manager.JoinableQueue(-1)
        status_queue = self.manager.JoinableQueue(-1)

        handles = list(self.__get_output_files())
        if self.json_conf["pick"]["run_options"]["shm"] is True:
            basetempdir = "/dev/shm"
        else:
            basetempdir = self.json_conf["pick"]["files"]["output_dir"]

        tempdirectory = tempfile.TemporaryDirectory(suffix="",
                                                    prefix="mikado_pick_tmp",
                                                    dir=basetempdir)
        tempdir = tempdirectory.name
        self.logger.info("Starting Mikado with multiple processes, temporary directory:\n\t%s",
                         tempdir)

        self.logger.debug("Creating the worker processes")

        working_processes = [LociProcesser(msgpack.dumps(self.json_conf),
                                           locus_queue,
                                           self.logging_queue,
                                           status_queue,
                                           _,
                                           tempdir)
                             for _ in range(1, self.procs+1)]
        # Start all processes
        [_.start() for _ in working_processes]

        self.logger.debug("Started all %d workers", self.procs)
        # No sense in keeping this data available on the main thread now

        try:
            mapper = self.__parse_multithreaded(locus_queue)
        except UnsortedInput:
            [_.terminate() for _ in working_processes]
            raise

        self.logger.debug("Joining children processes")

        percs = percentile(range(1, max(mapper["submit"]) + 1),
                           range(10, 101, 10))
        curr_perc = 0
        total = len(mapper["submit"])
        mapper["results"] = dict()

        while mapper["done"] != mapper["submit"]:
            counter, chrom, num_genes, loci, subloci, monoloci = status_queue.get()
            mapper["done"].add(counter)
            if counter in mapper["results"]:
                self.logger.fatal("%d double index found!", counter)
                raise KeyError

            mapper["results"][counter] = (chrom, num_genes, loci, subloci, monoloci)
            if len(mapper["done"]) > percs[curr_perc]:
                curr_perc += 1
                while len(mapper["done"]) > percs[curr_perc]:
                    curr_perc += 1
                real_perc = round(len(mapper["done"]) * 100 / total)
                self.logger.info("Done %s%% of loci (%s out of %s)", real_perc,
                                 len(mapper["done"]), total)
            chrom = mapper[counter]
            mapper[chrom]["done"].add(counter)
            if mapper[chrom]["done"] == mapper[chrom]["submit"]:
                self.logger.info("Finished with chromosome %s", chrom)

        [_.join() for _ in working_processes]
        self.logger.info("Joined children processes; starting to merge partial files")

        # Merge loci
        merge_loci(mapper,
                   handles,
                   total,
                   logger=self.logger,
                   source=self.json_conf["pick"]["output_format"]["source"])

        self.logger.info("Finished merging partial files")
        try:
            self.logger.debug("Cleaning up the temporary directory")
            tempdirectory.cleanup()
            self.logger.debug("Finished cleaning up")
            pass
        except (OSError, FileNotFoundError, FileExistsError) as exc:
            self.logger.warning("Failed to clean up the temporary directory %s, error: %s",
                                tempdir, exc)
        except KeyboardInterrupt:
            raise
        except Exception as exc:
            self.logger.exception("Failed to clean up the temporary directory %s, error: %s", exc)
        finally:
            return

    def __check_max_intron(self, current, invalids, row, max_intron):
        previous = None
        if hasattr(row, "feature"):
            chrom, feature, start, end, phase, tid = (row.chrom, row.feature, row.start,
                                                      row.end, row.phase, row.transcript)
        else:
            # return line, chrom, start, end, tid, is_transcript
            _, chrom, feature, start, end, phase, tid, _ = row

        for exon in reversed(current["transcripts"][tid]["exon_lines"]):
            if exon[2] == feature:
                previous = exon
                break

        current["transcripts"][tid]["exon_lines"].append((start, end, feature, phase))
        if previous:
            # I have to compare like with like.
            intron_length = (start - 1) - (previous[1] + 1) + 1
            if intron_length >= max_intron:
                self.logger.warning(
                    "%s has an intron (%s) greater than the maximum allowed (%s). Ignoring it.",
                    tid, intron_length, max_intron
                )
                del current["transcripts"][tid]
                invalids.add(tid)
                if current["transcripts"]:
                    current["start"] = min([current["transcripts"][trans]["start"]
                                            for trans in current["transcripts"]])
                    current["end"] = max([current["transcripts"][trans]["end"]
                                          for trans in current["transcripts"]])
                else:
                    current["start"], current["end"] = None, None
        return current, invalids

    transcript_gtf_pattern = re.compile(r"""transcript_id "([^"]*)\"""")

    def _parse_gtf_line(self, line):

        if not line or line[0] == "#":
            return None
        fields = line.split("\t")
        if len(fields) != 9:
            return None

        try:
            start, end = int(fields[3]), int(fields[4])
        except (ValueError, SystemError, TypeError):
            return None
        chrom = fields[0]
        is_exon = (GtfLine.exon_pattern.search(fields[2]) is not None)
        is_transcript = False
        if not is_exon:
            is_transcript = (GtfLine.transcript_pattern.search(fields[2]) is not None)

        if not (is_exon or is_transcript):
            return None

        tid = self.transcript_gtf_pattern.search(fields[-1])
        if tid is None:
            raise InvalidJson("Corrupt input GTF file, offending line:\n{}".format(line))
        tid = tid.groups()[0]
        if fields[7] in (None, ".", "?"):
            phase = None
        else:
            try:
                phase = int(fields[7])
            except (SystemError, TypeError, ValueError):
                return None
        return line, chrom, fields[2], start, end, phase, tid, is_transcript

    def __parse_multithreaded(self, locus_queue):
        counter = 0
        invalids = set()
        flank = self.json_conf["pick"]["clustering"]["flank"]
        mapper = dict()

        with self.define_input(multithreading=True) as input_annotation:
            current = {"chrom": None, "start": None, "end": None, "transcripts": dict()}
            max_intron = self.json_conf["prepare"]["max_intron_length"]
            for row in input_annotation:
                row = self._parse_gtf_line(row)
                if row is None:  # Header
                    continue
                line, chrom, feature, start, end, phase, tid, is_transcript = row
                if self.__regions and chrom not in self.__regions:
                    continue

                if is_transcript is False:
                    if tid in invalids:
                        continue
                    elif tid not in current["transcripts"]:
                        self.logger.fatal("Transcript %s is invalid", tid)
                        raise UnsortedInput
                        # invalids.add(tid)
                    else:
                        # Check max intron length. We presume that exons are sorted correctly.
                        current, invalids = self.__check_max_intron(current, invalids,
                                                                    row, max_intron=max_intron)
                elif is_transcript is True:
                    self.__test_sortedness((chrom, start, end),
                                           (current["chrom"], current["start"], current["end"]))
                    if self.__regions:
                        if not [_ for _ in self.__regions[chrom].find(start, end) if _.start <= start and end <= _.end]:
                            invalids.add(tid)
                            continue

                    if current["chrom"] != chrom or not current["start"]:
                        if current["chrom"] != chrom:
                            if current["chrom"] is not None and current["chrom"] != chrom:
                                self.logger.debug("Finished chromosome %s", current["chrom"])
                                if len(current["transcripts"]) > 0:
                                    self.logger.debug("Submitting locus # %d (%s), with transcripts:\n%s",
                                                      counter, "{}:{}-{}".format(current["chrom"],
                                                                                 current["start"], current["end"]),
                                                      ",".join(list(current["transcripts"].keys())))
                                    counter += 1
                                    mapper = self.add_to_index(current["transcripts"], counter, locus_queue, mapper)
                            self.logger.debug("Starting chromosome %s", chrom)

                        current["chrom"], current["start"], current["end"] = chrom, start, end
                        current["transcripts"] = dict()
                        current["transcripts"][tid] = dict()
                        current["transcripts"][tid]["chrom"] = chrom
                        current["transcripts"][tid]["start"] = start
                        current["transcripts"][tid]["end"] = end
                        current["transcripts"][tid]["definition"] = line
                        current["transcripts"][tid]["exon_lines"] = []
                    else:
                        if Superlocus.overlap((current["end"], current["start"]),
                                              (start, end), flank=flank) > 0:
                            # Add to the locus!
                            current["start"] = min(current["start"], start)
                            current["end"] = max(current["end"], end)
                        elif len(current["transcripts"]) > 0:
                            counter += 1
                            self.logger.debug("Submitting locus # %d (%s), with transcripts:\n%s",
                                              counter, "{}:{}-{}".format(current["chrom"],
                                                                         current["start"], current["end"]),
                                              ",".join(list(current["transcripts"].keys())))
                            mapper = self.add_to_index(current["transcripts"], counter, locus_queue, mapper)
                            current["start"], current["end"] = start, end
                            current["transcripts"] = dict()

                        current["transcripts"][tid] = dict()
                        current["transcripts"][tid]["chrom"] = chrom
                        current["transcripts"][tid]["definition"] = line
                        current["transcripts"][tid]["exon_lines"] = []
                        current["transcripts"][tid]["start"] = start
                        current["transcripts"][tid]["end"] = end

        if current["start"] is not None:
            if len(current["transcripts"]) > 0:
                counter += 1
                self.logger.debug("Submitting locus # %d (%s), with transcripts:\n%s",
                                  counter, "{}:{}-{}".format(current["chrom"],
                                                             current["start"], current["end"]),
                                  ",".join(list(current["transcripts"].keys())))
                mapper = self.add_to_index(current["transcripts"], counter, locus_queue, mapper)
            self.logger.debug("Finished chromosome %s", current["chrom"])

        locus_queue.put(("EXIT", None))

        return mapper

    def __submit_single_threaded(self):

        """
        Method to execute Mikado pick in single threaded mode.
        :return:
        """

        current_locus = None
        current_transcript = None

        handler = logging_handlers.QueueHandler(self.logging_queue)
        logger = logging.getLogger("queue_listener")
        logger.propagate = False
        logger.addHandler(handler)
        logger.setLevel(self.json_conf["log_settings"]["log_level"])
        logger.debug("Begun single-threaded run")

        intron_range = self.json_conf["pick"]["run_options"]["intron_range"]
        logger.debug("Intron range: %s", intron_range)

        handles = self.__get_output_files()

        locus_printer = functools.partial(print_locus,
                                          handles=handles,
                                          logger=logger,
                                          json_conf=self.json_conf)

        # last_printed = -1
        curr_chrom = None
        gene_counter = 0

        self.engine = dbutils.connect(json_conf=self.json_conf, logger=self.logger)

        submit_locus = functools.partial(self._submit_locus, **{"data_dict": None,
                                                                "engine": self.engine})

        counter = -1
        invalid = False
        max_intron = self.json_conf["prepare"]["max_intron_length"]
        skip_transcript = False
        with self.define_input() as input_annotation:
            for row in input_annotation:
                if row.is_exon is True and invalid is False and skip_transcript is False:
                    try:
                        current_transcript.add_exon(row)
                    except InvalidTranscript as exc:
                        self.logger.error("Transcript %s is invalid;\n%s",
                                          current_transcript.id,
                                          exc)
                        invalid = True
                elif row.is_transcript is True:
                    self.__test_sortedness(row, current_transcript)
                    current_locus, counter, gene_counter, curr_chrom = self.__check_transcript(
                        current_transcript, current_locus, counter, max_intron,
                        gene_counter, curr_chrom, locus_printer, submit_locus)
                    invalid = False
                    if self.__regions:
                        if row.chrom not in self.__regions:
                            skip_transcript = True
                            current_transcript = None
                        elif not [_ for _ in self.__regions[row.chrom].find(row.start, row.end)
                                  if _.start <= row.start and row.end <= _.end]:
                            skip_transcript = True
                            current_transcript = None
                        else:
                            current_transcript = Transcript(row, intron_range=intron_range, logger=logger)
                            skip_transcript = False
                    else:
                        current_transcript = Transcript(row, intron_range=intron_range, logger=logger)
                        skip_transcript = False
                    if skip_transcript is False and (current_transcript is None or row.chrom != current_transcript.chrom):
                        if current_transcript is not None:
                            self.logger.info("Finished chromosome %s",
                                             current_transcript.chrom)
                        self.logger.info("Starting chromosome %s", row.chrom)

        self.logger.info("Finished chromosome %s", current_locus.chrom)
        current_locus, counter, gene_counter, curr_chrom = self.__check_transcript(
            current_transcript, current_locus, counter, max_intron,
            gene_counter, curr_chrom, locus_printer, submit_locus)
        counter += 1
        self.logger.debug("Analysing locus # %d", counter)
        for stranded_locus in submit_locus(current_locus, counter):
            if stranded_locus.chrom != curr_chrom:
                curr_chrom = stranded_locus.chrom
                gene_counter = 0
            gene_counter = locus_printer(stranded_locus, gene_counter)
        # submit_locus(current_locus, counter)
        for group in handles:
            [_.close() for _ in group if _]
        logger.info("Final number of superloci: %d", counter)

    def __check_transcript(self, current_transcript, current_locus, counter, max_intron,
                           gene_counter, curr_chrom, locus_printer, submit_locus):

        if current_transcript is not None:
            current_transcript.finalize()
            if current_transcript.max_intron_length > max_intron:
                self.logger.warning(
                    "%s has an intron (%s) which is greater than the maximum allowed (%s). Removing it.",
                    current_transcript.id, current_transcript.max_intron_length, max_intron)
            elif current_locus and Superlocus.in_locus(
                            current_locus, current_transcript,
                            flank=self.json_conf["pick"]["clustering"]["flank"]) is True:
                    current_locus.add_transcript_to_locus(current_transcript, check_in_locus=False)
            else:
                counter += 1
                self.logger.debug("Analysing locus # %d", counter)
                try:
                    for stranded_locus in submit_locus(current_locus, counter):
                        if stranded_locus.chrom != curr_chrom:
                            curr_chrom = stranded_locus.chrom
                            gene_counter = 0
                        gene_counter = locus_printer(stranded_locus, gene_counter)
                except KeyboardInterrupt:
                    raise
                except Exception as exc:
                    self.logger.exception(
                        "Superlocus %s failed with exception: %s",
                        None if current_locus is None else current_locus.id, exc)

                current_locus = Superlocus(current_transcript, stranded=False, json_conf=self.json_conf,
                                           source=self.json_conf["pick"]["output_format"]["source"])
                if self.regressor is not None:
                    current_locus.regressor = self.regressor

        return current_locus, counter, gene_counter, curr_chrom

    def _parse_and_submit_input(self):

        """
        This method does the parsing of the input and submission of the loci to the
        _submit_locus method.
        :return: jobs (the list of all jobs already submitted)
        """

        single_thread = (self.json_conf["pick"]["run_options"]["single_thread"] or self.procs == 1
                         or self.json_conf["log_settings"]["log_level"] == "DEBUG")
        
        if single_thread is False:
            self.__submit_multi_threading()
        else:
            self.__submit_single_threaded()
        return

    def __cleanup(self):
        self.log_writer.stop()
        if self.queue_pool is not None:
            self.queue_pool.dispose()

        # Clean up the DB copied to SHM
        if self.json_conf["pick"]["run_options"]["shm"] is True:
            assert os.path.dirname(self.json_conf["db_settings"]["db"]) == os.path.join("/dev", "shm"), (
                self.json_conf["db_settings"]["db"], os.path.dirname(self.json_conf["db_settings"]["db"]),
                os.path.join("/dev", "shm"))
            self.main_logger.debug("Removing shared memory DB %s", self.json_conf["db_settings"]["db"])
            os.remove(self.json_conf["db_settings"]["db"])
        self.manager.shutdown()

    def __call__(self):

        """This method will activate the class and start the analysis of the input file."""

        # NOTE: Pool, Process and Manager must NOT become instance attributes!
        # Otherwise it will raise all sorts of mistakes

        self.logger.debug("Source: %s",
                          self.json_conf["pick"]["output_format"]["source"])
        if self.json_conf["db_settings"]["dbtype"] == "sqlite":
            self.queue_pool = sqlalchemy.pool.QueuePool(
                self.db_connection,
                pool_size=self.procs,
                max_overflow=0)

        try:
            self._parse_and_submit_input()
        except UnsortedInput as _:
            self.logger.error(
                "The input files were not properly sorted! Please run prepare and retry.")
        except Exception:
            self.__cleanup()
            raise

        self.__cleanup()
        self.main_logger.info("Finished analysis of %s", self.input_file)

        sys.exit(0)
# pylint: enable=too-many-instance-attributes
