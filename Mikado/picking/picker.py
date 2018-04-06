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
import logging
from logging import handlers as logging_handlers
import collections
import functools
import multiprocessing
from sqlalchemy.engine import create_engine  # SQLAlchemy/DB imports
from sqlalchemy.orm.session import sessionmaker
from sqlalchemy.pool import QueuePool as SqlPool
import sqlalchemy
import sqlalchemy.exc
from ..utilities import path_join
from ..utilities.log_utils import formatter
from ..parsers.GTF import GTF
from ..parsers.GFF import GFF3
from ..serializers.blast_serializer import Hit, Query
from ..serializers.junction import Junction, Chrom
from ..serializers.orf import Orf
from ..serializers.external import ExternalSource
from ..loci.superlocus import Superlocus, Transcript
from ..configuration.configurator import to_json, check_json  # Necessary for nosetests
from ..utilities import dbutils, merge_partial
from ..exceptions import UnsortedInput, InvalidJson, InvalidTranscript
from .loci_processer import analyse_locus, LociProcesser, merge_loci, print_locus
import multiprocessing.managers
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
import pickle
import warnings
# from math import floor
logging.captureWarnings(True)
warnings.simplefilter("always")
import sqlite3
try:
    import ujson as json
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
    def __init__(self, json_conf, commandline=""):

        """Constructor. It takes a single argument as input - the JSON/YAML configuration,
        prepared by the json_utils functions.

        :param json_conf: Either a configuration dictionary or the configuration file.

        :param commandline: optional, the commandline used to start the program
        :type commandline: str
        """

        # Mock variables
        self.formatter = self.main_logger = self.log_writer = self.log_handler = self.logger = None
        self.log_level = "WARN"

        # Things that have to be deleted upon serialisation
        self.not_pickable = ["queue_logger", "manager", "printer_process",
                             "log_process", "pool", "main_logger",
                             "log_handler", "log_writer", "logger", "engine"]

        # Now we start the real work
        self.commandline = commandline
        self.json_conf = json_conf

        self.__load_configuration()
        self.regressor = None

        self.procs = self.json_conf["pick"]["run_options"]["procs"]

        # Check the input file
        with self.define_input() as _:
            pass

        self.__create_output_handles()
        # pylint: disable=no-member
        multiprocessing.set_start_method(self.json_conf["multiprocessing_method"],
                                         force=True)

        # self.setup_logger()
        self.logger.info("Multiprocessing method: %s",
                         self.json_conf["multiprocessing_method"])

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
            if self.json_conf["pick"]["run_options"]["procs"] > 1:
                self.main_logger.warning("Reset number of threads to 1 as requested")
                self.procs = 1
        elif self.json_conf["pick"]["run_options"]["procs"] == 1:
            self.json_conf["pick"]["run_options"]["single_thread"] = True

        if self.json_conf["pick"]["run_options"]["preload"] is True and self.procs > 1:
            self.logger.warning(
                "Preloading using multiple threads can be extremely \
memory intensive, proceed with caution!")

        if self.locus_out is None:
            raise InvalidJson(
                "No output prefix specified for the final loci. Key: \"loci_out\"")

        # self.printer_process = Process(target=self.printer)
        self.queue_pool = None

    def define_input(self):
        """Function to check that the input file exists and is valid. It returns the parser."""

        if self.input_file.endswith(".gtf"):
            parser = GTF
        else:
            parser = GFF3

        verified = False
        with parser(self.input_file) as testing:
            for row in testing:
                if row.header is False:
                    verified = True
                    break
        if verified is False:
            raise InvalidJson(
                "Invalid input file: {0}".format(self.input_file))

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
            self.json_conf = check_json(self.json_conf, logger=self.logger)
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
                    """The \"{}\" property has now been moved to pick/{}. Please update your configuration files in the future.""".format(
                        key, new_home))
                self.logger.warn(warns)

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

        self.log_writer = logging_handlers.QueueListener(
            self.logging_queue, self.logger)
        self.log_writer.start()

        self.logger_queue_handler = logging_handlers.QueueHandler(self.logging_queue)
        self.queue_logger = logging.getLogger("parser")
        self.queue_logger.addHandler(self.logger_queue_handler)

        self.queue_logger.setLevel(logging.getLevelName(self.json_conf["log_settings"]["log_level"]))
        self.logger.warn("Current level for queue: %s", logging.getLevelName(self.queue_logger.level))

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

        metrics = Superlocus.available_metrics[3:]
        metrics.extend(["external.{}".format(_.source) for _ in session.query(ExternalSource.source).all()])
        metrics = Superlocus.available_metrics[:3] + sorted(metrics)

        if self.sub_out != '':
            assert isinstance(self.sub_out, str)
            sub_metrics_file = open(re.sub("$", ".metrics.tsv",
                                    re.sub(".gff.?$", "", self.sub_out)), "w")
            sub_scores_file = open(re.sub("$", ".scores.tsv",
                                   re.sub(".gff.?$", "", self.sub_out)), "w")
            sub_metrics = csv.DictWriter(
                sub_metrics_file,
                metrics,
                delimiter="\t")
            sub_metrics.writeheader()
            sub_scores = csv.DictWriter(
                sub_scores_file, score_keys, delimiter="\t")
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
                delimiter="\t")
            mono_metrics.writeheader()
            mono_scores = csv.DictWriter(
                mono_scores_file, score_keys, delimiter="\t")
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

        score_keys = ["source_score"]
        if self.regressor is None:
            score_keys += sorted(list(self.json_conf["scoring"].keys()))
        else:
            score_keys += self.regressor["scoring"].metrics

        score_keys = ["tid", "alias", "parent", "score"] + sorted(score_keys)
        # Define mandatory output files
        locus_metrics_file = open(re.sub("$", ".metrics.tsv", re.sub(
            ".gff.?$", "", self.locus_out)), "w")
        locus_scores_file = open(re.sub("$", ".scores.tsv", re.sub(
            ".gff.?$", "", self.locus_out)), "w")

        engine = create_engine("{0}://".format(self.json_conf["db_settings"]["dbtype"]),
                               creator=self.db_connection)
        session = sqlalchemy.orm.sessionmaker(bind=engine)()

        metrics = Superlocus.available_metrics[4:]
        metrics.extend(["external.{}".format(_.source) for _ in session.query(ExternalSource.source).all()])
        metrics = Superlocus.available_metrics[:4] + sorted(metrics)
        session.close()
        engine.dispose()

        locus_metrics = csv.DictWriter(
            locus_metrics_file,
            metrics,
            delimiter="\t")
        locus_metrics.writeheader()
        locus_scores = csv.DictWriter(locus_scores_file, score_keys, delimiter="\t")
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

    def __getstate__(self):

        state = self.__dict__.copy()
        for not_pickable in self.not_pickable:
            if not_pickable in state:
                del state[not_pickable]

        print(state)
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
        for hsp in engine.execute("select * from hsp where hsp_evalue <= {0}".format(
            self.json_conf["pick"]["chimera_split"]["blast_params"]["hsp_evalue"]
        )):
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
        dbutils.DBBASE.metadata.create_all(self.engine)
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
        slocus.load_all_transcript_data(engine=engine,
                                        data_dict=data_dict)
        # slocus_id = slocus.id
        if slocus.initialized is False:
            # This happens when we have removed all transcripts from the locus
            # due to errors which should have been caught and logged
            self.logger.warning(
                "%s had all transcripts failing checks, ignoring it",
                slocus.id)
            # Exit
            return []

        return analyse_locus(slocus=slocus,
                             counter=counter,
                             json_conf=self.json_conf,
                             printer_queue=None,
                             logging_queue=self.logging_queue,
                             data_dict=None,
                             engine=None)

    def __unsorted_interrupt(self, row, current_transcript):
        """
        Private method that brings the program to a screeching halt
         if the GTF/GFF is not properly sorted.
        :param row:
        :param current_transcript:
        :return:
        """

        # self.result_dict[counter] = [_]
        self.printer_queue.put_nowait(("EXIT", ))
        # self.printer_queue.put(("EXIT", counter + 1))
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

    def __test_sortedness(self, row, current_transcript):
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
            self.__unsorted_interrupt(row, current_transcript)

    @staticmethod
    def add_to_index(conn: sqlite3.Connection,
                     cursor: sqlite3.Cursor,
                     transcripts: [Transcript],
                     counter: int):

        """Method to create the simple indexed database for features."""

        transcripts = json.dumps([_.as_dict() for _ in transcripts])
        cursor.execute("INSERT INTO transcripts VALUES (?, ?)", (counter, transcripts))
        conn.commit()

        return

    @staticmethod
    def _create_temporary_store(tempdirectory):

        conn = sqlite3.connect(os.path.join(tempdirectory, "temp_store.db"))
        cursor = conn.cursor()
        # cursor.execute("PRAGMA journal_mode=wal")
        cursor.execute("CREATE TABLE transcripts (counter integer, json blob)")
        cursor.execute("CREATE INDEX tid_idx on transcripts(counter)")

        return conn, cursor

    def __submit_multi_threading(self, data_dict):

        """
        Method to execute Mikado pick in multi threaded mode.

        :param data_dict: The data dictionary
        :return:
        """

        intron_range = self.json_conf["pick"]["run_options"]["intron_range"]
        self.logger.info("Intron range: %s", intron_range)

        current_locus = None
        current_transcript = None

        locus_queue = multiprocessing.JoinableQueue(-1)

        handles = list(self.__get_output_files())
        [_.close() for _ in handles[0]]
        handles[0] = [_.name for _ in handles[0]]

        if handles[1][0] is not None:
            [_.close() for _ in handles[1]]
            handles[1] = [_.name for _ in handles[1]]
        if handles[2][0] is not None:
            [_.close() for _ in handles[2]]
            handles[2] = [_.name for _ in handles[2]]

        tempdirectory = tempfile.TemporaryDirectory(suffix="",
                                                    prefix="mikado_pick_tmp",
                                                    dir=self.json_conf["pick"]["files"]["output_dir"])
        tempdir = tempdirectory.name

        # tempdir = os.path.join(self.json_conf["pick"]["files"]["output_dir"], "mikado_pick_tmp")
        # os.makedirs(tempdir, exist_ok=True)

        self.logger.info("Creating the worker processes")
        conn, cursor = self._create_temporary_store(tempdir)
        working_processes = [LociProcesser(self.json_conf,
                                           data_dict,
                                           handles,
                                           locus_queue,
                                           self.logging_queue,
                                           _,
                                           tempdir)
                             for _ in range(1, self.procs+1)]
        # Start all processes
        [_.start() for _ in working_processes]
        self.logger.info("Started all %d workers", self.procs)
        # No sense in keeping this data available on the main thread now
        del data_dict

        counter = 0
        invalid = False
        with self.define_input() as input_annotation:
            for row in input_annotation:

                if row.is_exon is True and invalid is False:
                    try:
                        current_transcript.add_exon(row)
                    except InvalidTranscript as exc:
                        self.logger.error("Transcript %s is invalid;\n%s",
                                          current_transcript.id,
                                          exc)
                        invalid = True
                elif row.is_transcript is True:
                    if current_transcript is not None and invalid is False:
                        self.__test_sortedness(row, current_transcript)
                        if Superlocus.in_locus(
                                current_locus, current_transcript,
                                flank=self.json_conf["pick"]["clustering"]["flank"]) is True:
                            current_locus.add_transcript_to_locus(current_transcript,
                                                                  check_in_locus=False)
                        else:
                            if current_locus is not None:
                                counter += 1
                                self.logger.debug("Submitting locus # %d (%s), with transcripts:\n%s",
                                                  counter,
                                                  None if not current_locus else current_locus.id,
                                                  ",".join(list(current_locus.transcripts.keys())))
                                self.add_to_index(conn, cursor, current_locus.transcripts.values(), counter)
                                locus_queue.put((counter, ))
                            current_locus = Superlocus(
                                current_transcript,
                                stranded=False,
                                json_conf=self.json_conf,
                                source=self.json_conf["pick"]["output_format"]["source"])

                    if current_transcript is None or row.chrom != current_transcript.chrom:
                        if current_transcript is not None:
                            self.logger.info("Finished chromosome %s",
                                             current_transcript.chrom)
                        self.logger.info("Starting chromosome %s", row.chrom)

                    invalid = False
                    current_transcript = Transcript(
                        row,
                        intron_range=intron_range)

        if current_transcript is not None and invalid is False:
            if Superlocus.in_locus(
                    current_locus, current_transcript,
                    flank=self.json_conf["pick"]["clustering"]["flank"]) is True:
                current_locus.add_transcript_to_locus(
                    current_transcript, check_in_locus=False)
            else:
                if current_locus is not None:
                    counter += 1
                    self.logger.debug("Submitting locus #%d (%s)", counter,
                                      None if not current_locus else current_locus.id)
                    self.add_to_index(conn, cursor, current_locus.transcripts.values(), counter)
                    locus_queue.put((counter,))

                current_locus = Superlocus(
                    current_transcript,
                    stranded=False,
                    json_conf=self.json_conf,
                    source=self.json_conf["pick"]["output_format"]["source"])
                self.logger.debug("Created last locus %s",
                                  current_locus)
        elif invalid is True and current_locus is not None:
            counter += 1
            self.logger.debug("Submitting locus #%d (%s)", counter,
                              None if not current_locus else current_locus.id)
            self.add_to_index(conn, cursor, current_locus.transcripts.values(), counter)
            locus_queue.put((counter,))

        self.logger.info("Finished chromosome %s", current_locus.chrom)

        counter += 1
        self.add_to_index(conn, cursor, current_locus.transcripts.values(), counter)
        locus_queue.put((counter,))

        self.logger.debug("Submitting locus %s, counter %d, with transcripts:\n%s",
                          current_locus.id, counter,
                          ", ".join(list(current_locus.transcripts.keys())))
        locus_queue.put(("EXIT", ))
        self.logger.info("Joining children processes")
        [_.join() for _ in working_processes]
        conn.close()
        self.logger.info("Joined children processes; starting to merge partial files")

        # Merge loci
        merge_loci(self.procs,
                   handles[0],
                   prefix=self.json_conf["pick"]["output_format"]["id_prefix"],
                   tempdir=tempdir)

        for handle in handles[1]:
            if handle is not None:
                with open(handle, "a") as output:
                    partials = [os.path.join(tempdir,
                                             "{0}-{1}".format(os.path.basename(handle), _))
                                for _ in range(1, self.procs + 1)]
                    merge_partial(partials, output, logger=self.logger)

        for handle in handles[2]:
            if handle is not None:
                with open(handle, "a") as output:
                    partials = [os.path.join(tempdir,
                                             "{0}-{1}".format(os.path.basename(handle), _))
                                for _ in range(1, self.procs + 1)]
                    merge_partial(partials, output, logger=self.logger)

        self.logger.info("Finished merging partial files")
        try:
            # shutil.rmtree(tempdir)
            tempdirectory.cleanup()
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

    def __submit_single_threaded(self, data_dict):

        """
        Method to execute Mikado pick in single threaded mode.
        :param data_dict:
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
        logger.info("Intron range: %s", intron_range)

        handles = self.__get_output_files()

        locus_printer = functools.partial(print_locus,
                                          handles=handles,
                                          logger=logger,
                                          counter=None,
                                          json_conf=self.json_conf)

        # last_printed = -1
        curr_chrom = None
        gene_counter = 0

        if self.json_conf["pick"]["run_options"]["preload"] is False:
            self.engine = dbutils.connect(json_conf=self.json_conf, logger=self.logger)
        else:
            self.engine = None

        submit_locus = functools.partial(self._submit_locus, **{"data_dict": data_dict,
                                                                "engine": self.engine})

        counter = -1
        invalid = False
        with self.define_input() as input_annotation:
            for row in input_annotation:
                if row.is_exon is True and invalid is False:
                    try:
                        current_transcript.add_exon(row)
                    except InvalidTranscript as exc:
                        self.logger.error("Transcript %s is invalid;\n%s",
                                          current_transcript.id,
                                          exc)
                        invalid = True
                elif row.is_transcript is True:
                    if current_transcript is not None and invalid is False:
                        self.__test_sortedness(row, current_transcript)
                        if Superlocus.in_locus(
                                current_locus, current_transcript,
                                flank=self.json_conf["pick"]["clustering"]["flank"]) is True:
                            current_locus.add_transcript_to_locus(current_transcript,
                                                                  check_in_locus=False)
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

                            current_locus = Superlocus(
                                current_transcript,
                                stranded=False,
                                json_conf=self.json_conf,
                                source=self.json_conf["pick"]["output_format"]["source"])
                            if self.regressor is not None:
                                current_locus.regressor = self.regressor

                    if current_transcript is None or row.chrom != current_transcript.chrom:
                        if current_transcript is not None:
                            self.logger.info("Finished chromosome %s",
                                             current_transcript.chrom)
                        self.logger.info("Starting chromosome %s", row.chrom)

                    invalid = False
                    current_transcript = Transcript(
                        row,
                        intron_range=intron_range, logger=logger)

        if current_transcript is not None and invalid is False:
            if Superlocus.in_locus(
                    current_locus, current_transcript,
                    flank=self.json_conf["pick"]["clustering"]["flank"]) is True:
                current_locus.add_transcript_to_locus(
                    current_transcript, check_in_locus=False)
            else:
                counter += 1
                self.logger.debug("Analysing locus # %d", counter)
                for stranded_locus in submit_locus(current_locus, counter):
                    if stranded_locus.chrom != curr_chrom:
                        curr_chrom = stranded_locus.chrom
                        gene_counter = 0
                    gene_counter = locus_printer(stranded_locus, gene_counter)

                current_locus = Superlocus(
                    current_transcript,
                    stranded=False,
                    json_conf=self.json_conf,
                    source=self.json_conf["pick"]["output_format"]["source"])
                self.logger.debug("Created last locus %s",
                                  current_locus)
        elif current_transcript is not None and invalid is True:
            if current_locus is not None:
                counter += 1
                self.logger.debug("Analysing locus # %d", counter)
                for stranded_locus in submit_locus(current_locus, counter):
                    if stranded_locus.chrom != curr_chrom:
                        curr_chrom = stranded_locus.chrom
                        gene_counter = 0
                    gene_counter = locus_printer(stranded_locus, gene_counter)

        self.logger.info("Finished chromosome %s", current_locus.chrom)

        counter += 1
        self.logger.debug("Analysing locus # %d", counter)
        # if current_locus is not None:
        #     current_locus.load_all_transcript_data(pool=self.connection_pool,
        #                                            data_dict=data_dict)
        for stranded_locus in submit_locus(current_locus, counter):
            if stranded_locus.chrom != curr_chrom:
                curr_chrom = stranded_locus.chrom
                gene_counter = 0
            gene_counter = locus_printer(stranded_locus, gene_counter)
        # submit_locus(current_locus, counter)
        for group in handles:
            [_.close() for _ in group if _]
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
                pool_size=self.procs,
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

        sys.exit(0)
# pylint: enable=too-many-instance-attributes
