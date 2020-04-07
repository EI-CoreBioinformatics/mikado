"""
XML serialisation class.
"""

import os
import logging.handlers as logging_handlers
import logging
from sqlalchemy.orm.session import Session
from ...utilities.dbutils import DBBASE
import pysam
from ...utilities.dbutils import connect
from ...utilities.log_utils import create_null_logger, check_logger
from . import Query, Target, Hsp, Hit
from .utils import load_into_db
from .xml_utils import get_multipliers
from .xml_serialiser import _serialise_xmls
from .tab_serialiser import _serialise_tabular
import pandas as pd
import multiprocessing


__author__ = 'Luca Venturini'


# A serialisation class must have a ton of attributes ...
# pylint: disable=too-many-instance-attributes


class BlastSerializer:
    """This class has the role of taking in input a blast XML file and (partially)
    serialise it into a database. We are using SQLalchemy, so the database type
    could be any of SQLite, MySQL, PSQL, etc."""

    __name__ = "BlastSerializer"
    logger = create_null_logger(__name__)

    # This evaluates to characters A-z and |. Used to detect true matches
    # in the match line of blast alignments
    # Put here so it's calculated *once* inside the program
    __valid_matches = set([chr(x) for x in range(65, 91)] +
                          [chr(x) for x in range(97, 123)] +
                          ["|", "*"])

    def __init__(self, xml_name,
                 logger=None,
                 json_conf=None):
        """Initializing method. Arguments:

        :param xml_name: The XML(s) to parse.

        Arguments:

        :param json_conf: a configuration dictionary.
        :type json_conf: dict


        """

        if json_conf is None:
            raise ValueError("No configuration provided!")

        if logger is not None:
            self.logger = check_logger(logger)
        else:
            raise ValueError("No logger provided!")
        self.logger.debug("Started to serialise %s, log level: %s",
                         xml_name, self.logger.level)

        # Runtime arguments

        self.procs = json_conf["threads"]
        self.single_thread = json_conf["serialise"]["single_thread"]
        self.json_conf = json_conf
        # pylint: disable=unexpected-argument,E1123
        multiprocessing.set_start_method(self.json_conf["multiprocessing_method"],
                                         force=True)
        # pylint: enable=unexpected-argument,E1123
        self.logger.info("Number of dedicated workers: %d", self.procs)
        self.logging_queue = multiprocessing.Queue(-1)
        self.logger_queue_handler = logging_handlers.QueueHandler(self.logging_queue)
        self.queue_logger = logging.getLogger("parser")
        self.queue_logger.addHandler(self.logger_queue_handler)
        self.queue_logger.setLevel(self.json_conf["log_settings"]["log_level"])
        self.queue_logger.propagate = False
        self.log_writer = logging_handlers.QueueListener(self.logging_queue, self.logger)
        self.log_writer.start()

        self._max_target_seqs = json_conf["serialise"]["max_target_seqs"]
        self.maxobjects = json_conf["serialise"]["max_objects"]
        target_seqs = json_conf["serialise"]["files"]["blast_targets"]
        query_seqs = json_conf["serialise"]["files"]["transcripts"]

        self.header = None
        if xml_name is None:
            self.logger.warning("No BLAST XML provided. Exiting.")
            return

        self.engine = connect(json_conf)

        self._xml_debug = self.json_conf.get("serialise", dict()).get("blast_loading_debug", False)
        if self._xml_debug:
            self.logger.warning("Activating the XML debug mode")
            self.single_thread = True
            self.procs = 1

        # session = sessionmaker(autocommit=True)
        DBBASE.metadata.create_all(self.engine)  # @UndefinedVariable
        session = Session(bind=self.engine)
        self.session = session  # session()
        self.hit_i_string = str(Hit.__table__.insert(bind=self.engine).compile())
        self.hsp_i_string = str(Hsp.__table__.insert(bind=self.engine).compile())
        # Remove indices
        self.logger.debug("Created the session")
        # Load sequences if necessary
        self.__determine_sequences(query_seqs, target_seqs)
        self.xml = xml_name
        # Just a mock definition
        # self.get_query = functools.partial(self.__get_query_for_blast)
        self.not_pickable = ["manager", "printer_process",
                             "context", "logger_queue_handler", "queue_logger",
                             "log_writer",
                             "log_process", "pool", "main_logger",
                             "log_handler", "log_writer", "logger", "session",
                             "get_query", "engine", "query_seqs", "target_seqs"]

        self.queries, self.targets = dict(), dict()

        self.logger.debug("Finished __init__")

    def __getstate__(self):

        state = self.__dict__.copy()
        for not_pickable in self.not_pickable:
            if not_pickable in state:
                del state[not_pickable]

        return state

    def __determine_sequences(self, query_seqs, target_seqs):

        """Private method to assign the sequence file variables
        if necessary.
        :param query_seqs:
        :param target_seqs:
        :return:
        """

        if isinstance(query_seqs, str):
            assert os.path.exists(query_seqs)
            self.query_seqs = pysam.FastaFile(query_seqs)
        elif query_seqs is None:
            self.query_seqs = None
        else:
            self.logger.warn("Query type: %s", type(query_seqs))
            # assert "SeqIO.index" in repr(query_seqs)
            self.query_seqs = query_seqs

        self.target_seqs = []
        for target in target_seqs:
            if not os.path.exists(target):
                raise ValueError("{} not found!".format(target))
            self.target_seqs.append(pysam.FastaFile(target))

        return

    def __serialize_queries(self, queries):

        """Private method used to serialise the queries.
        Additional argument: a set containing all queries already present
        in the database."""

        counter = 0
        self.logger.info("Started to serialise the queries")
        objects = []
        for record, length in zip(self.query_seqs.references, self.query_seqs.lengths):
            if not record:
                continue
            if record in queries and queries[record][1] == length:
                continue
            elif record in queries:
                if queries[record][1] == length:
                    continue
                else:
                    raise KeyError("Discrepant length for query {} (original {}, new {})".format(
                        record, queries[record][1], length))
                # self.session.query(Query).filter(Query.query_name == record).update(
                #     {"query_length": length})
                # queries[record] = (queries[record][0], length)

            objects.append({
                "query_name": record,
                "query_length": length
            })

            if len(objects) >= self.maxobjects:
                self.logger.debug("Loading %d objects into the \"query\" table (total %d)",
                                 self.maxobjects, counter)

                # pylint: disable=no-member
                self.session.begin(subtransactions=True)
                # self.session.bulk_insert_mappings(Query, objects)
                self.engine.execute(Query.__table__.insert(), objects)
                # pylint: enable=no-member

                self.session.commit()
                counter += len(objects)
                objects = []

        if len(objects) > 0:
            self.logger.debug("Loading %d objects into the \"query\" table (total %d)",
                             len(objects), counter+len(objects))
            # pylint: disable=no-member
            counter += len(objects)
            # pylint: disable=no-member
            self.session.begin(subtransactions=True)
            self.engine.execute(Query.__table__.insert(), objects)
            # pylint: enable=no-member
            self.session.commit()
            # pylint: enable=no-member

        self.logger.info("Loaded %d objects into the \"query\" table", counter)

        if self.single_thread is True or self.procs == 1:
            _queries = pd.read_sql_table("query", self.engine, index_col="query_name",
                                        columns=["query_id", "query_length"])
            queries = dict((qname, int(qid)) for
                            qname, qid in zip(_queries.index, _queries["query_id"].values))
            self.logger.info("%d in queries", len(queries))
        else:
            queries = dict()

        return queries

    def __serialize_targets(self, targets):
        """
        This private method serialises all targets contained inside the target_seqs
        file into the database.
        :param targets: a cache dictionary which records whether the sequence
          is already present in the DB or not.
        """

        counter = 0
        objects = []
        self.logger.info("Started to serialise the targets")
        for target in self.target_seqs:
            for record, length in zip(target.references, target.lengths):
                if record in targets and targets[record][1] is True:
                    continue
                elif record in targets:
                    self.session.query(Target).filter(Target.target_name == record).update(
                        {"target_length": length})
                    targets[record] = (targets[record][0], True)
                    continue

                objects.append({
                    "target_name": record,
                    "target_length": length
                })
                counter += 1
                #
                # objects.append(Target(record, len(self.target_seqs[record])))
                if len(objects) >= self.maxobjects:
                    # counter += len(objects)
                    self.logger.debug("Loading %d objects into the \"target\" table",
                                     counter)
                    # self.session.bulk_insert_mappings(Target, objects)
                    self.session.begin(subtransactions=True)
                    self.engine.execute(Target.__table__.insert(), objects)
                    self.session.commit()
                    objects = []
        self.logger.debug("Loading %d objects into the \"target\" table, (total %d)",
                         len(objects), counter)
        # pylint: disable=no-member
        self.session.begin(subtransactions=True)
        self.engine.execute(Target.__table__.insert(), objects)
        # pylint: enable=no-member
        self.session.commit()

        self.logger.info("Loaded %d objects into the \"target\" table", counter)
        if self.single_thread is True or self.procs == 1:
            _targets = pd.read_sql_table("target", self.engine, index_col="target_name",
                                         columns=["target_id", "target_length"])
            targets = dict((qname, int(qid)) for
                            qname, qid in zip(_targets.index, _targets["target_id"].values))
            self.logger.debug("%d in targets", len(targets))
        else:
            targets = dict()
        return targets

    def __serialise_sequences(self):

        """ Private method called at the beginning of serialize. It is tasked
        with loading all necessary FASTA sequences into the DB and precaching the IDs.
        """

        targets = dict()
        queries = dict()
        self.logger.debug("Loading previous IDs")
        for query in self.session.query(Query):
            queries[query.query_name] = (query.query_id, query.query_length)
        for target in self.session.query(Target):
            targets[target.target_name] = (target.target_id, target.target_length)
        self.logger.debug("Loaded previous IDs; %d for queries, %d for targets",
                         len(queries), len(targets))

        self.logger.debug("Started the sequence serialisation")
        if self.target_seqs:
            targets = self.__serialize_targets(targets)
            # assert len(targets) > 0
        if self.query_seqs is not None:
            queries = self.__serialize_queries(queries)
            # assert len(queries) > 0

        return queries, targets

    def _load_into_db(self, hits, hsps, force=False):
        """

        :param hits:
        :param hsps:

        :param force: boolean flag. If set, data will be loaded no matter what.
        To be used at the end of the serialisation to load the final batch of data.
        :type force: bool

        :return:
        """
        load_into_db(self, hits, hsps, force=force)

    def serialize(self):

        """Method to serialize the BLAST XML file into a database
        provided with the __init__ method """
        self.queries, self.targets = self.__serialise_sequences()

        if isinstance(self.xml, str):
            if ".xml" in self.xml:
                self.__serialise_xmls()
            else:
                self.__serialise_tabular()
        else:
            assert isinstance(self.xml, (list, set))
            if all([".xml" in fname for fname in self.xml]):
                self.__serialise_xmls()
            else:
                self.__serialise_tabular()

    def __serialise_tabular(self):
        """Parser to perform the analysis of tabular BLAST files."""
        _serialise_tabular(self)

    def __serialise_xmls(self):
        _serialise_xmls(self)

    def __call__(self):
        """
        Alias for serialize
        """
        [idx.drop(bind=self.engine) for idx in Hit.__table__.indexes]
        [idx.drop(bind=self.engine) for idx in Hsp.__table__.indexes]
        self.engine.execute("PRAGMA foreign_keys=OFF")
        self.serialize()
        # Recreate the indices
        [idx.create(bind=self.engine) for idx in Hit.__table__.indexes]
        [idx.create(bind=self.engine) for idx in Hsp.__table__.indexes]
        self.engine.execute("PRAGMA foreign_keys=ON")

# pylint: enable=too-many-instance-attributes

    @staticmethod
    def get_multipliers(record, application=None):
        return get_multipliers(record, application=application)
