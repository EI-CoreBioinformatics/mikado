"""
XML serialisation class.
"""

import os
import logging.handlers as logging_handlers
import logging
import tempfile
import ujson as json
import sqlite3
import sqlalchemy
import sqlalchemy.exc
from sqlalchemy.orm.session import Session
from ...utilities.dbutils import DBBASE
import pysam
from ...utilities.dbutils import connect
from ...parsers.blast_utils import BlastOpener  # , XMLMerger
from ...utilities.log_utils import create_null_logger, check_logger
from . import Query, Target, Hsp, Hit, prepare_hit, InvalidHit
from xml.parsers.expat import ExpatError
import xml
import pandas as pd
# from queue import Empty
import multiprocessing


__author__ = 'Luca Venturini'


# A serialisation class must have a ton of attributes ...
# pylint: disable=too-many-instance-attributes


def _create_xml_db(filename):
    """Private method to create a DB for serialisation.
    :param filename: the name of the file to serialise
    :returns dbname, cursor: the name of the database and the SQLite cursor

    """

    directory = os.path.dirname(filename)
    try:
        dbname = tempfile.mktemp(suffix=".db", dir=directory)
        conn = sqlite3.connect(dbname)
    except (OSError, PermissionError, sqlite3.OperationalError):
        dbname = tempfile.mktemp(suffix=".db")
        conn = sqlite3.connect(dbname)
    cursor = conn.cursor()
    creation_string = "CREATE TABLE dump (query_counter integer, hits blob, hsps blob)"
    try:
        cursor.execute(  # TODO: change
            creation_string)
    except sqlite3.OperationalError:
        # Table already exists
        cursor.close()
        conn.close()
        os.remove(dbname)
        conn = sqlite3.connect(dbname)
        cursor = conn.cursor()
        cursor.execute(creation_string)
    cursor.execute("CREATE INDEX idx ON dump (query_counter)")
    return dbname, conn, cursor


def xml_pickler(json_conf, filename, default_header,
                max_target_seqs=10):
    valid, _, exc = BlastOpener(filename).sniff(default_header=default_header)
    engine = connect(json_conf, strategy="threadlocal")
    session = Session(bind=engine)

    if not valid:
        err = "Invalid BLAST file: %s" % filename
        raise TypeError(err)
    dbname, conn, cursor = _create_xml_db(filename)
    try:
        with BlastOpener(filename) as opened:
            try:
                for query_counter, record in enumerate(opened, start=1):
                    hits, hsps = objectify_record(session, record, [], [],
                                                  max_target_seqs=max_target_seqs)

                    cursor.execute("INSERT INTO dump VALUES (?, ?, ?)",
                                   (query_counter, json.dumps(hits), json.dumps(hsps))
                                   )
            except ExpatError as err:
                # logger.error("%s is an invalid BLAST file, sending back anything salvageable", filename)
                raise ExpatError("{} is an invalid BLAST file, sending back anything salvageable.\n{}".format(filename,
                                                                                                              err))
    except xml.etree.ElementTree.ParseError as err:
        # logger.error("%s is an invalid BLAST file, sending back anything salvageable", filename)
        raise xml.etree.ElementTree.ParseError(
            "{} is an invalid BLAST file, sending back anything salvageable.\n{}".format(filename, err))
    except ValueError as err:
        # logger.error("Invalid BLAST entry")
        raise ValueError(
            "{} is an invalid BLAST file, sending back anything salvageable.\n{}".format(filename, err))

    cursor.close()
    conn.commit()
    conn.close()
    return dbname


class XmlSerializer:
    """This class has the role of taking in input a blast XML file and (partially)
    serialise it into a database. We are using SQLalchemy, so the database type
    could be any of SQLite, MySQL, PSQL, etc."""

    __name__ = "XMLSerializer"
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
        self.logging_queue = multiprocessing.Queue(-1)
        self.logger_queue_handler = logging_handlers.QueueHandler(self.logging_queue)
        self.queue_logger = logging.getLogger("parser")
        self.queue_logger.addHandler(self.logger_queue_handler)
        self.queue_logger.setLevel(self.json_conf["log_settings"]["log_level"])
        self.queue_logger.propagate = False
        self.log_writer = logging_handlers.QueueListener(self.logging_queue, self.logger)
        self.log_writer.start()

        self.__max_target_seqs = json_conf["serialise"]["max_target_seqs"]
        self.maxobjects = json_conf["serialise"]["max_objects"]
        target_seqs = json_conf["serialise"]["files"]["blast_targets"]
        query_seqs = json_conf["serialise"]["files"]["transcripts"]

        self.header = None
        if xml_name is None:
            self.logger.warning("No BLAST XML provided. Exiting.")
            return

        self.engine = connect(json_conf)

        # session = sessionmaker(autocommit=True)
        DBBASE.metadata.create_all(self.engine)  # @UndefinedVariable
        session = Session(bind=self.engine, autocommit=False, autoflush=False, expire_on_commit=False)
        self.session = session  # session()
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
            if record in queries and queries[record][1] is not None:
                continue
            elif record in queries:
                self.session.query(Query).filter(Query.query_name == record).update(
                    {"query_length": length})
                queries[record] = (queries[record][0], length)
                continue

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
            self.session.begin()
            self.engine.execute(Query.__table__.insert(), objects)
            # pylint: enable=no-member
            self.session.commit()
            # pylint: enable=no-member
            self.logger.debug("Loaded %d objects into the \"query\" table", counter)

        _queries = pd.read_sql_table("query", self.engine, index_col="query_name",
                                    columns=["query_id", "query_length"])
        _queries_d = _queries.to_dict("list")
        queries = dict((key, (qid, qlen)) for (key, qid, qlen) in
                          zip(_queries.index, _queries_d["query_id"], _queries_d["query_length"]))
        self.logger.info("%d in queries", len(queries))
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
        _targets = pd.read_sql_table("target", self.engine, index_col="target_name",
                                     columns=["target_id", "target_length"])
        _targets_d = _targets.to_dict("list")
        targets = dict((key, (qid, qlen)) for (key, qid, qlen) in
                       zip(_targets.index, _targets_d["target_id"], _targets_d["target_length"]))
        self.logger.debug("%d in targets", len(targets))
        return targets

    def __serialise_sequences(self):

        """ Private method called at the beginning of serialize. It is tasked
        with loading all necessary FASTA sequences into the DB and precaching the IDs.
        """

        targets = dict()
        queries = dict()
        self.logger.debug("Loading previous IDs")
        for query in self.session.query(Query):
            queries[query.query_name] = (query.query_id, (query.query_length is not None))
        for target in self.session.query(Target):
            targets[target.target_name] = (target.target_id, (target.target_length is not None))
        self.logger.debug("Loaded previous IDs; %d for queries, %d for targets",
                         len(queries), len(targets))

        self.logger.debug("Started the sequence serialisation")
        if self.target_seqs:
            targets = self.__serialize_targets(targets)
            assert len(targets) > 0
        if self.query_seqs is not None:
            queries = self.__serialize_queries(queries)
            assert len(queries) > 0

        return queries, targets

    def __load_into_db(self, hits, hsps, force=False):
        """

        :param hits:
        :param hsps:

        :param force: boolean flag. If set, data will be loaded no matter what.
        To be used at the end of the serialisation to load the final batch of data.
        :type force: bool

        :return:
        """

        self.logger.debug("Checking whether to load %d hits and %d hsps",
                          len(hits), len(hsps))

        tot_objects = len(hits) + len(hsps)
        if len(hits) == 0:
            self.logger.debug("No hits to serialise. Exiting")
            return hits, hsps

        if tot_objects >= self.maxobjects or force:
            # Bulk load
            self.logger.debug("Loading %d BLAST objects into database", tot_objects)

            try:
                # pylint: disable=no-member
                self.session.begin(subtransactions=True)
                self.engine.execute(Hit.__table__.insert(), hits)
                self.engine.execute(Hsp.__table__.insert(), hsps)
                # pylint: enable=no-member
                self.session.commit()
            except sqlalchemy.exc.IntegrityError as err:
                self.logger.critical("Failed to serialise BLAST!")
                self.logger.exception(err)
                raise err
            self.logger.debug("Loaded %d BLAST objects into database", tot_objects)
            hits, hsps = [], []
        return hits, hsps

    def serialize(self):

        """Method to serialize the BLAST XML file into a database
        provided with the __init__ method """

        # Load sequences in DB, precache IDs
        self.queries, self.targets = self.__serialise_sequences()
        if isinstance(self.xml, str):
            self.xml = [self.xml]
        else:
            assert isinstance(self.xml, (list, set))

        hits, hsps = [], []
        hit_counter, record_counter = 0, 0

        for filename in self.xml:
            valid, header, exc = BlastOpener(filename).sniff()
            if valid is True:
                self.header = header
            else:
                self.logger.error(exc)
                self.xml.remove(filename)

        if self.single_thread is True or self.procs == 1:
            for filename in self.xml:
                valid, _, exc = BlastOpener(filename).sniff(default_header=self.header)
                if not valid:
                    self.logger.error(exc)
                    continue
                try:
                    self.logger.debug("Analysing %s", filename)
                    with BlastOpener(filename) as opened:
                        for record in opened:
                            record_counter += 1
                            if record_counter > 0 and record_counter % 10000 == 0:
                                self.logger.info("Parsed %d queries", record_counter)
                            hits, hsps = objectify_record(self.session,
                                                          record, hits, hsps,
                                                          max_target_seqs=self.__max_target_seqs, logger=self.logger)

                            hits, hsps = load_into_db(self, hits, hsps, force=False)
                    self.logger.debug("Finished %s", filename)
                except ExpatError:
                    self.logger.error("%s is an invalid BLAST file, saving what's available", filename)
            _, _ = self.__load_into_db(hits, hsps, force=True)

        else:
            self.logger.debug("Creating a pool with %d processes",
                              min(self.procs, len(self.xml)))

            pool = multiprocessing.Pool(self.procs)
            results = []
            for num, filename in enumerate(self.xml):
                args = (self.json_conf, filename, self.header)
                kwds = {"max_target_seqs": self.__max_target_seqs}
                pool.apply_async(xml_pickler, args=args, kwds=kwds, callback=results.append)
            pool.close()
            pool.join()

            for dbfile in results:
                conn = sqlite3.connect(dbfile)
                cursor = conn.cursor()
                for query_counter, __hits, __hsps in cursor.execute("SELECT * FROM dump"):
                    record_counter += 1
                    __hits = json.loads(__hits)
                    __hsps = json.loads(__hsps)
                    hit_counter += len(__hits)
                    hits.extend(__hits)
                    hsps.extend(__hsps)
                    hits, hsps = load_into_db(self, hits, hsps, force=False)
                    if record_counter > 0 and record_counter % 10000 == 0:
                        self.logger.debug("Parsed %d queries", record_counter)
                cursor.close()
                conn.close()
                os.remove(dbfile)

            self.logger.debug("Finished sending off the data for serialisation")
            _, _ = self.__load_into_db(hits, hsps, force=True)

        self.logger.info("Loaded %d alignments for %d queries",
                         hit_counter, record_counter)

        self.logger.info("Finished loading blast hits")
        if hasattr(self, "logging_queue"):
            self.logging_queue.close()

    def __call__(self):
        """
        Alias for serialize
        """
        self.serialize()

    @staticmethod
    def get_multipliers(record):

        """
        Private quick method to determine the multipliers for a BLAST alignment
        according to the application present in the record.
        :param record:
        :type record: Bio.Blast.Record.Blast
        :return:
        """

        q_mult, h_mult = 1, 1

        # application = record.application.upper()
        application = record.application.upper()

        if application in ("BLASTN", "TBLASTX", "BLASTP"):
            q_mult = 1
            h_mult = 1
        elif application == "BLASTX":
            q_mult = 3
            h_mult = 1
        elif application == "TBLASTN":
            q_mult = 1
            h_mult = 3

        return q_mult, h_mult

# pylint: enable=too-many-instance-attributes


def load_into_db(self, hits, hsps, force=False):
    """
    :param hits:
    :param hsps:

    :param force: boolean flag. If set, data will be loaded no matter what.
    To be used at the end of the serialisation to load the final batch of data.
    :type force: bool

    :return:
    """

    self.logger.debug("Checking whether to load %d hits and %d hsps",
                      len(hits), len(hsps))

    tot_objects = len(hits) + len(hsps)
    if len(hits) == 0:
        self.logger.debug("No hits to serialise. Exiting")
        return hits, hsps

    if tot_objects >= self.maxobjects or force:
        # Bulk load
        self.logger.debug("Loading %d BLAST objects into database", tot_objects)

        try:
            # pylint: disable=no-member
            self.session.begin(subtransactions=True)
            self.engine.execute(Hit.__table__.insert(), hits)
            self.engine.execute(Hsp.__table__.insert(), hsps)
            # pylint: enable=no-member
            self.session.commit()
        except sqlalchemy.exc.IntegrityError as err:
            self.logger.critical("Failed to serialise BLAST!")
            self.logger.exception(err)
            raise err
        self.logger.debug("Loaded %d BLAST objects into database", tot_objects)
        hits, hsps = [], []
    return hits, hsps


def _get_query_for_blast(session: sqlalchemy.orm.session.Session, record):
    """ This private method formats the name of the query
    recovered from the BLAST hit. It will cause an exception if the target is not
    present in the dictionary.
    :param record:
    :return: current_query (ID in the database), name
    """

    got = session.query(Query).filter(sqlalchemy.or_(
        Query.query_name == record.query,
        Query.query_name == record.query.split()[0],
    )).one()
    return got.query_id, got.query_name


def _get_target_for_blast(session, alignment):
    """ This private method retrieves the correct target_id
    key for the target of the BLAST. If the entry is not present
    in the database, it will be created on the fly.
    The method returns the index of the current target and
    and an updated target dictionary.
    :param alignment: an alignment child of a BLAST record object
    :return: current_target (ID in the database), targets
    """

    got = session.query(Target).filter(sqlalchemy.or_(
        Target.target_name == alignment.accession,
        Target.target_name == alignment.hit_id)).one()

    # current_target = targets[accession][0]
    return got.target_id


def objectify_record(session, record, hits, hsps, max_target_seqs=10000, logger=create_null_logger()):
    """
    Private method to serialise a single record into the DB.

    :param record: The BLAST record to load into the DB.
    :type record: Bio.Blast.Record.Blast
    :param hits: Cache of hits to load into the DB.
    :type hits: list

    :param hsps: Cache of hsps to load into the DB.
    :type hsps: list

    :returns: hits, hsps
    :rtype: (list, list)
    """

    if len(record.alignments) == 0:
        return hits, hsps

    current_query, name = _get_query_for_blast(session, record)

    current_evalue = -1
    current_counter = 0

    # for ccc, alignment in enumerate(record.alignments):
    for ccc, alignment in enumerate(record.alignments):
        if ccc + 1 > max_target_seqs:
            break

        logger.debug("Started the hit %s vs. %s", name, record.alignments[ccc].hit_id)
        current_target = _get_target_for_blast(session, alignment)

        hit_dict_params = dict()
        (hit_dict_params["query_multiplier"],
         hit_dict_params["target_multiplier"]) = XmlSerializer.get_multipliers(record)
        hit_evalue = min(_.expect for _ in record.alignments[ccc].hsps)
        hit_bs = max(_.score for _ in record.alignments[ccc].hsps)
        if current_evalue < hit_evalue:
            current_counter += 1
            current_evalue = hit_evalue

        hit_dict_params["hit_number"] = current_counter
        hit_dict_params["evalue"] = hit_evalue
        hit_dict_params["bits"] = hit_bs

        # Prepare for bulk load
        try:
            hit, hit_hsps = prepare_hit(alignment, current_query,
                                        current_target, **hit_dict_params)
        except InvalidHit as exc:
            logger.error(exc)
            continue
        hits.append(hit)
        hsps.extend(hit_hsps)

    return hits, hsps
