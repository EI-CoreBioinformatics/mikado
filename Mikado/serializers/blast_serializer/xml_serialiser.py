"""
XML serialisation class.
"""

import os
import functools
import logging.handlers as logging_handlers
import logging
import tempfile
try:
    import ujson as json
except ImportError:
    import json
import sqlite3
# import pickle
import sqlalchemy
import sqlalchemy.exc
from sqlalchemy.orm.session import Session
from ...utilities.dbutils import DBBASE
import pyfaidx
from ...utilities.dbutils import connect
from ...parsers.blast_utils import BlastOpener  # , XMLMerger
from ...utilities.log_utils import create_null_logger, check_logger
from . import Query, Target, Hsp, Hit, prepare_hit, InvalidHit
from xml.parsers.expat import ExpatError
import xml
from queue import Empty
import multiprocessing


__author__ = 'Luca Venturini'


# A serialisation class must have a ton of attributes ...
# pylint: disable=too-many-instance-attributes

class _XmlPickler(multiprocessing.Process):

    def __init__(self,
                 queries,
                 targets,
                 filequeue: multiprocessing.Queue,
                 returnqueue,
                 default_header,
                 identifier,
                 logging_queue,
                 level="WARN",
                 max_target_seqs=10,
                 maxobjects=20000):

        super().__init__()
        self.queries = queries
        self.targets = targets
        self.level = level
        self.logging_queue = logging_queue
        self.handler = logging_handlers.QueueHandler(logging_queue)
        self.handler.setLevel(self.level)
        self.__identifier = identifier
        self.name = self._name = "_XmlPickler-{0}".format(self.identifier)
        self.logger = logging.getLogger(self.name)
        self.logger.addHandler(self.handler)
        self.logger.setLevel(self.level)
        self.filequeue = filequeue
        self.returnqueue = returnqueue
        self.default_header = default_header
        self.maxobjects = maxobjects
        self.__max_target_seqs = max_target_seqs
        self.logger.debug("Started %s", self.name)

    def __getstate__(self):

        state = self.__dict__.copy()

        state["logger"].removeHandler(state["handler"])
        state["handler"].close()
        state["handler"] = None

        # state["_pickler"] = None
        state["logger"] = None
        return state

    def __setstate__(self, state):

        self.__dict__.update(state)
        self.handler = logging_handlers.QueueHandler(self.logging_queue)
        self.handler.setLevel(self.level)
        self.logger = logging.getLogger(self.name)
        self.logger.addHandler(self.handler)
        self.logger.setLevel(self.level)

    def _create_db(self, filename):

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
            self.logger.error(
                "Temporary db %s already exists (maybe from a previous aborted run?), dropping its contents",
                dbname)
            cursor.close()
            conn.close()
            os.remove(dbname)
            conn = sqlite3.connect(dbname)
            cursor = conn.cursor()
            cursor.execute(creation_string)
        cursor.execute("CREATE INDEX idx ON dump (query_counter)")
        self.logger.debug("Created tables for shelf %s", dbname)

        return dbname, conn, cursor

    def _pickler(self, filename):

        # Check the header is alright
        valid, _, exc = BlastOpener(filename).sniff(default_header=self.default_header)
        if not valid:
            self.logger.warning("Invalid BLAST file: %s", filename)
            return []

        self.logger.debug("Starting to pickle %s", filename)
        # hits, hsps = [], []

        pickle_count = 0
        query_counter = 0
        # Create the database
        dbname, conn, cursor = self._create_db(filename)

        try:
            with BlastOpener(filename) as opened:
                try:
                    for query_counter, record in enumerate(opened, start=1):
                        hits, hsps = objectify_record(self, record, [], [],
                                                      max_target_seqs=self.__max_target_seqs)

                        cursor.execute("INSERT INTO dump VALUES (?, ?, ?)",
                                       (query_counter, json.dumps(hits), json.dumps(hsps))
                                       )
                        if query_counter % self.maxobjects and query_counter > 0:
                            conn.commit()
                            cursor.close()
                            conn.close()
                            yield dbname
                            dbname, conn, cursor = self._create_db(filename)
                            pickle_count += 1

                except ExpatError:
                    self.logger.error("%s is an invalid BLAST file, sending back anything salvageable",
                                      filename)
                    raise
        except xml.etree.ElementTree.ParseError:
            self.logger.error("%s is an invalid BLAST file, sending back anything salvageable",
                              filename)
            raise
        except ValueError:
            self.logger.error("Invalid BLAST entry")
            raise

        self.logger.debug("Finished serialising %s in %s subsection", filename, pickle_count)
        cursor.close()
        conn.commit()
        conn.close()
        yield dbname
        # del records

    def run(self):
        """
        While running, the process will get the filenames to analyse from the first queue
        and return them through the second one.
        """

        while True:
            try:
                number, filename = self.filequeue.get(timeout=10)
            except multiprocessing.TimeoutError:
                self.logger.error(
                    "Something has gone awry in %s, no data received from the queue after waiting 10s. Aborting.",
                    self._name)
                # self.filequeue.put("EXIT")
                # return 0
                raise
            except Empty:
                continue

            if filename == "EXIT":
                self.logger.debug("Process %s received EXIT signal, terminating",
                                 self._name)
                self.filequeue.put((number, filename))
                return 0
            try:
                for pickled in self._pickler(filename):
                    self.logger.debug("Sending serialised information in {}".format(pickled))
                    self.returnqueue.put((number, [pickled]))
            except Exception as exc:
                self.logger.error("Error encountered in %s, blocking the program.", self.identifier)
                self.returnqueue.put((number, "FINISHED"))
                self.logger.exception(exc)
                raise

            self.returnqueue.put((number, "FINISHED"))

    @property
    def identifier(self):
        return self.__identifier


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

        self.procs = json_conf["serialise"]["procs"]
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
        session = Session(bind=self.engine, autocommit=True, autoflush=True)
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
            self.query_seqs = pyfaidx.Fasta(query_seqs)
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
            self.target_seqs.append(pyfaidx.Fasta(target))

        return

    def __serialize_queries(self, queries):

        """Private method used to serialise the queries.
        Additional argument: a set containing all queries already present
        in the database."""

        counter = 0
        self.logger.info("Started to serialise the queries")
        objects = []
        for record in self.query_seqs.records:
            if not record:
                continue
            if record in queries and queries[record][1] is not None:
                continue
            elif record in queries:
                self.session.query(Query).filter(Query.query_name == record).update(
                    {"query_length": len(self.query_seqs[record])})
                queries[record] = (queries[record][0], len(self.query_seqs[record]))
                continue

            objects.append({
                "query_name": record,
                "query_length": len(self.query_seqs[record])
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
        for query in self.session.query(Query):
            queries[query.query_name] = (query.query_id, query.query_length)
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
            for record in target.records:
                if record in targets and targets[record][1] is True:
                    continue
                elif record in targets:
                    self.session.query(Target).filter(Target.target_name == record).update(
                        {"target_length": len(self.target_seqs[record])})
                    targets[record] = (targets[record][0], True)
                    continue

                objects.append({
                    "target_name": record,
                    "target_length": len(target[record])
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
                    # pylint: disable=no-member
                    # pylint: enable=no-member
                    # self.logger.info("Loaded %d objects into the \"target\" table",
                    #                  len(objects))
                    # objects = []
        self.logger.debug("Loading %d objects into the \"target\" table, (total %d)",
                         len(objects), counter)
        # pylint: disable=no-member
        self.session.begin(subtransactions=True)
        self.engine.execute(Target.__table__.insert(), objects)
        # pylint: enable=no-member
        self.session.commit()

        self.logger.info("Loaded %d objects into the \"target\" table", counter)
        for target in self.session.query(Target):
            targets[target.target_name] = (target.target_id, target.target_length is not None)
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

        # Create the function that will retrieve the query_id given the name
        # self.get_query = functools.partial(self.__get_query_for_blast,
        #                                    **{"queries": self.queries})

        hits, hsps = [], []
        hit_counter, record_counter = 0, 0

        for filename in self.xml:
            valid, header, exc = BlastOpener(filename).sniff()
            if valid is True:
                self.header = header
            else:
                self.logger.error(exc)
                self.xml.remove(filename)

        if self.single_thread is True:
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
                            hits, hsps = objectify_record(self, record, hits, hsps,
                                                          max_target_seqs=self.__max_target_seqs)

                            hits, hsps = load_into_db(self, hits, hsps, force=False)
                    self.logger.debug("Finished %s", filename)
                except ExpatError:
                    self.logger.error("%s is an invalid BLAST file, saving what's available",
                                      filename)

            _, _ = self.__load_into_db(hits, hsps, force=True)

        else:
            self.logger.debug("Creating a pool with %d processes",
                             min(self.procs, len(self.xml)))

            filequeue = multiprocessing.Queue(-1)
            returnqueue = multiprocessing.Queue(-1)

            procs = [_XmlPickler(
                self.queries,
                self.targets,
                filequeue,
                returnqueue,
                self.header,
                _,
                logging_queue=self.logging_queue,
                # level=self.logger.level,
                level=self.json_conf["log_settings"]["log_level"],
                maxobjects=int(self.maxobjects/self.procs),
                max_target_seqs=self.__max_target_seqs
            )
                for _ in range(min([self.procs, len(self.xml)]))
                ]

            self.logger.debug("Starting to pickle and serialise %d files", len(self.xml))
            [_.start() for _ in procs]  # Start processes
            for number, xml_name in enumerate(self.xml):
                filequeue.put((number, xml_name))

            self.logger.debug("Finished sending off the data for serialisation")

            filequeue.put((None, "EXIT"))
            returned = []
            while len(returned) != len(self.xml):
                number, result = returnqueue.get()
                if result == "FINISHED":
                    self.logger.debug("Finished receiving pickles for %d", number)
                    returned.append(number)
                    continue
                if result == "EXIT":
                    continue
                for dbfile in result:
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
            [_.join() for _ in procs]  # Wait for processes to join
            self.logger.info("All %d children finished", len(procs))
            del procs

            _, _ = self.__load_into_db(hits, hsps, force=True)
            returnqueue.close()
            filequeue.close()

        self.logger.info("Loaded %d alignments for %d queries",
                         hit_counter, record_counter)

        self.logger.info("Finished loading blast hits")
        # [_.close() for _ in self.logger.handlers]
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
        :return:
        """

        q_mult, h_mult = 1, 1

        # application = record.application.upper()
        application = record.program.upper()

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


def _get_query_for_blast(self, record):
    """ This private method formats the name of the query
    recovered from the BLAST hit. It will cause an exception if the target is not
    present in the dictionary.
    :param record:
    :return: current_query (ID in the database), name
    """

    if record.id in self.queries:
        name = record.id
    else:
        name = record.id.split()[0]
        if name not in self.queries:
            raise KeyError("{} not found in the queries!".format(record))

    self.logger.debug("Started with %s", name)

    if self.queries[name][1] is False:
        raise KeyError("{} not found in the queries!".format(record))

    current_query = self.queries[name][0]
    return current_query, name


def _get_target_for_blast(self, alignment):
    """ This private method retrieves the correct target_id
    key for the target of the BLAST. If the entry is not present
    in the database, it will be created on the fly.
    The method returns the index of the current target and
    and an updated target dictionary.
    :param alignment: an alignment child of a BLAST record object
    :return: current_target (ID in the database), targets
    """

    if alignment.accession in self.targets:
        accession = alignment.accession
    elif alignment.id in self.targets:
        accession = alignment.id
    else:
        raise KeyError("{} not found in the targets!".format(alignment.accession))

    current_target = self.targets[accession][0]
    return current_target


def objectify_record(self, record, hits, hsps, max_target_seqs=10000):
    """
    Private method to serialise a single record into the DB.

    :param record: The BLAST record to load into the DB.
    :param hits: Cache of hits to load into the DB.
    :type hits: list

    :param hsps: Cache of hsps to load into the DB.
    :type hsps: list

    :returns: hits, hsps
    :rtype: (list, list)
    """

    if len(record.hits) == 0:
        return hits, hsps

    current_query, name = _get_query_for_blast(self, record)

    current_evalue = -1
    current_counter = 0

    # for ccc, alignment in enumerate(record.alignments):
    for ccc, alignment in enumerate(record.hits):
        if ccc + 1 > max_target_seqs:
            break

        self.logger.debug("Started the hit %s vs. %s",
                          # name, record.alignments[ccc].accession)
                          name, record.hits[ccc].id)
        current_target = _get_target_for_blast(self, alignment)

        hit_dict_params = dict()
        (hit_dict_params["query_multiplier"],
         hit_dict_params["target_multiplier"]) = XmlSerializer.get_multipliers(record)
        hit_evalue = min(_.evalue for _ in record.hits[ccc].hsps)
        hit_bs = max(_.bitscore for _ in record.hits[ccc].hsps)
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
            self.logger.error(exc)
            continue
        hits.append(hit)
        hsps.extend(hit_hsps)

    return hits, hsps
