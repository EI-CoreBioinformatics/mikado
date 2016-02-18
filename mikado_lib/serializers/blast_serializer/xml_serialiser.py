"""
XML serialisation class.
"""

import os
import functools
import logging.handlers as logging_handlers
import logging
import tempfile
import pickle
import itertools
import sqlalchemy
import sqlalchemy.exc
from Bio.Blast.NCBIXML import parse as xparser
from sqlalchemy.orm.session import sessionmaker
from ...utilities.dbutils import DBBASE
import pyfaidx
from ...utilities.dbutils import connect
from ...parsers.blast_utils import create_opener  # , XMLMerger
from ...parsers import HeaderError
from ...utilities.log_utils import create_null_logger, check_logger
from . import Query, Target, Hsp, Hit, prepare_hit
from xml.parsers.expat import ExpatError
import multiprocessing

__author__ = 'Luca Venturini'


# A serialisation class must have a ton of attributes ...
# pylint: disable=too-many-instance-attributes

def _pickle_xml(filename, default_header, maxobjects, logging_queue):

    """
    Private method to load the records from an XML file into a pickled file,
    for faster loading.
    :param filename:
    :param logging_queue: the queue to be used for logging

    :return: a list of the
    """

    records = []
    pfiles = []

    handler = logging_handlers.QueueHandler(logging_queue)
    logger = logging.getLogger("{0}:{1}-{2}".format(filename))
    logger.addHandler(handler)
    logger.setLevel("DEBUG")

    # Check the header is alright
    assert default_header is not None
    valid, _, exc = _sniff(filename, default_header=default_header)
    if not valid:
        logger.warning("Invalid BLAST file: %s", filename)
        return pfiles

    logger.info("Starting to pickle %s", filename)
    try:
        for record in xparser(create_opener(filename)):
            if len(record.descriptions) > 0:
                records.append(record)
            if len(records) > maxobjects:
                pickle_temp = tempfile.mkstemp(suffix=".pickle",
                                               dir=os.path.dirname(filename))
                with open(pickle_temp[1], "wb") as pickled:
                    pickle.dump(records, pickled)
                pfiles.append(pickle_temp[1])
                records = []
    except ExpatError:
        logger.error("%s is an invalid BLAST file, sending back anything salvageable",
                     filename)

    pickle_temp = tempfile.mkstemp(suffix=".pickle")
    with open(pickle_temp[1], "wb") as pickled:
        pickle.dump(records, pickled)
    pfiles.append(pickle_temp[1])
    logger.info("Finished pickling %s", filename)
    del records
    return pfiles


def _sniff(filename, default_header=None):

    """
    Function that either derives the default XML header for the instance (if undefined)
    or checks that the given file is compatible with it.
    :param filename: The filename to check for consistency.
    :return: boolean (passed or not passed)
    :rtype: (bool, list, str)
    """

    handle = create_opener(filename)
    header = []
    exc = None
    valid = True
    while True:
        line = next(handle)

        if "<Iteration>" in line:
            break
        line = line.rstrip()
        if not line:
            exc = HeaderError("Invalid header for {0}:\n\n{1}".format(
                filename,
                "\n".join(header)
            ))
            break
        if len(header) > 10**3:
            exc = HeaderError("Abnormally long header ({0}) for {1}:\n\n{2}".format(
                len(header),
                filename,
                "\n".join(header)
            ))
            break
        header.append(line)
    if not any(iter(True if "BlastOutput" in x else False for x in header)):
        exc = HeaderError("Invalid header for {0}:\n\n{1}".format(
            filename, "\n".join(header)))

    if default_header is not None and exc is None:
            checker = [header_line for header_line in header if
                       "BlastOutput_query" not in header_line]
            previous_header = [header_line for header_line in default_header if
                               "BlastOutput_query" not in header_line]
            if checker != previous_header:
                exc = HeaderError("BLAST XML header does not match for {0}".format(
                    filename))
    elif exc is None:
        default_header = header
    handle.close()

    if exc is not None:
        valid = False

    return valid, default_header, exc


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

    def __init__(self, xml,
                 logger=None,
                 json_conf=None):
        """Initializing method. Arguments:

        :param xml: The XML to parse.

        Arguments:

        :param json_conf: a configuration dictionary.
        :type json_conf: dict


        """

        if json_conf is None:
            raise ValueError("No configuration provided!")

        if logger is not None:
            self.logger = check_logger(logger)

        # Runtime arguments
        self.threads = json_conf["serialise"]["threads"]
        self.single_thread = json_conf["serialise"]["single_thread"]
        if self.threads > 1 and self.single_thread is False:
            self.context = multiprocessing.get_context()
            self.manager = self.context.Manager()
            self.logging_queue = self.manager.Queue(-1)
            self.logger_queue_handler = logging_handlers.QueueHandler(self.logging_queue)
            self.queue_logger = logging.getLogger("parser")
            self.queue_logger.addHandler(self.logger_queue_handler)
            self.queue_logger.setLevel(json_conf["log_settings"]["log_level"])
            self.queue_logger.propagate = False
            self.log_writer = logging_handlers.QueueListener(self.logging_queue, self.logger)
            self.log_writer.start()

        self.discard_definition = json_conf["serialise"]["discard_definition"]
        self.__max_target_seqs = json_conf["serialise"]["max_target_seqs"]
        self.maxobjects = json_conf["serialise"]["max_objects"]
        target_seqs = json_conf["serialise"]["files"]["blast_targets"]
        query_seqs = json_conf["serialise"]["files"]["transcripts"]

        self.header = None
        if xml is None:
            self.logger.warning("No BLAST XML provided. Exiting.")
            return

        self.engine = connect(json_conf)

        session = sessionmaker()
        session.configure(bind=self.engine)
        DBBASE.metadata.create_all(self.engine)  # @UndefinedVariable
        self.session = session()
        self.logger.debug("Created the session")
        # Load sequences if necessary
        self.__determine_sequences(query_seqs, target_seqs)
        self.xml = xml
        # Just a mock definition
        self.get_query = functools.partial(self.__get_query_for_blast)
        self.not_pickable = ["manager", "printer_process",
                             "context", "logger_queue_handler", "queue_logger",
                             "log_writer",
                             "log_process", "pool", "main_logger",
                             "log_handler", "log_writer", "logger", "session",
                             "get_query", "engine", "query_seqs", "target_seqs"]

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

        if isinstance(target_seqs, str):
            assert os.path.exists(target_seqs)
            self.target_seqs = pyfaidx.Fasta(target_seqs)
        elif target_seqs is None:
            self.target_seqs = None
        else:
            self.logger.warn("Target (%s) type: %s",
                             target_seqs,
                             type(target_seqs))
            # assert "SeqIO.index" in repr(target_seqs)
            self.target_seqs = target_seqs
        return

    def __serialize_queries(self, queries):

        """Private method used to serialise the queries.
        Additional argument: a set containing all queries already present
        in the database."""

        counter = 0
        self.logger.info("Started to serialise the queries")
        objects = []
        for record in self.query_seqs.records:
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
            #
            # objects.append(Target(record, len(self.target_seqs[record])))
            if len(objects) >= self.maxobjects:
                self.logger.info("Loading %d objects into the \"query\" table (total %d)",
                                 self.maxobjects, counter)
                # self.session.bulk_insert_mappings(Query, objects)
                # pylint: disable=no-member
                self.engine.execute(Query.__table__.insert(), objects)
                # pylint: enable=no-member

                self.session.commit()
                counter += len(objects)
                objects = []
                # pylint: disable=no-member
                # # pylint: enable=no-member
                # self.logger.info("Loaded %d objects into the \"target\" table",
                #                  len(objects))
                # objects = []
        self.logger.info("Loading %d objects into the \"query\" table (total %d)",
                         len(objects), counter+len(objects))
        # pylint: disable=no-member
        # self.engine.execute(Target.__table__.insert(),
        #                     [{"target_name": obj.target_name,
        #                       "target_length": obj.target_length} for obj in objects])
        counter += len(objects)
        # pylint: disable=no-member
        self.engine.execute(Query.__table__.insert(), objects)
        # pylint: enable=no-member
        # self.session.bulk_insert_mappings(Query, objects)
        self.session.commit()
        # pylint: enable=no-member
        self.logger.info("Loaded %d objects into the \"query\" table", counter)
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
        for record in self.target_seqs.records:
            if record in targets and targets[record][1] is True:
                continue
            elif record in targets:
                self.session.query(Target).filter(Target.target_name == record).update(
                    {"target_length": len(self.target_seqs[record])})
                targets[record] = (targets[record][0], True)
                continue

            objects.append({
                "target_name": record,
                "target_length": len(self.target_seqs[record])
            })
            counter += 1
            #
            # objects.append(Target(record, len(self.target_seqs[record])))
            if len(objects) >= self.maxobjects:
                # counter += len(objects)
                self.logger.info("Loading %d objects into the \"target\" table",
                                 counter)
                # self.session.bulk_insert_mappings(Target, objects)
                self.engine.execute(Target.__table__.insert(), objects)
                self.session.commit()
                objects = []
                # pylint: disable=no-member
                # pylint: enable=no-member
                # self.logger.info("Loaded %d objects into the \"target\" table",
                #                  len(objects))
                # objects = []
        self.logger.info("Loading %d objects into the \"target\" table, (total %d)",
                         len(objects), counter)
        # pylint: disable=no-member
        self.engine.execute(Target.__table__.insert(), objects)
        # pylint: enable=no-member
        self.session.commit()

        self.logger.info("Loaded %d objects into the \"target\" table", counter)
        for target in self.session.query(Target):
            targets[target.target_name] = (target.target_id, target.target_length is not None)
        self.logger.info("%d in targets", len(targets))
        return targets

    def __serialise_sequences(self):

        """ Private method called at the beginning of serialize. It is tasked
        with loading all necessary FASTA sequences into the DB and precaching the IDs.
        """

        targets = dict()
        queries = dict()
        self.logger.info("Loading previous IDs")
        for query in self.session.query(Query):
            queries[query.query_name] = (query.query_id, (query.query_length is not None))
        for target in self.session.query(Target):
            targets[target.target_name] = (target.target_id, (target.target_length is not None))
        self.logger.info("Loaded previous IDs; %d for queries, %d for targets",
                         len(queries), len(targets))

        self.logger.info("Started the sequence serialisation")
        if self.target_seqs is not None:
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
        if tot_objects >= self.maxobjects or force:
            # Bulk load
            self.logger.info("Loading %d BLAST objects into database", tot_objects)

            try:
                # pylint: disable=no-member
                self.engine.execute(Hit.__table__.insert(), hits)
                self.engine.execute(Hsp.__table__.insert(), hsps)
                # pylint: enable=no-member
                self.session.commit()
            except sqlalchemy.exc.IntegrityError as err:
                self.logger.critical("Failed to serialise BLAST!")
                self.logger.exception(err)
                raise err
            self.logger.info("Loaded %d BLAST objects into database", tot_objects)
            hits, hsps = [], []
        return hits, hsps

    def __serialise_record(self, record, hits, hsps, targets):
        """
        Private method to serialise a single record into the DB.

        :param record: The BLAST record to load into the DB.
        :param hits: Cache of hits to load into the DB.
        :type hits: list

        :param hsps: Cache of hsps to load into the DB.
        :type hsps: list

        :param targets: dictionary which holds the relationship target ID/name.
         It can be updated in place if a target has not been serialised already.

        :returns: hits, hsps, hit_counter, targets
        :rtype: (list, list, int, dict)
        """

        hit_counter = 0

        if len(record.descriptions) == 0:
            return hits, hsps, hit_counter, targets

        current_query, name = self.get_query(record)

        current_evalue = -1
        current_counter = 0

        for ccc, alignment in enumerate(record.alignments):
            if ccc + 1 > self.__max_target_seqs:
                break

            self.logger.debug("Started the hit %s vs. %s",
                              name, record.alignments[ccc].accession)
            try:
                current_target, targets = self.__get_target_for_blast(alignment,
                                                                      targets)
            except sqlalchemy.exc.IntegrityError as exc:
                self.session.rollback()
                self.logger.exception(exc)
                continue

            hit_counter += 1

            hit_dict_params = dict()
            (hit_dict_params["query_multiplier"],
             hit_dict_params["target_multiplier"]) = self.__get_multipliers(record)
            if current_evalue < record.descriptions[ccc].e:
                current_counter += 1
                current_evalue = record.descriptions[ccc].e

            hit_dict_params["hit_number"] = current_counter
            hit_dict_params["evalue"] = record.descriptions[ccc].e
            hit_dict_params["bits"] = record.descriptions[ccc].bits

            # Prepare for bulk load
            hit, hit_hsps = prepare_hit(alignment, current_query, current_target,
                                        **hit_dict_params)
            hits.append(hit)
            hsps.extend(hit_hsps)

            hits, hsps = self.__load_into_db(hits, hsps, force=False)

        return hits, hsps, hit_counter, targets

    def serialize(self):

        """Method to serialize the BLAST XML file into a database
        provided with the __init__ method """

        # Load sequences in DB, precache IDs
        queries, targets = self.__serialise_sequences()
        if isinstance(self.xml, str):
            self.xml = [self.xml]
            # self.xml_parser = xparser(create_opener(self.xml))
        else:
            assert isinstance(self.xml, (list, set))
            # if len(self.xml) < 1:
            #     raise ValueError("No input file provided!")
            # elif len(self.xml) == 1:
            #     self.xml_parser = xparser(create_opener(list(self.xml)[0]))
            # else:
            #     self.xml_parser = xparser(XMLMerger(self.xml))  # Merge in memory

        # Create the function that will retrieve the query_id given the name
        self.get_query = functools.partial(self.__get_query_for_blast,
                                           **{"queries": queries})

        hits, hsps = [], []
        hit_counter, record_counter = 0, 0

        for filename in self.xml:
            valid, header, exc = _sniff(filename)
            if valid is True:
                self.header = header
            else:
                self.logger.error(exc)
                self.xml.remove(filename)

        if self.threads == 1 or self.single_thread is True:
            for filename in self.xml:
                valid, _, exc = _sniff(filename, default_header=self.header)
                if not valid:
                    self.logger.error(exc)
                    continue
                try:
                    for record in xparser(create_opener(filename)):
                        record_counter += 1
                        if record_counter > 0 and record_counter % 10000 == 0:
                            self.logger.info("Parsed %d queries", record_counter)

                        hits, hsps, partial_hit_counter, targets = self.__serialise_record(record,
                                                                                           hits,
                                                                                           hsps,
                                                                                           targets)
                        hit_counter += partial_hit_counter
                        if hit_counter > 0 and hit_counter % 10000 == 0:
                            self.logger.info("Serialized %d alignments", hit_counter)
                except ExpatError:
                    self.logger.error("%s is an invalid BLAST file, saving what's available",
                                      filename)

            _, _ = self.__load_into_db(hits, hsps, force=True)

        else:
            self.logger.info("Creating a pool with %d processes",
                             min(self.threads, len(self.xml)))

            pool = multiprocessing.Pool(min(self.threads, len(self.xml)))
            args = zip(
                self.xml,
                [self.header] * len(self.xml),
                [self.maxobjects] * len(self.xml),
                [self.logging_queue] * len(self.xml)
            )

            pickle_results = pool.starmap_async(_pickle_xml, args)
            self.logger.info("Starting to pickle and serialise %d files", len(self.xml))
            for pickle_file in itertools.chain(*pickle_results.get()):
                with open(pickle_file, "rb") as pickled:
                    for record in pickle.load(pickled):
                        record_counter += 1
                        if record_counter > 0 and record_counter % 10000 == 0:
                            self.logger.info("Parsed %d queries", record_counter)

                        hits, hsps, partial_hit_counter, targets = self.__serialise_record(record,
                                                                                           hits,
                                                                                           hsps,
                                                                                           targets)
                        hit_counter += partial_hit_counter
                        if hit_counter > 0 and hit_counter % 10000 == 0:
                            self.logger.info("Serialized %d alignments", hit_counter)
                os.remove(pickle_file)

            _, _ = self.__load_into_db(hits, hsps, force=True)

        self.logger.info("Loaded %d alignments for %d queries",
                         hit_counter, record_counter)

        self.logger.info("Finished loading blast hits")

    def __call__(self):
        """
        Alias for serialize
        """
        self.serialize()

    @staticmethod
    def __get_multipliers(record):

        """
        Private quick method to determine the multipliers for a BLAST alignment
        according to the application present in the record.
        :param record:
        :return:
        """

        q_mult, h_mult = 1, 1

        if record.application in ("BLASTN", "TBLASTX", "BLASTP"):
            q_mult = 1
            h_mult = 1
        elif record.application == "BLASTX":
            q_mult = 3
            h_mult = 1
        elif record.application == "TBLASTN":
            q_mult = 1
            h_mult = 3

        return q_mult, h_mult

    def __get_query_for_blast(self, record, queries):

        """ This private method formats the name of the query
        recovered from the BLAST hit, verifies whether it is present or not
        in the DB, and moreover whether the data can be updated (e.g.
        by adding the query length)
        :param record:
        :param queries:
        :return: current_query (ID in the database), name
        """

        if self.discard_definition is False:
            name = record.query.split()[0]
        else:
            name = record.query_id
        self.logger.debug("Started with %s", name)

        if name in queries:
            try:
                current_query = queries[name][0]
            except TypeError as exc:
                raise TypeError("{0}, {1}".format(exc, name))
            if queries[name][1] is False:
                self.session.query(Query).filter(Query.query_name == name).update(
                    {"query_length": record.query_length})
                self.session.commit()
        else:
            self.logger.warn("%s not found among queries, adding to the DB now",
                             name)
            current_query = Query(name, record.query_length)
            self.session.add(current_query)
            self.session.commit()
            queries[name] = (current_query.query_id, True)
            current_query = current_query.query_id
        return current_query, name

    def __get_target_for_blast(self, alignment, targets):

        """ This private method retrieves the correct target_id
        key for the target of the BLAST. If the entry is not present
        in the database, it will be created on the fly.
        The method returns the index of the current target and
        and an updated target dictionary.
        :param alignment: an alignment child of a BLAST record object
        :param targets: dictionary caching the known targets
        :return: current_target (ID in the database), targets
        """

        if alignment.accession in targets:
            current_target = targets[alignment.accession][0]
            if targets[alignment.accession][1] is False:
                self.session.query(Target).filter(
                    Target.target_name == alignment.accession).\
                    update({"target_length": alignment.length})
                self.session.commit()
                targets[alignment.accession] = (targets[alignment.accession][0],
                                                True)
        else:
            current_target = Target(alignment.accession,
                                    alignment.length)
            self.logger.warn("%s not found among targets, adding to the DB now",
                             alignment.accession)
            self.session.add(current_target)
            self.session.commit()
            assert isinstance(current_target.target_id, int)
            targets[alignment.accession] = (current_target.target_id, True)
            current_target = current_target.target_id
        return current_target, targets

    def load_into_db(self, objects):
        """
        :param objects: Objects to be loaded into the database
        :type objects: list

        Method to perform the bulk loading of objects into the SQL database.

        """

        try:
            self.session.bulk_save_objects(objects, return_defaults=False)
            self.session.commit()
        except sqlalchemy.exc.IntegrityError as err:
            self.logger.error('Database corrupted')
            self.logger.error(err)
            self.logger.error('Dropping and reloading')
            self.session.rollback()
            self.session.query(Hsp).delete()
            self.session.query(Hit).delete()
            self.session.bulk_save_objects(objects, return_defaults=False)
            self.session.commit()
# pylint: enable=too-many-instance-attributes
