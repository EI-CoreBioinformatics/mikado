#!/usr/bin/env python3
# coding:utf-8

"""
This module contains the necessary classes for serialising and querying ORF data.
"""

import os
import sqlite3
import pysam
from sqlalchemy import Column, String, Integer, ForeignKey, CHAR, Index, Float, Boolean
import sqlalchemy.exc
from sqlalchemy.orm import relationship, backref, column_property
from sqlalchemy.orm.session import sessionmaker
from sqlalchemy import select
from ..utilities.dbutils import DBBASE, Inspector, connect
from ..parsers import bed12  # , GFF
from .blast_serializer import Query
from ..utilities.log_utils import create_null_logger, check_logger
import pandas as pd
from ..exceptions import InvalidSerialization
import logging
import logging.handlers as logging_handlers
import multiprocessing as mp
import msgpack
import zlib


# This is a serialization class, it must have a ton of attributes ...
# pylint: disable=too-many-instance-attributes


class Orf(DBBASE):

    """
    Serialization class for ORFs derived from BED12 files.
    """

    __tablename__ = "orf"

    orf_id = Column(Integer, primary_key=True)
    query_id = Column(Integer, ForeignKey(Query.query_id), unique=False)
    start = Column(Integer, nullable=False)
    end = Column(Integer, nullable=False)
    orf_name = Column(String(200))
    strand = Column(CHAR)
    thick_start = Column(Integer, nullable=False)
    thick_end = Column(Integer, nullable=False)
    score = Column(Float)
    has_start_codon = Column(Boolean, nullable=True)
    has_stop_codon = Column(Boolean, nullable=True)
    cds_len = Column(Integer)
    phase = Column(Integer, nullable=False)

    __table_args__ = (Index("orf_index", "query_id", "thick_start", "thick_end"),
                      Index("query_index", "query_id"))

    query_object = relationship(Query, uselist=False,
                                backref=backref("orfs"), lazy="immediate")

    query = column_property(select([Query.query_name]).where(
        Query.query_id == query_id))

    def __init__(self, bed12_object, query_id):
        if not isinstance(bed12_object, bed12.BED12):
            raise TypeError("Invalid data type!")
        self.query_id = query_id
        self.start = bed12_object.start
        self.end = bed12_object.end
        self.thick_start = bed12_object.thick_start
        self.thick_end = bed12_object.thick_end
        self.orf_name = bed12_object.name
        self.strand = bed12_object.strand
        self.score = bed12_object.score
        self.has_start_codon = bed12_object.has_start_codon
        self.has_stop_codon = bed12_object.has_stop_codon
        self.cds_len = bed12_object.cds_len
        self.phase = bed12_object.phase

    def __str__(self):
        return "{chrom}\t{start}\t{end}".format(
            chrom=self.query,
            start=self.start,
            end=self.end
        )

    def as_dict(self):

        return {
            "start": self.start,
            "end": self.end,
            "orf_name": self.orf_name,
            "strand": self.strand,
            "thick_start": self.thick_start,
            "thick_end": self.thick_end,
            "score": self.score,
            "has_start_codon": self.has_start_codon,
            "has_stop_codon": self.has_stop_codon,
            "cds_len": self.cds_len,
            "phase": self.phase
        }

    @staticmethod
    def create_dict(bed12_object, query_id):
        if bed12_object.header is False and bed12_object.start is None:
            raise ValueError("Invalid BED! {}".format(bed12_object))

        obj = bed12_object.as_simple_dict()
        obj["query_id"] = query_id
        return obj
        # return {
        #     "query_id": query_id,
        #     "start": bed12_object.start,
        #     "end": bed12_object.end,
        #     "orf_name": bed12_object.name,
        #     "strand": bed12_object.strand,
        #     "thick_start": bed12_object.thick_start,
        #     "thick_end": bed12_object.thick_end,
        #     "score": bed12_object.score,
        #     "has_start_codon": bed12_object.has_start_codon,
        #     "has_stop_codon": bed12_object.has_stop_codon,
        #     "cds_len": bed12_object.cds_len,
        #     "phase": bed12_object.phase
        # }

    @classmethod
    def as_bed12_static(cls, state, query_name):
        """Class method to transform the mapper into a BED12 object.
        Usable from outside the class.

        :param state: the original state derived from the mapping.

        :param query_name: the name of the query, retrieved from the Query associated object/table.
        """
        __bed12 = bed12.BED12()

        __bed12.header = False
        __bed12.query = __bed12.chrom = query_name
        __bed12.start = state.start
        __bed12.end = state.end
        __bed12.name = state.orf_name
        __bed12.score = state.score
        __bed12.strand = state.strand
        __bed12.thick_start = state.thick_start
        __bed12.thick_end = state.thick_end
        __bed12.rgb = 0
        __bed12.block_count = 1
        __bed12.block_sizes = [state.end]
        __bed12.block_starts = [0]
        __bed12.transcriptomic = True
        __bed12.phase = state.phase

        # Verbose block, but it is necessary as raw extraction from SQL
        # yields 0/1 instead of True/False
        if state.has_start_codon:
            __bed12.has_start_codon = True
        else:
            __bed12.has_start_codon = False
        if state.has_stop_codon:
            __bed12.has_stop_codon = True
        else:
            __bed12.has_stop_codon = False

        return __bed12

    def as_bed12(self):
        """Method to transform the mapper into a BED12 object."""

        return self.as_bed12_static(self, self.query)


def line_parser_func(handle, fai, send_queue):
    fai = pysam.FastaFile(fai)
    for num, line in enumerate(open(handle)):
        if line[0] == "#":
            send_queue.put((num, line, None))
        else:
            _f = line.split("\t")
            if _f[0] not in fai:
                seq = None
            else:
                seq = zlib.compress(fai[line.split("\t")[0]].encode(), 1)
            send_queue.put_nowait((num, line, seq))

    send_queue.put("EXIT")


class OrfSerializer:
    """
    This class has the purpose of automating the loading of ORF information into the SQL database.
    """

    logger = create_null_logger("__orf_serializer__")

    def __init__(self,
                 handle,
                 json_conf=None,
                 logger=None):

        """Constructor function. Arguments:
        - handle         the BED12 file
        - db             Output DB
        - fasta_index    A SeqIO-like index of sequence records.
        Alternatively, the path to the FASTA file. REQUIRED.
        - maxobjects    Integer. Indicates how big should the cache be
        for objects to be loaded inside the DB

        It is HIGHLY RECOMMENDED to provide the fasta index,
        as it will make the population of the Query table much faster.

        :param handle: the input BED12 file
        :type handle: (io.TextIOWrapper|str)

        :param json_conf: a configuration dictionary
        :type json_conf: dict

        """

        if logger is not None:
            self.logger = check_logger(logger)

        fasta_index = json_conf["serialise"]["files"]["transcripts"]
        self._max_regression = json_conf["serialise"]["max_regression"]
        self._table = json_conf["serialise"]["codon_table"]
        self.procs = json_conf["threads"]
        self.single_thread = json_conf["serialise"]["single_thread"]
        self.adjust_start = json_conf["serialise"]["start_adjustment"]
        if self.single_thread:
            self.procs = 1

        if isinstance(fasta_index, str):
            assert os.path.exists(fasta_index)
            self.fasta_index = pysam.FastaFile(fasta_index)
        elif fasta_index is None:
            exc = ValueError("A fasta index is needed for the serialization!")
            self.logger.exception(exc)
            return
        else:
            assert isinstance(fasta_index, pysam.FastaFile)
            self.fasta_index = fasta_index

        if isinstance(handle, str):
            self.is_bed12 = (".bed12" in handle or ".bed" in handle)
        else:
            self.is_bed12 = (".bed12" in handle.name or ".bed" in handle.name.endswith)

        self.engine = connect(json_conf, logger)

        self._handle = handle
        Session = sessionmaker(bind=self.engine, autocommit=False, autoflush=False, expire_on_commit=False)
        session = Session()
        # session.configure(bind=self.engine)

        inspector = Inspector.from_engine(self.engine)
        if Orf.__tablename__ not in inspector.get_table_names():
            DBBASE.metadata.create_all(self.engine)
        self.session = session
        self.maxobjects = json_conf["serialise"]["max_objects"]
        self.log_level = json_conf["log_settings"]["log_level"]

    def load_fasta(self):

        """
        Private method to load data from the FASTA file into the database.
        :param cache: a dictionary which memoizes the IDs.
        :type cache: dict

        :return: the updated cache
        :rtype: dict

        """

        objects = []

        cache = pd.read_sql_table("query", self.engine, index_col="query_name", columns=["query_name", "query_id"])
        cache = cache.to_dict()["query_id"]

        if self.fasta_index is not None:
            done = 0
            self.logger.debug("%d entries already present in db, %d in the index",
                             len([fasta_key for fasta_key in self.fasta_index.references if
                                  fasta_key not in cache]),
                             self.fasta_index.nreferences)
            found = set()
            for ref, length in zip(self.fasta_index.references, self.fasta_index.lengths):
                if ref in cache:
                    continue
                objects.append({"query_name": ref, "query_length": length})
                # objects.append(Query(ref, length))
                assert ref not in found, ref
                found.add(ref)
                if len(objects) >= self.maxobjects:
                    done += len(objects)
                    self.engine.execute(
                        Query.__table__.insert(),
                        objects
                    )
                    # self.session.bulk_save_objects(objects)
                    self.session.commit()
                    self.logger.debug(
                        "Loaded %d transcripts into query table", done)
                    objects = []

            done += len(objects)
            self.engine.execute(
                Query.__table__.insert(),
                objects
            )
            # self.session.begin(subtransactions=True)
            # self.session.bulk_save_objects(objects)
            self.session.commit()
            self.logger.debug("Finished loading %d transcripts into query table", done)
        return

    def __serialize_single_thread(self):

        self.bed12_parser = bed12.Bed12Parser(self._handle,
                                              fasta_index=self.fasta_index,
                                              logger=self.logger,
                                              is_gff=(not self.is_bed12),
                                              transcriptomic=True,
                                              max_regression=self._max_regression,
                                              table=self._table)
        objects = []
        done = 0

        not_found = set()
        for row in self.bed12_parser:
            if row.header is True:
                continue
            if row.invalid is True:
                self.logger.warning("Invalid entry, reason: %s\n%s",
                                    row.invalid_reason,
                                    row)
                continue
            if row.id in self.query_cache:
                current_query = self.query_cache[row.id]
            elif not self.initial_cache:
                current_query = Query(row.id, row.end)
                not_found.add(row.id)
                self.session.add(current_query)
                self.session.commit()
                self.query_cache[current_query.query_name] = current_query.query_id
                current_query = current_query.query_id
            else:
                self.logger.critical(
                    "The provided ORFs do not match the transcripts provided and already present in the database.\
Please check your input files.")
                raise InvalidSerialization

            # current_junction = Orf(row, current_query)
            obj = Orf.create_dict(row, current_query)
            if obj["start"] is None or not isinstance(obj["start"], int):
                raise ValueError("Invalid object: {}".format(obj))
                # continue
            objects.append(Orf.create_dict(row, current_query))
            if len(objects) >= self.maxobjects:
                done += len(objects)
                self.session.begin(subtransactions=True)
                # self.session.bulk_save_objects(objects)
                self.engine.execute(
                    Orf.__table__.insert(),
                    objects
                )
                self.session.commit()
                self.logger.debug("Loaded %d ORFs into the database", done)
                objects = []

        done += len(objects)
        # self.session.begin(subtransactions=True)
        # self.session.bulk_save_objects(objects, update_changed_only=False)
        self.engine.execute(
            Orf.__table__.insert(),
            objects
        )
        self.session.commit()
        self.session.close()
        self.logger.info("Finished loading %d ORFs into the database", done)

        orfs = pd.read_sql_table("orf", self.engine, index_col="query_id")
        if orfs.shape[0] != done:
            raise ValueError("I should have serialised {} ORFs, but {} are present!".format(done, orfs.shape[0]))

    def __serialize_multiple_threads(self):
        """"""

        manager = mp.Manager()
        send_queue = manager.Queue(-1)
        return_queue = manager.JoinableQueue(-1)
        self.logging_queue = mp.Queue(-1)
        self.logger_queue_handler = logging_handlers.QueueHandler(self.logging_queue)
        self.queue_logger = logging.getLogger("parser")
        self.queue_logger.addHandler(self.logger_queue_handler)
        self.queue_logger.setLevel(self.log_level)
        self.queue_logger.propagate = False
        self.log_writer = logging_handlers.QueueListener(self.logging_queue, self.logger)
        self.log_writer.start()

        line_parser = mp.Process(target=line_parser_func,
                                 args=(self._handle, self.fasta_index.filename, send_queue))
        line_parser.start()

        parsers = [bed12.Bed12ParseWrapper(
            identifier=index,
            rec_queue=send_queue,
            log_queue=self.logging_queue,
            level=self.log_level,
            return_queue=return_queue,
            fasta_index=None,
            is_gff=(not self.is_bed12),
            transcriptomic=True,
            max_regression=self._max_regression,
            table=self._table) for index in range(self.procs)]
        [_.start() for _ in parsers]

        not_found = set()
        done = 0
        objects = []
        procs_done = 0
        while True:
            num = return_queue.get()
            if num in ("FINISHED", b"FINISHED"):
                procs_done += 1
                if procs_done == self.procs:
                    break
                else:
                    continue
            num, obj = num
            try:
                object = msgpack.loads(obj, raw=False)
            except TypeError:
                raise TypeError(obj)

            if object["id"] in self.query_cache:
                current_query = self.query_cache[object["id"]]
            elif not self.initial_cache:
                current_query = Query(object["id"], object["end"])
                not_found.add(object["id"])
                self.session.add(current_query)
                self.session.commit()
                self.query_cache[current_query.query_name] = current_query.query_id
                current_query = current_query.query_id
            else:
                self.logger.critical(
                    "The provided ORFs do not match the transcripts provided and already present in the database.\
This could be due to having called the ORFs on a FASTA file different from `mikado_prepared.fasta`, the output of \
mikado prepare. If this is the case, please use mikado_prepared.fasta to call the ORFs and then restart \
`mikado serialise` using them as input.")
                raise InvalidSerialization

            object["query_id"] = current_query
            objects.append(object)
            if len(objects) >= self.maxobjects:
                done += len(objects)
                self.session.begin(subtransactions=True)
                self.engine.execute(
                    Orf.__table__.insert(),
                    objects
                )
                self.session.commit()
                self.logger.debug("Loaded %d ORFs into the database", done)
                objects = []

        [proc.join() for proc in parsers]
        done += len(objects)
        # self.session.begin(subtransactions=True)
        # self.session.bulk_save_objects(objects, update_changed_only=False)
        if objects:
            self.engine.execute(
                Orf.__table__.insert(),
                objects
            )
        self.session.commit()
        self.session.close()
        self.logger.info("Finished loading %d ORFs into the database", done)

        manager.shutdown()
        orfs = pd.read_sql_table("orf", self.engine, index_col="query_id")
        if orfs.shape[0] != done:
            raise ValueError("I should have serialised {} ORFs, but {} are present!".format(done, orfs.shape[0]))

    def serialize(self):
        """
        This method performs the parsing of the ORF file and the
        loading into the SQL database.
        """

        self.load_fasta()
        self.query_cache = pd.read_sql_table("query", self.engine, index_col="query_name", columns=["query_name", "query_id"])
        self.query_cache = self.query_cache.to_dict()["query_id"]
        self.initial_cache = (len(self.query_cache) > 0)

        if self.procs == 1:
            self.__serialize_single_thread()
        else:
            try:
                self.__serialize_multiple_threads()
            finally:
                pass

    def __call__(self):
        """
        Alias for serialize
        """

        try:
            [idx.drop(bind=self.engine) for idx in Orf.__table__.indexes]
        except (sqlalchemy.exc.IntegrityError, sqlite3.IntegrityError) as exc:
            self.logger.debug("Corrupt table found, deleting and restarting")
            self.session.query(Orf).delete()
        try:
            self.serialize()
        except (sqlalchemy.exc.IntegrityError, sqlite3.IntegrityError) as exc:
            self.logger.error("DB corrupted, reloading data. Error: %s",
                              exc)
            self.session.query(Query).delete()
            self.session.query(Orf).delete()
            try:
                self.serialize()
            except InvalidSerialization:
                raise
        finally:
            [idx.create(bind=self.engine) for idx in Orf.__table__.indexes]
