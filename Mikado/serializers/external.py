# coding: utf-8

"""
This module is necessary to serialise external data.
The initial model is very simple - each external data
will have a tag, and internally the files will be TAB-delimited
2-column files, one for the TID and one for the score.
"""

import os
import pyfaidx
import sqlite3
from sqlalchemy import Column, String, Integer, ForeignKey, CHAR, Index, Float, Boolean
from sqlalchemy.sql.schema import PrimaryKeyConstraint
from sqlalchemy.ext.hybrid import hybrid_property
import sqlalchemy.exc
from sqlalchemy.orm import relationship, backref, column_property
from sqlalchemy.orm.session import Session  # sessionmaker
from sqlalchemy import select
from ..utilities.dbutils import DBBASE, Inspector, connect
from .blast_serializer import Query
from ..utilities.log_utils import create_null_logger, check_logger
from ..utilities.log_utils import check_logger, create_default_logger
from csv import DictReader


class ExternalSource(DBBASE):

    __tablename__ = "external_sources"

    source_id = Column(Integer, primary_key=True)
    source = Column(String)

    def __init__(self, source):

        self.source = source


class External(DBBASE):

    """This class serialises transcript ids, source of the score, and score itself."""

    __tablename__ = "external"

    query_id = Column(Integer, ForeignKey(Query.query_id), unique=False)
    source_id = Column(Integer, ForeignKey(ExternalSource.source_id), unique=False)
    ext_constraint = PrimaryKeyConstraint("query_id", "source_id", name="source_key")
    source = column_property(select([ExternalSource.source]).where(
        ExternalSource.source_id == source_id))
    score = Column(Float)

    query = column_property(select([Query.query_name]).where(
        Query.query_id == query_id))

    __table_args__ = (ext_constraint, )

    def __init__(self, query_id, source_id, score):

        self.query_id = query_id
        self.source_id = source_id
        self.score = score


class ExternalSerializer:

    def __init__(self, handle,
                 json_conf=None,
                 logger=None):

        """
        :param handle: the file to be serialized.
        :type handle: str | io.IOBase | io.TextIOWrapper

        :param json_conf: Optional configuration dictionary with db connection parameters.
        :type json_conf: dict | None
        """

        if logger is not None:
            logger = check_logger(logger)
            self.logger = logger
        else:
            self.logger = create_default_logger("external")

        self.handle = None

        if handle is None:
            self.logger.warning("No input file specified. Exiting.")
            self.close()
            return

        self.handle = open(handle)
        fasta_index = json_conf["serialise"]["files"]["transcripts"]
        if isinstance(fasta_index, str):
            assert os.path.exists(fasta_index)
            self.fasta_index = pyfaidx.Fasta(fasta_index)
            # self.fasta_index = SeqIO.index(fasta_index, "fasta")
        elif fasta_index is None:
            exc = ValueError("A fasta index is needed for the serialization!")
            self.logger.exception(exc)
            return
        else:
            assert isinstance(fasta_index, pyfaidx.Fasta)
            self.fasta_index = fasta_index

        self.parser = DictReader(self.handle, delimiter="\t")
        if not len(self.parser.fieldnames) > 1:
            error = TypeError("Not enough fields specified for the external file! Header: {}".format(
                self.parser.fieldnames))
            self.logger.exception(error)
            raise error
        elif self.parser.fieldnames[0].lower() != "tid":
            error = ValueError("The first column of the external file does not specify transcript IDs!")
            self.logger.exception(error)
            raise error

        self.engine = connect(json_conf, logger=logger)

        session = Session(bind=self.engine, autocommit=True, autoflush=True)
        inspector = Inspector.from_engine(self.engine)
        if External.__tablename__ not in inspector.get_table_names():
            DBBASE.metadata.create_all(self.engine)  # @UndefinedVariable

        self.session = session
        if json_conf is not None:
            self.maxobjects = json_conf["serialise"]["max_objects"]
        else:
            self.maxobjects = 10000

    def serialize(self):
        """This function will load the external scores into the database."""

        # First we have to load the sources inside the first table

        sources = dict()
        self.session.begin(subtransactions=True)
        for source in self.parser.fieldnames[1:]:
            source = ExternalSource(source)
            self.session.add(source)
        self.session.commit()

        # Now retrieve the values from the dictionary
        cache = dict()
        for source in self.session.query(ExternalSource):
            sources[source.source] = source.source_id
            continue
        for record in self.session.query(Query):
            cache[record.query_name] = record.query_id

        done = 0
        cache = self.load_fasta(cache)

        # Now we can serialise the content of the table
        tid_column = self.parser.fieldnames[0]
        objects = []
        for row in self.parser:
            tid = cache[row[tid_column]]
            for source in sources:
                objects.append(External(tid, sources[source], float(row[source])))
            if len(objects) >= self.maxobjects:
                self.session.begin(subtransactions=True)
                done += len(objects)
                self.session.bulk_save_objects(objects)
                self.logger.warning("Serialised %d values", done)
                self.session.commit()
                objects = []
        self.session.begin(subtransactions=True)
        self.session.bulk_save_objects(objects)
        done += len(objects)
        self.session.commit()
        self.logger.warning("Finished serialising %d values", done)
        self.close()

    def __call__(self):
        """
        Alias for serialize()
        """

        self.serialize()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def close(self):
        """
        Closing method for with. It will close the handles, if they are defined.
        :return:
        """

        if self.handle is not None:
            self.handle.close()
        return

    def load_fasta(self, cache):
        """
        Private method to load data from the FASTA file into the database.
        :param cache: a dictionary which memoizes the IDs.
        :type cache: dict

        :return: the updated cache
        :rtype: dict

        """

        objects = []

        if self.fasta_index is not None:
            done = 0
            self.logger.info("%d entries already present in db, %d in the index",
                             len([fasta_key for fasta_key in self.fasta_index if
                                  fasta_key not in cache]),
                             len(self.fasta_index.keys()))
            found = set()
            for record in self.fasta_index.keys():
                if record in cache:
                    continue
                objects.append(Query(record, len(self.fasta_index[record])))
                assert record not in found, record
                found.add(record)
                if len(objects) >= self.maxobjects:
                    done += len(objects)
                    self.session.begin(subtransactions=True)
                    self.session.bulk_save_objects(objects)
                    self.session.commit()
                    self.logger.info(
                        "Loaded %d transcripts into query table", done)
                    objects = []

            done += len(objects)
            self.session.begin(subtransactions=True)
            self.session.bulk_save_objects(objects)
            self.session.commit()
            self.logger.info(
                "Finished loading %d transcripts into query table", done)
        return cache
