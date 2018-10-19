# coding: utf-8

"""
This module is necessary to serialise external data.
The initial model is very simple - each external data
will have a tag, and internally the files will be TAB-delimited
2-column files, one for the TID and one for the score.
"""

import os
import pyfaidx
from sqlalchemy import Column, String, Integer, ForeignKey, Float
from sqlalchemy.sql.schema import PrimaryKeyConstraint
from sqlalchemy.orm import column_property
from sqlalchemy.orm.session import Session  # sessionmaker
from sqlalchemy import select
from ..utilities.dbutils import DBBASE, Inspector, connect
from .blast_serializer import Query
from ..utilities.log_utils import check_logger, create_default_logger
from csv import DictReader
import numbers
import sqlalchemy.exc
import pandas as pd


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
        if not isinstance(score, numbers.Number):
            raise sqlalchemy.exc.ArgumentError("This class only accepts numeric scores")
        self.score = score


class ExternalSerializer:

    def __init__(self, handle,
                 json_conf=None,
                 logger=None, delimiter="\t"):

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

        self.handle = handle

        if handle is None:
            self.logger.warning("No input file specified. Exiting.")
            self.close()
            return

        fasta_index = json_conf["serialise"]["files"]["transcripts"]
        if isinstance(fasta_index, str):
            assert os.path.exists(fasta_index)
            self.fasta_index = pyfaidx.Fasta(fasta_index)
            # self.fasta_index = SeqIO.index(fasta_index, "fasta")
        elif fasta_index is None:
            self.logger.debug("No fasta index provided, we presume transcripts are already in the DB.")
            # exc = ValueError("A fasta index is needed for the serialization!")
            # self.logger.exception(exc)
            self.fasta_index = None
        else:
            assert isinstance(fasta_index, pyfaidx.Fasta)
            self.fasta_index = fasta_index

        try:
            self.data = pd.read_csv(self.handle, delimiter=delimiter, index_col=["tid"])
        except ValueError:
            self.data = pd.read_csv(self.handle, delimiter=delimiter)
            if self.data.index.name is not None:
                exc = ValueError("Invalid index name: {}. It should be tid or not specified.".format(
                    self.data.index.name))
                self.logger.exception(exc)
                return
            else:
                self.data.index.name = "tid"
        except KeyboardInterrupt:
            raise KeyboardInterrupt
        except Exception as exc:
            self.logger.exception(exc)
            raise

        if len(self.data.columns) == 0:
            error = TypeError("Not enough fields specified for the external file! Header: {}".format(
                self.data.columns))
            self.logger.exception(error)
            raise error

        self.data.fillna(0, inplace=True)
        for column in self.data.columns:
            try:
                self.data[column].astype("float")
            except ValueError:
                exc = ValueError("Invalid non-numeric values in external table, for column {}. Aborting".format(
                    column
                ))
                self.logger.critical(exc)
                raise

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
        for source in self.data.columns:
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
        objects = []
        decimal_place = 1
        for source in sources:
            for idx in self.data.index:
                objects.append(External(query_id=cache[idx],
                                        source_id=sources[source],
                                        score=self.data[source][idx]))
            if len(objects) >= self.maxobjects:
                self.session.begin(subtransactions=True)
                done += len(objects)
                self.session.bulk_save_objects(objects)
                self.logger.debug("Serialised %d values", done)
                self.session.commit()
                objects = []
                if done >= self.data.shape[1] * 0.1 * decimal_place:
                    self.logger.info("Serialised %s%% of the external values",
                                     decimal_place * 10)

        self.session.begin(subtransactions=True)
        self.session.bulk_save_objects(objects)
        done += len(objects)
        self.session.commit()
        self.logger.info("Finished serialising %d values", done)
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

        # Now update cache
        cache.update(
            dict((_.query_name, _.query_id) for _ in iter(self.session.query(Query)))
        )

        assert len(cache) > 0

        return cache
