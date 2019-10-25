# coding: utf-8

"""
This module is necessary to serialise external data.
The initial model is very simple - each external data
will have a tag, and internally the files will be TAB-delimited
2-column files, one for the TID and one for the score.
"""

import os
import pyfaidx
from sqlalchemy import Column, String, Integer, ForeignKey, Float, Boolean
from sqlalchemy.sql.schema import PrimaryKeyConstraint
from sqlalchemy.orm import column_property
from sqlalchemy.orm.session import Session  # sessionmaker
from sqlalchemy import select
from ..utilities.dbutils import DBBASE, Inspector, connect
from .blast_serializer import Query
from ..utilities.log_utils import check_logger, create_default_logger
import numbers
import pandas as pd
import numpy as np
import re
from collections import Counter
import sqlalchemy.exc


class ExternalSource(DBBASE):

    __tablename__ = "external_sources"

    source_id = Column(Integer, primary_key=True)
    source = Column(String, unique=True)
    rtype = Column(String, unique=False)
    valid_raw = Column(Boolean)

    def __init__(self, source, rtype, valid_raw):

        self.source = source
        if valid_raw not in (True, False):
            raise ValueError("\"Valid raw\" flags must be boolean!")
        if np.dtype("bool") == rtype:
            rtype = "bool"
        elif np.dtype("int") == rtype:
            rtype = "int"
        elif np.dtype("float") == rtype:
            rtype = "float"
        elif np.dtype("complex") == rtype:
            rtype = "complex"
        else:
            raise ValueError("Invalid source rtype for {}: {}".format(source, rtype))

        self.rtype = rtype
        self.valid_raw = valid_raw


class External(DBBASE):

    """This class serialises transcript ids, source of the score, and score itself."""

    __tablename__ = "external"

    query_id = Column(Integer, ForeignKey(Query.query_id), unique=False)
    source_id = Column(Integer, ForeignKey(ExternalSource.source_id), unique=False)
    ext_constraint = PrimaryKeyConstraint("query_id", "source_id", name="source_key")
    source = column_property(select([ExternalSource.source]).where(
        ExternalSource.source_id == source_id))
    score = Column(String, nullable=False)

    query = column_property(select([Query.query_name]).where(
        Query.query_id == query_id))

    valid_raw = column_property(select([ExternalSource.valid_raw]).where(
        ExternalSource.source_id == source_id))

    rtype = column_property(select([ExternalSource.rtype]).where(
        ExternalSource.source_id == source_id))

    __table_args__ = (ext_constraint, )

    def __init__(self, query_id, source_id, score):

        self.query_id = query_id
        self.source_id = source_id
        if not isinstance(score, numbers.Number):
            raise sqlalchemy.exc.ArgumentError("Invalid score for external values: {}".format(type(score)))
        score = str(score)
        assert score.strip()
        self.score = str(score)


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
        if isinstance(fasta_index, (str, bytes)):
            if isinstance(fasta_index, bytes):
                fasta_index = fasta_index.decode()
            if not os.path.exists(fasta_index):
                error = """I cannot find the mikado prepared FASTA file with the transcripts to analyse.
                Please run mikado serialise in the folder with the correct files, and/or modify the configuration
                or the command line options."""
                self.logger.critical(error)
                raise AssertionError(error)
            self.fasta_index = pyfaidx.Fasta(fasta_index)
            # self.fasta_index = SeqIO.index(fasta_index, "fasta")
        elif fasta_index is None:
            self.logger.debug("No fasta index provided, we presume transcripts are already in the DB.")
            self.fasta_index = None
        elif isinstance(fasta_index, pyfaidx.Fasta):
            self.fasta_index = fasta_index
        else:
            error = "Unkwnown FASTA index: {}. I will presume that the transcripts are in the database".format(
                type(fasta_index))
            self.logger.warning(error)

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

        self.engine = connect(json_conf, logger=logger)

        session = Session(bind=self.engine, autocommit=False, autoflush=False, expire_on_commit=False)
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

        # Check columns
        cols = []
        for col in self.data.columns:
            cols.append(re.sub(r"\.[0-9]*", '', str(col)))

        cols = Counter(cols)

        if cols.most_common()[0][1] > 1:  # IE the most common element is present more than one time
            raise IndexError("Duplicated values in the external table: {}".format(
                ",".join([_[0] for _ in cols.most_common() if _[1] > 1])
            ))

        for source in self.data.columns:

            if ((self.data[source].dtype == np.dtype("float") or
                    self.data[source].dtype == np.dtype("int")) and
                    0 <= self.data[source].min() <= self.data[source].max() <= 1):
                valid_raw = True
            else:
                valid_raw = False

            rtype = self.data[source].dtype
            if rtype == np.dtype("bool"):
                # We have to cast booleans as integers otherwise the conversion after extraction will fail
                self.data[source] = self.data[source].astype(int)

            source = ExternalSource(source, rtype=rtype, valid_raw=valid_raw)
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
                done += len(objects)
                self.session.bulk_save_objects(objects)
                self.logger.debug("Serialised %d values", done)
                self.session.commit()
                objects = []
                if done >= self.data.shape[1] * self.data.shape[0] * 0.1 * decimal_place:
                    self.logger.info("Serialised %s%% of the external values",
                                     decimal_place * 10)
                    decimal_place += 1

        self.session.bulk_save_objects(objects, update_changed_only=False)
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

        if hasattr(self, "fasta_index") and self.fasta_index is not None:
            self.fasta_index.close()
        if hasattr(self, "session"):
            self.session.close()
        if hasattr(self, "engine"):
            self.engine.dispose()
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
