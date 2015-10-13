#!/usr/bin/env python3
# coding:utf-8

"""
This module contains the necessary classes for serialising and querying ORF data.
"""

import os
import sqlite3

from Bio import SeqIO
from sqlalchemy import Column, String, Integer, ForeignKey, CHAR, Index, Float, Boolean
import sqlalchemy.exc
from sqlalchemy.orm import relationship, backref, column_property
from sqlalchemy.orm.session import sessionmaker
from sqlalchemy import select
from mikado_lib.utilities.dbutils import DBBASE, Inspector, connect

from mikado_lib.parsers import bed12
from mikado_lib.serializers.blast_serializer import Query
from mikado_lib.utilities.log_utils import create_null_logger, check_logger


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

    def __str__(self):
        return "{chrom}\t{start}\t{end}".format(
            chrom=self.query,
            start=self.start,
            end=self.end
        )

    @classmethod
    def as_bed12_static(cls, state, query_name):
        """Class method to transform the mapper into a BED12 object.
        Usable from outside the class."""
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

    # @hybrid_property
    # def query(self):
    #     """
    #     This property returns the name column value of the corresponding Query object.
    #     """
    #
    #     return self.query_object.query_name
    #
    # @query.expression
    # def query(cls):
    #     return Query.query_name


class OrfSerializer:
    """
    This class has the purpose of automating the loading of ORF information into the SQL database.
    """

    logger = create_null_logger("__orf_serializer__")

    def __init__(self, handle, fasta_index=None,
                 maxobjects=1000000, json_conf=None, logger=None):

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
        :type handle: io.IOBase
        :type handle: str

        :param fasta_index: a dictionary-like BioPython index object
        :type fasta_index: Bio.File._IndexedSeqFileDict

        :param maxobjects: Maximum number of entries to cache before bulk loading. Default: 10^6
        :type maxobjects: int

        :param json_conf: a configuration dictionary
        :type json_conf: dict

        """

        if logger is not None:
            self.logger = check_logger(logger)

        if isinstance(fasta_index, str):
            assert os.path.exists(fasta_index)
            self.fasta_index = SeqIO.index(fasta_index, "fasta")
        elif fasta_index is None:
            exc = ValueError("A fasta index is needed for the serialization!")
            self.logger.exception(exc)
            return
        else:
            assert "SeqIO.index" in repr(fasta_index)
            self.fasta_index = fasta_index

        self.bed12_parser = bed12.Bed12Parser(handle,
                                              fasta_index=fasta_index,
                                              transcriptomic=True)
        self.engine = connect(json_conf, logger)

        session = sessionmaker()
        session.configure(bind=self.engine)

        inspector = Inspector.from_engine(self.engine)
        if Orf.__tablename__ not in inspector.get_table_names():
            DBBASE.metadata.create_all(self.engine)
        self.session = session()
        self.maxobjects = maxobjects

    def serialize(self):
        """
        This method performs the parsing of the ORF file and the
        loading into the SQL database.
        """

        objects = []
        # Dictionary to hold the data before bulk loading into the database
        cache = dict()

        for record in self.session.query(Query):
            cache[record.query_name] = record.query_id

        done = 0
        if self.fasta_index is not None:
            self.logger.info("%d entries to load", len(self.fasta_index))
            self.logger.info("%d entries already present in db",
                             len([fasta_key for fasta_key in self.fasta_index if
                                  fasta_key not in cache]))
            found = set()
            for record in self.fasta_index:
                if record in cache:
                    continue
                objects.append(Query(record, len(self.fasta_index[record])))
                self.logger.debug("Appended %s", record)
                assert record not in found, record
                found.add(record)
                if len(objects) >= self.maxobjects:
                    done += len(objects)
                    self.session.bulk_save_objects(objects)
                    self.logger.info(
                        "Loaded %d transcripts into query table", done)
                    objects = []

            done += len(objects)
            self.session.bulk_save_objects(objects)

            self.session.commit()
            self.logger.info(
                "Finished loading %d transcripts into query table", done)
            objects = []
            done = 0

        self.logger.info("Loading IDs into the cache")
        for record in self.session.query(Query):
            cache[record.query_name] = record.query_id
        self.logger.info("Finished loading IDs into the cache")

        for row in self.bed12_parser:
            if row.header is True:
                continue
            if row.invalid is True:
                self.logger.warn("Invalid entry: %s", row)
                continue
            if row.id in cache:
                current_query = cache[row.id]
            else:
                current_query = Query(row.id, row.end)
                self.session.add(current_query)
                self.session.commit()
                cache[current_query.query_name] = current_query.query_id
                current_query = current_query.query_id

            current_junction = Orf(row, current_query)
            objects.append(current_junction)
            if len(objects) >= self.maxobjects:
                done += len(objects)
                self.session.bulk_save_objects(objects)
                self.logger.info("Loaded %d ORFs into the database", done)
                objects = []

        done += len(objects)
        self.session.bulk_save_objects(objects)
        self.logger.info("Finished loading %d ORFs into the database", done)
        self.session.commit()

    def __call__(self):
        """
        Alias for serialize
        """

        try:
            self.serialize()
        except (sqlalchemy.exc.IntegrityError, sqlite3.IntegrityError) as exc:
            self.logger.error("DB corrupted, reloading data. Error: %s",
                              exc)
            self.session.query(Query).delete()
            self.session.query(Orf).delete()
            self.serialize()
