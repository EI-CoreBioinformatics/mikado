#!/usr/bin/env python3
# coding:utf-8

"""
This module contains the necessary classes for serialising and querying ORF data.
"""

import os
import logging
from Bio import SeqIO
from sqlalchemy import Column, String, Integer, ForeignKey, CHAR, Index, Float, Boolean
from sqlalchemy.ext.hybrid import hybrid_property
from mikado_lib.parsers import bed12
from sqlalchemy.orm import relationship, backref
from sqlalchemy.engine import create_engine
from sqlalchemy.orm.session import sessionmaker
from mikado_lib.serializers.dbutils import dbBase, Inspector, connect
from mikado_lib.serializers.blast_utils import Query


class Orf(dbBase):

    """
    Serialization class for ORFs derived from BED12 files.


    :param id: numerical key index
    :type orf_id: int

    :param query_id: foreign key id for the blast Query table
    :type query_id: int

    :param start: start position of the ORF. It should always be 1
    :type start: int

    :param end: end position of the ORF. It is equivalent to the transcript length.
    :type end: int

    :param name: Name of the ORF.
    :type name: str

    :param strand: Strand of the ORF. One of "+","-"
    :type strand: str

    :param thickStart: ORF start position on the transcript.
    :type thickStart: int

    :param thickEnd: ORF end position on the transcript.
    :type thickEnd: int

    :param score: Score assigned to this ORF
    :type score: float

    :param has_start_codon: boolean flag that indicates whether a start codon is present at the beginning of the ORF.
    :type has_start_codon: bool

    :param has_stop_codon: boolean flag that indicates whether a termination codon is present at the end of the ORF.
    :type has_stop_codon: bool

    :param cds_len: length of the ORF
    :type cds_len: int

    :param query_object: a Query object which is referenced inside the BLAST tables.
    :type query_object: mikado_lib.serializers.blast_utils.Query

    """

    __tablename__ = "orf"

    orf_id = Column(Integer, primary_key=True)
    query_id = Column(Integer, ForeignKey(Query.query_id), unique=False)
    start = Column(Integer, nullable=False)
    end = Column(Integer, nullable=False)
    orf_name = Column(String(200))
    strand = Column(CHAR)
    thickStart = Column(Integer, nullable=False)
    thickEnd = Column(Integer, nullable=False)
    score = Column(Float)
    has_start_codon = Column(Boolean, nullable=True)
    has_stop_codon = Column(Boolean, nullable=True)
    cds_len = Column(Integer)

    __table_args__ = (Index("orf_index", "query_id", "thickStart", "thickEnd"), Index("query_index", "query_id"))

    query_object = relationship(Query, uselist=False, backref=backref("orfs"), lazy="immediate")

    def __init__(self, bed12_object, query_id):
        if type(bed12_object) is not bed12.BED12:
            raise TypeError("Invalid data type!")
        self.query_id = query_id
        self.start = bed12_object.start
        self.end = bed12_object.end
        self.thickStart = bed12_object.thickStart
        self.thickEnd = bed12_object.thickEnd
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
        b = bed12.BED12()

        b.header = False
        b.query = b.chrom = query_name
        b.start = state.start
        b.end = state.end
        b.name = state.orf_name
        b.score = state.score
        b.strand = state.strand
        b.thickStart = state.thickStart
        b.thickEnd = state.thickEnd
        b.rgb = 0
        b.blockCount = 1
        b.blockSizes = [state.end]
        b.blockStarts = [0]

        # Verbose block, but it is necessary as raw extraction from SQL
        # yields 0/1 instead of True/False
        if state.has_start_codon:
            b.has_start_codon = True
        else:
            b.has_start_codon = False
        if state.has_stop_codon:
            b.has_stop_codon = True
        else:
            b.has_stop_codon = False

        return b

    def as_bed12(self):
        """Method to transform the mapper into a BED12 object."""

        return self.as_bed12_static(self, self.query)

    @hybrid_property
    def query(self):
        """
        This property returns the name column value of the corresponding Query object.
        """

        return self.query_object.query_name


class OrfSerializer:
    """
    This class has the purpose of automating the loading of ORF information into the SQL database.
    """

    def __init__(self, handle, db=None, fasta_index=None, maxobjects=1000000, json_conf=None):

        """Constructor function. Arguments:
        - handle         the BED12 file
        - db             Output DB
        - fasta_index    A SeqIO-like index of sequence records. Alternatively, the path to the FASTA file. REQUIRED.
        - maxobjects    Integer. Indicates how big should the cache be for objects to be loaded inside the DB
        
        It is HIGHLY RECOMMENDED to provide the fasta index, as it will make the population of the Query
        table much faster.

        :param handle: the input BED12 file
        :type handle: io.IOBase
        :type handle: str

        :param db: the database to be used as output
        :type db: str
        :type db: None

        :param fasta_index: a dictionary-like BioPython index object
        :type fasta_index: Bio.File._IndexedSeqFileDict

        :param maxobjects: Maximum number of entries to cache before bulk loading. Default: 10^6
        :type maxobjects: int

        :param json_conf: a configuration dictionary
        :type json_conf: dict

        """
        self.logger = logging.getLogger("main")
        self.logger.setLevel(logging.INFO)
        self.handler = logging.StreamHandler()
        self.formatter = logging.Formatter("{asctime} - {levelname} - {message}", style='{')
        self.handler.setFormatter(self.formatter)
        self.logger.addHandler(self.handler)

        if type(fasta_index) is str:
            assert os.path.exists(fasta_index)
            self.fasta_index = SeqIO.index(fasta_index, "fasta")
        elif fasta_index is None:
            raise ValueError("A fasta index is needed for the serialization!")
        else:
            assert "SeqIO.index" in repr(fasta_index)
            self.fasta_index = fasta_index

        self.BED12 = bed12.Bed12Parser(handle, fasta_index=fasta_index, transcriptomic=True)
        if json_conf is not None:
            self.engine = connect(json_conf)
        else:
            self.engine = create_engine("sqlite:///{0}".format(db))

        session = sessionmaker()
        session.configure(bind=self.engine)

        inspector = Inspector.from_engine(self.engine)
        if Orf.__tablename__ not in inspector.get_table_names():
            dbBase.metadata.create_all(self.engine)  # @UndefinedVariable
        self.session = session()
        self.maxobjects = maxobjects

    def serialize(self):
        """
        This method performs the parsing of the ORF file and the loading into the SQL database.
        """

        objects = []
        cache = dict()  # Dictionary to hold the data before bulk loading into the database

        for record in self.session.query(Query):
            cache[record.query_name] = record.query_id

        done = 0
        if self.fasta_index is not None:
            self.logger.info("{0} entries to load".format(len(self.fasta_index)))
            self.logger.info("{0} entries already present in db".format(
                len(list(filter(lambda rr: rr not in cache, self.fasta_index.keys())))
            ))
            found = set()
            for record in self.fasta_index:
                objects.append(Query(record, len(self.fasta_index[record])))
                self.logger.debug("Appended {0}".format(record))
                assert record not in found, record
                found.add(record)
                if len(objects) >= self.maxobjects:
                    done += len(objects)
                    self.logger.info("Loaded {0} transcripts into query table".format(done))
                    self.session.bulk_save_objects(objects)
                    objects = []

            done += len(objects)
            self.logger.info("Finished loading {0} transcripts into query table".format(done))
            self.session.bulk_save_objects(objects)
            self.session.commit()
            objects = []
            done = 0

        self.logger.info("Loading IDs into the cache")
        for record in self.session.query(Query):
            cache[record.query_name] = record.query_id
        self.logger.info("Finished loading IDs into the cache")

        for row in self.BED12:
            if row.header is True:
                continue
            if row.invalid is True:
                self.logger.warn("Invalid entry: {0}".format(row))
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
                self.logger.info("Loaded {0} ORFs into the database".format(done))
                self.session.bulk_save_objects(objects)
                objects = []

        done += len(objects)
        self.logger.info("Finished loading {0} ORFs into the database".format(done))

        self.session.bulk_save_objects(objects)
        self.session.commit()

    def __call__(self):
        """
        Alias for serialize
        """

        self.serialize()
