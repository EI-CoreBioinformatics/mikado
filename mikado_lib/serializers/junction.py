# coding: utf-8

"""
This module is necessary to serialise the junction information
provided in BED12 format by programs such as TopHat or portcullis.
Please convert any other (custom) format into BED12
before loading.
"""

import io
import os
from mikado_lib.parsers import bed12
from mikado_lib.serializers.dbutils import dbBase
from mikado_lib.serializers.dbutils import Inspector
from sqlalchemy import Column
from sqlalchemy import String
from sqlalchemy import Integer
from sqlalchemy import ForeignKey
from sqlalchemy import CHAR
from sqlalchemy import Index
from sqlalchemy import Float
from sqlalchemy.orm import relationship, backref
from sqlalchemy.engine import create_engine
from sqlalchemy.orm.session import sessionmaker
from sqlalchemy.ext.hybrid import hybrid_method,hybrid_property


class Chrom(dbBase):
    """
    Simple serialization class for chromosomes.

    :param id: numerical id (table key)
    :type id: int

    :param name: chromosome name
    :type name: str

    :param length: length of the chromosome
    :type length: None
    :type length: int
    """

    __tablename__ = "chrom"
    __table_args__ = {"extend_existing": True}

    id = Column(Integer, primary_key=True)
    name = Column(String(200))
    length = Column(Integer, nullable=True)

    def __init__(self, name, length=None):
        self.name = name
        if length is not None:
            assert type(length) is int
        self.length = length


class Junction(dbBase):
    """
    Class that describes the junction table in the database.

    :param id: numerical id
    :type id: int

    :param chrom_id: numerical foreign id for the Chrom table
    :type chrom_id: int

    :param start: start position of the junction (1-based)
    :type start: int

    :param end: end position of the junction (1-based)
    :type end: int

    :param name: Name of the junction
    :type name: str

    :param strand: Strand of the junction. One of "+","-"
    :type strand: str

    :param junctionStart: Start internal position of the junction i.e. intron start (1-based)
    :type junctionStart: int

    :param junctionEnd: End internal position of the junction i.e. intron end (1-based)
    :type junctionEnd: int

    :param score: Score of the junction.
    :type score: Float

    :param chrom_object: a reference to the Chrom object over which the junction is positioned.
    :type chrom_object: Chrom

    """

    __tablename__ = "junctions"

    id = Column(Integer, primary_key=True)
    chrom_id = Column(Integer, ForeignKey(Chrom.id), unique=False)
    start = Column(Integer, nullable=False)
    end = Column(Integer, nullable=False)
    name = Column(String(200))
    strand = Column(CHAR)
    junctionStart = Column(Integer, nullable=False)
    junctionEnd = Column(Integer, nullable=False)
    score = Column(Float)
    __table_args__ = (Index("junction_index", "chrom_id", "junctionStart", "junctionEnd"), {"extend_existing": True})

    chrom_object = relationship(Chrom, uselist=False, backref=backref("junctions"), lazy="immediate")

    def __init__(self, bed12_object, chrom_id):
        """
        Serialization initializer.
        :param bed12_object: a BED12-like class instance.
        :type bed12_object: mikado_lib.parsers.bed12.BED12

        :param chrom_id: the numerical foreign key for the Chrom table.
        :type chrom_id: int
        """

        if type(bed12_object) is not bed12.BED12:
            raise TypeError("Invalid data type!")
        self.chrom_id = chrom_id
        self.start = bed12_object.start
        self.end = bed12_object.end
        self.junctionStart = bed12_object.thickStart
        self.junctionEnd = bed12_object.thickEnd
        self.name = bed12_object.name
        self.strand = bed12_object.strand
        self.score = bed12_object.score

    def __str__(self):
        return "{chrom}\t{start}\t{end}".format(
            chrom=self.chrom,
            start=self.start,
            end=self.end
        )

    @hybrid_property
    def chrom(self):
        """
        This property returns the name of the chromosome upon which the junction is located.
        """

        return self.chrom_object.name

    @hybrid_method
    def is_equal(self, chrom, start, end, strand):
        """
        Function to verify whether a set of coordinates is equal to those in the DB.

        :param chrom: chromosome
        :type chrom: str

        :param start: start position (1-offset)
        :type start: int

        :param end: end position (1-offset)
        :type end: int

        :param strand: strand of the junction
        :type strand: str
        """
        return (self.chrom == chrom) and (self.start == start) and \
               (self.end == end) and (self.strand == strand)


class JunctionSerializer:
    """
    This class is used to serialize a junction BED12 file into an SQL database.
    """

    def __init__(self, handle, db, fai=None, maxobjects=10000, json_conf=None):

        """
        :param handle: the file to be serialized.
        :type handle: str | io.IOBase

        :param db: the database to connect to. Can be None.
        :type db: None | str

        :param fai: an optional FAI file to be used for generating the database.
        :type fai: str | io.IOBase

        :param maxobjects: maximal number of objects to cache in memory before bulk loading. Default 10^4
        :type maxobjects: int

        :param json_conf: Optional configuration dictionary with db connection parameters.
        :type json_conf: dict | None
        """

        self.BED12 = bed12.Bed12Parser(handle)
        if json_conf is not None:
            if json_conf["dbtype"] == "sqlite":
                self.engine = create_engine("sqlite:///{0}".format(json_conf["db"]))
            else:
                self.engine = create_engine("{dbtype}://{dbuser}:{dbpasswd}@{dbhost}/{db}".format(
                    dbtype=json_conf["dbtype"],
                    dbuser=json_conf["dbuser"],
                    dbpasswd=json_conf["dbpasswd"],
                    dbhost=json_conf["dbhost"],
                    db=json_conf["db"]))
        else:
            self.engine = create_engine("sqlite:///{0}".format(db))

        session = sessionmaker()
        session.configure(bind=self.engine)

        inspector = Inspector.from_engine(self.engine)
        if Junction.__tablename__ not in inspector.get_table_names():
            dbBase.metadata.create_all(self.engine)  # @UndefinedVariable

        self.session = session()
        self.maxobjects = maxobjects

        self.fai = fai
        if type(fai) is str:
            assert os.path.exists(fai)
            self.fai = open(self.fai)
        else:
            if fai is not None:
                assert type(fai) is io.TextIOWrapper
            self.fai = fai

    def serialize(self):
        """
        Workhorse of the class. It parses the input file and loads it into the database.
        """

        sequences = dict()
        if self.fai is not None:
            for line in self.fai:
                name, length = line.rstrip().split()[:2]
                current_chrom = Chrom(name, length=int(length))
                self.session.add(current_chrom)
                sequences[current_chrom.name] = current_chrom.id
            self.session.commit()
            for query in self.session.query(Chrom):
                sequences[query.name] = query.id

        objects = []

        for row in self.BED12:
            if row.header is True:
                continue
            if row.chrom in sequences:
                current_chrom = sequences[row.chrom]
            else:
                current_chrom = Chrom(row.chrom)
                self.session.add(current_chrom)
                self.session.commit()
                sequences[current_chrom.name] = current_chrom.id
                current_chrom = current_chrom.id

            current_junction = Junction(row, current_chrom)
            objects.append(current_junction)
            if len(objects) >= self.maxobjects:
                self.session.bulk_save_objects(objects)
                objects = []
        self.session.bulk_save_objects(objects)
        self.session.commit()

    def __call__(self):
        """
        Alias for serialize()
        """

        self.serialize()
