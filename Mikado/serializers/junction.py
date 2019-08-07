# coding: utf-8

"""
This module is necessary to serialise the junction information
provided in BED12 format by programs such as TopHat or portcullis.
Please convert any other (custom) format into BED12
before loading.
"""

import io
import os
from sqlalchemy import Column
from sqlalchemy import String
from sqlalchemy import Integer
from sqlalchemy import ForeignKey
from sqlalchemy import CHAR
from sqlalchemy import Index
from sqlalchemy import Float
from sqlalchemy import select
from sqlalchemy.orm import relationship, column_property
from sqlalchemy.orm.session import Session
from sqlalchemy.ext.hybrid import hybrid_method
from ..utilities.dbutils import DBBASE, Inspector, connect
from ..parsers import bed12
from ..utilities.log_utils import check_logger, create_default_logger
import pyfaidx


# pylint: disable=too-few-public-methods
class Chrom(DBBASE):
    """
    Simple serialization class for chromosomes.
    """

    __tablename__ = "chrom"
    __table_args__ = {"extend_existing": True}

    chrom_id = Column(Integer, primary_key=True)
    name = Column(String(200), unique=True)
    length = Column(Integer, nullable=True)

    def __init__(self, name, length=None):
        self.name = name
        if length is not None:
            assert isinstance(length, int)
        self.length = length


# A serialisation class must have a ton of attributes ...
# pylint: disable=too-many-instance-attributes
class Junction(DBBASE):
    """
    Class that describes the junction table in the database.
    """

    __tablename__ = "junctions"

    # pylint: disable=invalid-name
    id = Column(Integer, primary_key=True)
    # pylint: enable=invalid-name
    chrom_id = Column(Integer, ForeignKey(Chrom.chrom_id), unique=False)
    start = Column(Integer, nullable=False)
    end = Column(Integer, nullable=False)
    name = Column(String(200))
    strand = Column(CHAR)
    junction_start = Column(Integer, nullable=False)
    junction_end = Column(Integer, nullable=False)
    score = Column(Float)
    __table_args__ = (Index(
        "junction_index", "chrom_id", "junction_start", "junction_end"), {"extend_existing": True})

    chrom_object = relationship(Chrom, uselist=False)
    chrom = column_property(select([Chrom.name]).where(chrom_id == Chrom.chrom_id))

    def __init__(self, bed12_object, chrom_id):
        """
        Serialization initializer.
        :param bed12_object: a BED12-like class instance.
        :type bed12_object: Mikado.parsers.bed12.BED12

        :param chrom_id: the numerical foreign key for the Chrom table.
        :type chrom_id: int
        """

        if not isinstance(bed12_object, bed12.BED12):
            raise TypeError("Invalid data type!")
        self.chrom_id = chrom_id
        self.start = bed12_object.start
        self.end = bed12_object.end
        self.junction_start = bed12_object.thick_start
        self.junction_end = bed12_object.thick_end
        self.name = bed12_object.name
        self.strand = bed12_object.strand
        self.score = bed12_object.score

    def __str__(self):
        return "{chrom}\t{start}\t{end}".format(
            chrom=self.chrom,
            start=self.start,
            end=self.end
        )

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
        # print(self.chrom_object.chrom_id)

        return ((self.chrom == chrom) and (self.start == start) and
                (self.end == end) and (self.strand == strand))
# pylint: enable=too-many-instance-attributes


class JunctionSerializer:

    """
    This class is used to serialize a junction BED12 file into an SQL database.
    """

    def __init__(self, handle,
                 json_conf=None,
                 logger=None):

        """
        :param handle: the file to be serialized.
        :type handle: str | io.IOBase | io.TextIOWrapper

        :param json_conf: Optional configuration dictionary with db connection parameters.
        :type json_conf: dict | None
        """

        self.bed12_parser = None
        self.fai = None

        if logger is not None:
            logger = check_logger(logger)
            self.logger = logger
        else:
            self.logger = create_default_logger("junction")

        if handle is None:
            self.logger.warning("No input file specified. Exiting.")
            self.close()
            return

        self.bed12_parser = bed12.Bed12Parser(handle)
        self.engine = connect(json_conf, logger=logger)

        session = Session(bind=self.engine, autocommit=False, autoflush=False, expire_on_commit=False)
        inspector = Inspector.from_engine(self.engine)
        if Junction.__tablename__ not in inspector.get_table_names():
            DBBASE.metadata.create_all(self.engine)  # @UndefinedVariable

        self.session = session
        if json_conf is not None:
            self.maxobjects = json_conf["serialise"]["max_objects"]
        else:
            self.maxobjects = 10000

        if "genome_fai" not in json_conf["reference"] or not json_conf["reference"]["genome_fai"]:
            _ = pyfaidx.Fasta(json_conf["reference"]["genome"])
            self.fai = _.faidx.indexname
            _.close()
        else:
            self.fai = json_conf["reference"]["genome_fai"]

        if isinstance(self.fai, str):
            assert os.path.exists(self.fai), self.fai
            # noinspection PyTypeChecker
            self.fai = open(self.fai)
        else:
            if self.fai is not None:
                assert isinstance(self.fai, io.TextIOWrapper)

    def serialize(self):
        """
        Workhorse of the class. It parses the input file and loads it into the database.
        """

        sequences = dict()
        if self.bed12_parser is None:
            self.logger.warning("No input file specified. Exiting.")
            return

        if self.fai is not None:
            self.session.begin(subtransactions=True)
            for line in self.fai:
                name, length = line.rstrip().split()[:2]
                if self.session.query(Chrom).filter(Chrom.name == name).all():
                    continue
                try:
                    current_chrom = Chrom(name, length=int(length))
                except ValueError:
                    raise ValueError(line)
                self.session.add(current_chrom)
                sequences[current_chrom.name] = current_chrom.chrom_id
            self.session.commit()
            for query in self.session.query(Chrom):
                sequences[query.name] = query.chrom_id
            self.fai.close()

        objects = []

        for row in self.bed12_parser:
            if row.header is True:
                continue
            if row.chrom in sequences:
                current_chrom = sequences[row.chrom]
            else:
                current_chrom = Chrom(row.chrom)
                self.session.begin(subtransactions=True)
                self.session.add(current_chrom)
                self.session.commit()
                sequences[current_chrom.name] = current_chrom.chrom_id
                current_chrom = current_chrom.chrom_id

            current_junction = Junction(row, current_chrom)
            objects.append(current_junction)
            if len(objects) >= self.maxobjects:
                self.logger.debug("Serializing %d objects", len(objects))
                self.session.begin(subtransactions=True)
                self.session.bulk_save_objects(objects)
                self.session.commit()
                objects = []

        self.session.bulk_save_objects(objects)
        self.session.commit()
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

        if self.bed12_parser is not None:
            self.bed12_parser.close()
        if self.fai is not None:
            self.fai.close()

# pylint: enable=too-few-public-methods
