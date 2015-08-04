# coding: utf-8

"""This module is used to serialise BLAST objects into a database.
It consists of various different classes:

- Query
- Target
- Hit
- HSP

The module also contains helper functions such as mean().
"""

import os
import sqlalchemy
import gzip
import subprocess
from Bio import SeqIO
import operator
import sqlalchemy.exc
from sqlalchemy.ext.hybrid import hybrid_property  # hybrid_method
from sqlalchemy import Column, String, Integer, Float, ForeignKey, Index
from sqlalchemy.sql.schema import PrimaryKeyConstraint, UniqueConstraint
from sqlalchemy.orm import relationship, backref
from Bio.Blast.NCBIXML import parse as xparser
import io
from sqlalchemy import create_engine
from sqlalchemy.orm.session import sessionmaker
from mikado_lib.serializers.dbutils import dbBase
import logging


def mean(l: list):
    """
    :param l: a list of integers
    :type l: list(int)

    Very simple function to calculate the mean of a numeric array.
    """

    if len(l) == 0:
        raise ZeroDivisionError
    return sum(l) / len(l)


def merge(intervals: [(int, int)]):
    """
    :param intervals: a list of integer duplexes
    :type intervals: list

    This function is used to merge together intervals, which have to be supplied as a list
    of duplexes - (start,stop). The function will then merge together overlapping tuples and
    return a list of non-overlapping tuples.
    If the list is composed by only one element, the function returns immediately.
    """

    # Assume tuple of the form (start,end)
    # And return 0- and 1-length intervals
    new_intervals = []
    for interval in intervals:
        new_intervals.append(tuple(sorted(interval)))

    intervals = new_intervals[:]
    if len(intervals) < 2:
        return intervals

    # Sort according to start, end
    intervals = sorted(intervals, key=operator.itemgetter(0, 1))
    final_list = [intervals[0]]

    for start, end in intervals[1:]:
        if start > final_list[-1][1]:
            final_list.append(tuple([start, end]))
        elif end > final_list[-1][1]:
            final_list[-1] = tuple([final_list[-1][0], end])
    return final_list


class Query(dbBase):
    """
    Very simple serialization class for Query objects.

    :param id: integer key
    :type id: int

    :param name: name of the queries
    :type name: str

    :param length: length of the queries
    :type length: int
    """

    __tablename__ = "query"
    query_id = Column(Integer, primary_key=True)
    query_name = Column(String(200), unique=True, index=True)
    query_length = Column(Integer, nullable=True)  # This so we can load data also from the orf class

    def __init__(self, name, length):
        self.query_name = name
        self.query_length = length


class Target(dbBase):

    """
    Very simple serialization class for Target objects.
    """

    __tablename__ = "target"

    target_id = Column(Integer, primary_key=True)
    target_name = Column(String(200), unique=True, index=True)
    target_length = Column(Integer)

    def __init__(self, target_name, target_length):
        """
        Constructor method.
        :param target_name: name of the targets
        :type target_name: str

        :param target_length: length of the targets
        :type target_length: int
        """

        self.target_name = target_name
        self.target_length = target_length


class Hsp(dbBase):
    """
    This class serializes and stores into the DB the various HSPs.
    It is directly connected to the Hit table, through the "hit_id" 
    reference key.
    The Hit reference can be accessed through the hit_object attribute;
    back-reference (Hit to Hsps) is given by the "hsps" attribute.
    
    Keys:

    :param hit_id: Reference key for the Hit table
    :type hit_id: int

    :param counter: It indicates the progressive number of the HSP for the hit
    :type counter: int

    :param query_hsp_start: Start position on the query
    :type query_hsp_start; int

    :param query_hsp_end: End position on the query
    :type query_hsp_end: int

    :param target_hsp_start: Start position on the target
    :type target_hsp_start: int

    :param target_hsp_end: End position on the target
    :type target_hsp_end: int

    :param hsp_evalue: Evalue of the HSP
    :type hsp_evalue: float

    :param hsp_bits: Bit-score of the HSP
    :type hsp_bits: float

    :param hsp_identity: Identity (in %) of the alignment
    :type  hsp_identity: float

    :param hsp_length: Length of the HSP
    :type hsp_length: int
    
    An HSP row has the following constraints:
    - Counter,hit_id must be unique (and are primary keys)
    - The combination ("Hit_id","query_hsp_start","query_hsp_end", "target_hsp_start", "target_hsp_end") must be unique

    Moreover, the following properties are also present:

    :param query_object: The referenced Query
    :type query_object: Query

    :param target_object: The reference Target
    :type target_object: Target

    """

    __tablename__ = "hsp"
    counter = Column(Integer)  # Indicates the number of the HSP inside the hit
    query_id = Column(Integer, ForeignKey(Query.query_id), unique=False)
    target_id = Column(Integer, ForeignKey(Target.target_id), unique=False)
    pk_constraint = PrimaryKeyConstraint("counter", "query_id", "target_id", name="hsp_constraint")
    query_index = Index("hsp_query_idx", "query_id", unique=False)
    target_index = Index("hsp_target_idx", "target_id", unique=False)
    query_hsp_start = Column(Integer)
    query_hsp_end = Column(Integer)
    target_hsp_start = Column(Integer)
    target_hsp_end = Column(Integer)
    uni_constraint = UniqueConstraint("query_id", "target_id", "query_hsp_start", "query_hsp_end", "target_hsp_start",
                                      "target_hsp_end")
    hsp_evalue = Column(Float)
    hsp_bits = Column(Float)
    hsp_identity = Column(Float)
    hsp_length = Column(Integer)

    query_object = relationship(Query, uselist=False)
    target_object = relationship(Target, uselist=False)

    __table_args__ = (pk_constraint,)

    def __init__(self, hsp, counter, query_id, target_id):

        """

        :param hsp: an hsp object from the serialized XML

        :param counter: integer which indicates the HSP position in the HSP list for the hit
        :type counter: int

        :param query_id: Foreign key for the Query table
        :type query_id: int

        :param target_id: Foreign key for the Target table
        :type target_id: int
        """

        self.counter = counter
        self.query_hsp_start = hsp.query_start
        self.query_hsp_end = hsp.query_end
        self.target_hsp_start = hsp.sbjct_start
        self.target_hsp_end = hsp.sbjct_end
        self.hsp_identity = float(hsp.identities) / hsp.align_length * 100
        self.hsp_length = hsp.align_length
        self.hsp_bits = hsp.bits
        self.hsp_evalue = hsp.expect
        self.query_id = query_id
        self.target_id = target_id

    def __str__(self):
        """Simple printing function."""
        line = [self.query, self.target, self.query_hsp_start, self.query_hsp_end, self.target_hsp_start,
                self.target_hsp_end, self.hsp_evalue]
        return "\t".join([str(x) for x in line])

    # @profile
    def as_dict(self):
        """Method to return a dict representation of the object. Necessary for storing.
        This method returns a dictionary *without any attribute that requires joined data*.
        It is meant to be used only by the method as_dict of the Hit class."""

        keys = [
            "query_hsp_start",
            "query_hsp_end",
            "target_hsp_start",
            "target_hsp_end",
            "hsp_evalue",
            "hsp_bits"
        ]

        state = dict().fromkeys(keys)

        for key in keys:
            state[key] = getattr(self, key)

        return state

    def as_full_dict(self):
        """Method to return a dict representation of the object.
        This method also checks query name and hit name, so it is slower than as_dict and used
        when it is necessary to retrieve data independently from Hit.
        """

        state = self.as_dict()
        state["query"] = self.query
        state["target"] = self.target
        state["query_hsp_cov"] = self.query_hsp_cov
        state["target_hsp_cov"] = self.target_hsp_cov

        return state

    @hybrid_property
    def query(self):
        """Returns the name of the query sequence, through a nested SQL query."""
        return self.query_object.query_name

    @hybrid_property
    def target(self):
        """Returns the name of the target sequence, through a nested SQL query."""
        return self.target_object.target_name

    @hybrid_property
    def query_hsp_cov(self):
        """This property returns the percentage of the query which is covered by the HSP."""
        return (self.query_hsp_end - self.query_hsp_start + 1) / self.query_object.query_length

    @hybrid_property
    def target_hsp_cov(self):
        """This property returns the percentage of the target which is covered by the HSP."""
        return (self.target_hsp_end - self.target_hsp_start + 1) / self.target_object.target_length


class Hit(dbBase):
    """This class is used to serialise and store in a DB a BLAST hit.
    Stored attributes:
    
    - id                Indexing key
    - query_id            Foreign ID key for the query table
    - target_id            Foreign ID key for the target table
    - qt_constrating    
    
    """

    __tablename__ = "hit"
    query_id = Column(Integer, ForeignKey(Query.query_id), unique=False)
    target_id = Column(Integer, ForeignKey(Target.target_id), unique=False)
    qt_constraint = PrimaryKeyConstraint("query_id", "target_id", name="hit_id")
    query_index = Index("hit_query_idx", "query_id", unique=False)
    target_index = Index("hit_target_idx", "query_id", unique=False)
    evalue = Column(Float)
    bits = Column(Float)
    global_identity = Column(Float)
    query_start = Column(Integer)
    query_end = Column(Integer)
    target_start = Column(Integer)
    target_end = Column(Integer)
    hit_number = Column(Integer)
    query_multiplier = Column(Float)  # Probably I should move this to a separate table!
    target_multiplier = Column(Float)
    query_aligned_length = Column(Integer)
    target_aligned_length = Column(Integer)

    query_object = relationship(Query, uselist=False, lazy="immediate", backref=backref("hits"))
    target_object = relationship(Target, uselist=False, lazy="immediate", backref=backref("hits"))

    hsps = relationship(Hsp, uselist=True, lazy="subquery", backref=backref("hit_object", uselist=False),
                        foreign_keys=[query_id, target_id],
                        primaryjoin="and_(Hit.query_id==Hsp.query_id, Hit.target_id==Hsp.target_id)")

    __table_args__ = (qt_constraint,)

    def __init__(self, query_id, target_id, alignment, evalue, bits, hit_number=1, query_multiplier=1,
                 target_multiplier=1):
        """This function takes as input the id of a target, the id of the query, and a hit-object from the XML.
        The multiplier keyword is used to calculate the ratio between the query and the target.

        :param query_id: reference key for the Query table
        :type query_id: int

        :param target_id: reference key for the Target table
        :type target_id: int

        :param alignment: The BLAST Hit object

        :param evalue: Evalue of the hit (recovered from the "description" object)
        :type evalue: float

        :param bits: BitScore of the hit (recovered from the "description" object)
        :type bits: float

        :param hit_number: progressive index that indicates the priority of the hit in the database
        :type hit_number: int

        :param query_multiplier: either 1 or 3, it depends on the type of BLAST
        :type query_multiplier: int

        :param target_multiplier: either 1 or 3, it depends on the type of BLAST
        :type target_multiplier: int
        """

        self.query_id = query_id
        self.target_id = target_id
        self.query_multiplier = query_multiplier
        self.target_multiplier = target_multiplier
        self.hit_number = hit_number
        self.evalue = evalue
        self.bits = bits

        self.global_identity = mean([hsp.identities / hsp.align_length * 100 for hsp in alignment.hsps])
        q_intervals = [tuple([hsp.query_start, hsp.query_end]) for hsp in alignment.hsps]
        q_merged_intervals = sorted(merge(q_intervals), key=operator.itemgetter(0, 1))
        q_aligned = sum([tup[1] - tup[0] + 1 for tup in q_merged_intervals])
        self.query_aligned_length = q_aligned
        self.query_start = q_merged_intervals[0][0]
        self.query_end = q_merged_intervals[-1][1]

        t_intervals = [tuple([hsp.sbjct_start, hsp.sbjct_end]) for hsp in alignment.hsps]
        t_merged_intervals = sorted(merge(t_intervals), key=operator.itemgetter(0, 1))
        t_aligned = sum([tup[1] - tup[0] + 1 for tup in t_merged_intervals])
        self.target_aligned_length = t_aligned
        self.target_start = t_merged_intervals[0][0]
        self.target_end = t_merged_intervals[-1][1]

    def __str__(self):
        line = [self.query, self.target, self.evalue, self.bits, self.query_start, self.query_end, self.target_start,
                self.target_end, self.query_len, self.target_len]

        return "\t".join(str(x) for x in line)

    # @profile
    def as_dict(self):
        """Method to return a dict representation of the object.
        Necessary for storing."""

        keys = [
            "evalue",
            "bits",
            "global_identity",
            "query_start",
            "query_end",
            "target_start",
            "target_end",
            "hit_number",
            "query_multiplier",
            "target_multiplier",
            "query_aligned_length",
            "target_aligned_length",
        ]

        special_keys = ["query", "target", "hsps", "query_len",
                        "target_len", "query_hit_ratio", "hit_query_ratio"]  # Keys we want to set directly

        state = dict().fromkeys(keys + special_keys)

        for key in keys:
            state[key] = getattr(self, key)

        state["query"] = self.query
        state["target"] = self.target
        state["query_len"] = self.query_len
        state["target_len"] = self.target_len
        state["query_hit_ratio"] = self.query_hit_ratio
        state["hit_query_ratio"] = self.hit_query_ratio

        state["hsps"] = []
        for hsp in self.hsps:
            h = hsp.as_dict()
            h["query_hsp_cov"] = (h["query_hsp_end"] - h["query_hsp_start"] + 1) / (state["query_len"])
            h["target_hsp_cov"] = (h["target_hsp_end"] - h["target_hsp_start"] + 1) / (state["target_len"])
            state["hsps"].append(hsp.as_dict())

        return state

    @hybrid_property
    def query_len(self):
        """
        This method retrieves the length of the query from the entry in the Query table.
        :rtype int
        """
        return self.query_object.query_length

    @hybrid_property
    def query(self):
        """
        This method retrieves the name of the query from the entry in the Query table.
        :rtype str
        """
        return self.query_object.query_name

    @hybrid_property
    def target(self):
        """
        This method retrieves the name of the target from the entry in the Target table.
        :rtype str
        """

        return self.target_object.target_name

    @hybrid_property
    def target_len(self):
        """
        This method retrieves the length of the target from the entry in the Target table.
        :rtype int
        """

        return self.target_object.target_length

    @hybrid_property
    def query_hit_ratio(self):
        """
        This property returns the quotient (Query Length)/(Target Length)
        """

        return self.query_len * self.query_multiplier / (self.target_len * self.target_multiplier)

    @hybrid_property
    def hit_query_ratio(self):
        """
        This property returns the quotient (Target Length)/(Query Length)
        """

        return self.target_len * self.target_multiplier / (self.query_len * self.query_multiplier)


class XmlSerializer:
    """This class has the role of taking in input a blast XML file and (partially) serialise it into
    a database. We are using SQLalchemy, so the database type could be any of SQLite, MySQL, PSQL, etc."""

    def __init__(self, db, xml, max_target_seqs=float("Inf"),
                 target_seqs=None,
                 query_seqs=None,
                 keep_definition=False, maxobjects=10000,
                 json_conf=None
                 ):
        """Initializing method. Arguments:

        :param db: the name of the SQLite database
        :type db: str

        :param xml: The XML to parse.

        Optional arguments:

        :param target_seqs: either a BioPython index of a FASTA file or the file itself.

        :param query_seqs: either a BioPython index of a FASTA file or the file itself.

        :param keep_definition: flag. If set, the "definition" field in the XML
        instead of "id" will be used for serializing.
        :type keep_definition: bool

        :param maxobjects: maximum number of objects to keep in memory before bulk loading. Default: 10^4
        :type maxobjects: int

        :param json_conf: a configuration dictionary.
        :type json_conf: dict


        """

        self.logger = logging.getLogger("main")
        self.logger.setLevel(logging.INFO)
        self.handler = logging.StreamHandler()
        self.formatter = logging.Formatter("{asctime} - {levelname} - {message}", style='{')
        self.handler.setFormatter(self.formatter)
        self.logger.addHandler(self.handler)

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
        dbBase.metadata.create_all(self.engine)  # @UndefinedVariable
        self.session = session()
        self.logger.info("Created the session")
        if type(xml) is not xparser:
            if type(xml) is str:
                if not os.path.exists(xml):
                    raise OSError("Invalid file address.")
                if xml.endswith("xml"):
                    xml = open(xml)
                elif xml.endswith(".xml.gz"):
                    xml = gzip.open(xml, mode="rt")  # The mode must be textual, otherwise xparser dies
                elif xml.endswith(".asn.gz"):
                    zcat = subprocess.Popen(["zcat", xml], shell=False,
                                            stdout=subprocess.PIPE)  # I cannot seem to make it work with gzip.open
                    blast_formatter = subprocess.Popen(['blast_formatter', '-outfmt', '5',
                                                        '-archive', '-'], shell=False,
                                                       stdin=zcat.stdout,
                                                       stdout=subprocess.PIPE)
                    xml = io.TextIOWrapper(blast_formatter.stdout, encoding="UTF-8")
            assert type(xml) in (gzip.GzipFile, io.TextIOWrapper), type(xml)
            self.xml = xml
            self.xml_parser = xparser(xml)
        else:
            self.xml_parser = xml  # This is the BLAST object we will serialize
        # Runtime arguments
        self.keep_definition = keep_definition

        if type(query_seqs) is str:
            assert os.path.exists(query_seqs)
            self.query_seqs = SeqIO.index(query_seqs, "fasta")
        elif query_seqs is None:
            self.query_seqs = None
        else:
            assert "SeqIO.index" in repr(query_seqs)
            self.query_seqs = query_seqs

        if type(target_seqs) is str:
            assert os.path.exists(target_seqs)
            self.target_seqs = SeqIO.index(target_seqs, "fasta")
        elif target_seqs is None:
            self.target_seqs = None
        else:
            assert "SeqIO.index" in repr(target_seqs)
            self.target_seqs = target_seqs

        self.max_target_seqs = max_target_seqs
        self.maxobjects = maxobjects

    def serialize(self):

        """Method to serialize the BLAST XML file into a database provided with the __init__ method """

        q_mult = 1  # Query multiplier: 3 for BLASTX (nucleotide vs. protein), 1 otherwise
        h_mult = 1  # Target multiplier: 3 for TBLASTN (protein vs. nucleotide), 1 otherwise

        objects = []

        targets = dict()
        queries = dict()
        self.logger.info("Loading previous IDs")
        for query in self.session.query(Query):
            queries[query.query_name] = query.query_id
        for query in self.session.query(Target):
            targets[query.target_name] = query.target_id
        self.logger.info("Loaded previous IDs")

        self.logger.info("Started the serialisation")
        if self.target_seqs is not None:
            self.logger.info("Started to serialise the targets")
            for record in self.target_seqs:
                if record in targets:
                    continue
                objects.append(Target(record, len(self.target_seqs[record])))
                if len(objects) >= self.maxobjects:
                    self.logger.info("Loading {0} objects into the \"target\" table".format(len(objects)))
                    self.engine.execute(Target.__table__.insert(),
                                        [{"target_name": obj.target_name,
                                          "target_length": obj.target_length} for obj in objects]
                                        )
                    #                         self.session.bulk_save_objects(objects, return_defaults=False)
                    self.logger.info("Loaded {0} objects into the \"target\" table".format(len(objects)))
                    objects = []
            self.logger.info("Loading {0} objects into the \"target\" table".format(len(objects)))
            self.engine.execute(Target.__table__.insert(),
                                [{"target_name": obj.target_name,
                                  "target_length": obj.target_length} for obj in objects]
                                )
            #             self.session.bulk_save_objects(objects)
            self.logger.info("Loaded {0} objects into the \"target\" table".format(len(objects)))
            objects = []
            self.logger.info("Loaded targets")
            self.session.commit()

        if self.query_seqs is not None:
            self.logger.info("Started to serialise the queries")
            for record in self.query_seqs:
                if record in queries:
                    continue
                objects.append(Query(record, len(self.query_seqs[record])))
                if len(objects) >= self.maxobjects:
                    self.logger.info("Loading {0} objects into the \"query\" table".format(len(objects)))
                    self.engine.execute(Query.__table__.insert(),
                                        [{"query_name": obj.query_name,
                                          "query_length": obj.query_length} for obj in objects]
                                        )
                    #                     self.session.bulk_save_objects(objects, return_defaults=False)
                    self.logger.info("Loaded {0} objects into the \"query\" table".format(len(objects)))
                    objects = []
            self.logger.info("Loading {0} objects into the \"query\" table".format(len(objects)))
            self.engine.execute(Query.__table__.insert(),
                                [{"query_name": obj.query_name, "query_length": obj.query_length} for obj in objects]
                                )
            #             self.session.bulk_save_objects(objects)
            self.session.commit()
            self.logger.info("Loaded {0} objects into the \"query\" table".format(len(objects)))
            objects = []
            self.logger.info("Queries serialised")

        self.logger.info("Loading all IDs")
        for query in self.session.query(Query):
            queries[query.query_name] = (query.query_id, query.query_length is not None)
        for query in self.session.query(Target):
            targets[query.target_name] = (query.target_id, query.target_length is not None)
        self.logger.info("Loaded all IDs")

        # Memorize the mapping ... it is just faster

        query_counter = 0
        for record in self.xml_parser:
            if record.application == "BLASTN":
                q_mult = 1
                h_mult = 1
            elif record.application == "BLASTX":
                q_mult = 3
                h_mult = 1
            elif record.application == "TBLASTN":
                q_mult = 1
                h_mult = 3
            if len(record.descriptions) == 0:
                continue
            query_counter += 1
            if self.keep_definition is True:
                name = record.query.split()[0]
            else:
                name = record.query_id
            self.logger.debug("Started with {0}".format(name))

            if name in queries:
                current_query = queries[name][0]
                if queries[name][1] is False:
                    self.session.query(Query).filter(Query.query_name == name).update(
                        {"query_length": record.query_length})
                    self.session.commit()
            else:
                self.logger.warn("Adding {0} to the db".format(name))
                current_query = Query(name, record.query_length)
                self.session.add(current_query)
                self.session.commit()
                queries[name] = (current_query.query_id, True)
                current_query = current_query.query_id

            for ccc, alignment in filter(lambda x: x[0] <= self.max_target_seqs, enumerate(record.alignments)):

                self.logger.debug("Started the hit {0}-{1}".format(name, record.alignments[ccc].accession))
                evalue = record.descriptions[ccc].e
                bits = record.descriptions[ccc].bits
                alignment = record.alignments[ccc]
                hit_num = ccc + 1
                if record.alignments[ccc].accession in targets:
                    current_target = targets[record.alignments[ccc].accession][0]
                    if targets[record.alignments[ccc].accession][1] is False:
                        self.session.query(Target).filter(Target.target_name == record.alignments[ccc].accession).\
                            update({"target_length": record.query_length})

                else:
                    current_target = Target(record.alignments[ccc].accession, record.alignments[ccc].length)
                    self.session.add(current_target)
                    try:
                        self.session.commit()
                        assert type(current_target.target_id) is int
                        targets[record.alignments[ccc].accession] = (current_target.target_id, True)
                        current_target = current_target.target_id
                    except sqlalchemy.exc.IntegrityError:
                        self.session.rollback()
                        continue

                current_hit = Hit(current_query, current_target, alignment, evalue, bits,
                                  hit_number=hit_num,
                                  query_multiplier=q_mult,
                                  target_multiplier=h_mult)
                objects.append(current_hit)

                for counter, hsp in enumerate(alignment.hsps):
                    current_hsp = Hsp(hsp, counter, current_query, current_target)
                    objects.append(current_hsp)

            if len(objects) >= self.maxobjects:
                self.logger.info("Loading {0} objects into the hit, hsp tables; done {1} queries".format(len(objects),
                                                                                                         query_counter))
                self.load_into_db(objects)
                self.logger.info(
                    "Loaded {0} objects into the hit, hsp tables; done {1} queries".format(len(objects), query_counter))
                objects = []

        self.logger.info("Loading {0} objects into the hit, hsp tables".format(len(objects)))
        self.load_into_db(objects)

        self.logger.info("Loaded {0} objects into the hit, hsp tables".format(len(objects)))
        self.session.commit()
        self.logger.info("Finished loading blast hits")

    def __call__(self):
        """
        Alias for serialize
        """
        self.serialize()

    def load_into_db(self, objects):
        """
        :param objects: Objects to be loaded into the database
        :type objects: list

        Method to perform the bulk loading of objects into the SQL database.

        """

        try:
            self.session.bulk_save_objects(objects, return_defaults=False)
            self.session.commit()
        except Exception as err:
            self.logger.error('Database corrupted'.format(err))
            self.logger.error(err)
            self.logger.error('Dropping and reloading')
            self.session.rollback()
            self.session.query(Hsp).delete()
            self.session.query(Hit).delete()
            self.session.bulk_save_objects(objects, return_defaults=False)
            self.session.commit()
