# coding: utf-8

"""This module is used to serialise BLAST objects into a database.
It consists of various different classes:

- Query
- Target
- Hit
- HSP

The module also contains helper functions such as mean().
"""

import collections
import os
import sqlalchemy
import functools
from Bio import SeqIO
import operator
import sqlalchemy.exc
from sqlalchemy.ext.hybrid import hybrid_property  # hybrid_method
from sqlalchemy import Column, String, Integer, Float, ForeignKey, Index
from sqlalchemy.sql.schema import PrimaryKeyConstraint, UniqueConstraint
from sqlalchemy.orm import relationship, backref
from Bio.Blast.NCBIXML import parse as xparser
from sqlalchemy.orm.session import sessionmaker
from mikado_lib.serializers.dbutils import DBBASE
from mikado_lib.serializers.dbutils import connect
from mikado_lib.parsers.blast_utils import XMLMerger, create_opener, merge
from mikado_lib.log_utils import create_null_logger, check_logger
# pylint: disable=no-name-in-module
from numpy import mean
# pylint: enable=no-name-in-module


# These two classes are OK like this, they do not need more public methods!
# pylint: disable=too-few-public-methods
class Query(DBBASE):
    """
    Very simple serialization class for Query objects.

    :return id: integer key
    :rtype id: int

    :return name: name of the queries
    :rtype name: str

    :return length: length of the queries
    :rtype length: int
    """

    __tablename__ = "query"
    query_id = Column(Integer, primary_key=True)
    query_name = Column(String(200), unique=True, index=True)
    # This so we can load data also from the orf class
    query_length = Column(Integer, nullable=True)

    named_tup = collections.namedtuple("Query",
                                       ["query_id", "query_name", "query_length"])

    def __init__(self, name, length):
        self.query_name = name
        self.query_length = length

    def as_tuple(self):
        """Quick function to convert the SQLalchemy object into
        a named tuple with the same fields"""
        return self.named_tup(self.query_id, self.query_name, self.query_length)


class Target(DBBASE):

    """
    Very simple serialization class for Target objects.
    """

    __tablename__ = "target"

    target_id = Column(Integer, primary_key=True)
    target_name = Column(String(200), unique=True, index=True)
    target_length = Column(Integer)
    named_tup = collections.namedtuple("Target",
                                       ["target_id", "target_name", "target_length"])

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

    def as_tuple(self):
        """Quick function to convert the SQLalchemy object
        into a named tuple with the same fields"""
        return self.named_tup(
            self.target_id, self.target_name, self.target_length)
# pylint: enable=too-few-public-methods


class Hsp(DBBASE):
    """
    This class serializes and stores into the DB the various HSPs.
    It is directly connected to the Hit table, through the "hit_id"
    reference key.
    The Hit reference can be accessed through the hit_object attribute;
    back-reference (Hit to Hsps) is given by the "hsps" attribute.

    Keys:

    :return hit_id: Reference key for the Hit table
    :rtype hit_id: int

    :return counter: It indicates the progressive number of the HSP for the hit
    :rtype counter: int

    :return query_hsp_start: Start position on the query
    :rtype query_hsp_start; int

    :return query_hsp_end: End position on the query
    :rtype query_hsp_end: int

    :return target_hsp_start: Start position on the target
    :rtype target_hsp_start: int

    :return target_hsp_end: End position on the target
    :rtype target_hsp_end: int

    :return hsp_evalue: Evalue of the HSP
    :rtype hsp_evalue: float

    :return hsp_bits: Bit-score of the HSP
    :rtype hsp_bits: float

    :return hsp_identity: Identity (in %) of the alignment
    :rtype  hsp_identity: float

    :return hsp_length: Length of the HSP
    :rtype hsp_length: int

    An HSP row has the following constraints:
    - Counter,hit_id must be unique (and are primary keys)
    - The combination ("Hit_id","query_hsp_start","query_hsp_end",
    "target_hsp_start", "target_hsp_end") must be unique

    Moreover, the following properties are also present:

    :return query_object: The referenced Query
    :rtype query_object: Query

    :return target_object: The reference Target
    :rtype target_object: Target

    """

    __tablename__ = "hsp"
    counter = Column(Integer)  # Indicates the number of the HSP inside the hit
    query_id = Column(Integer, ForeignKey(Query.query_id), unique=False)
    target_id = Column(Integer, ForeignKey(Target.target_id), unique=False)
    pk_constraint = PrimaryKeyConstraint("counter",
                                         "query_id",
                                         "target_id",
                                         name="hsp_constraint")
    query_index = Index("hsp_query_idx", "query_id", unique=False)
    target_index = Index("hsp_target_idx", "target_id", unique=False)
    combined_index = Index("hsp_combined_idx",
                           "query_id",
                           "target_id", unique=False)
    full_index = Index("hsp_full_idx",
                       "counter", "query_id", "target_id", unique=True)
    query_hsp_start = Column(Integer)
    query_hsp_end = Column(Integer)
    target_hsp_start = Column(Integer)
    target_hsp_end = Column(Integer)
    uni_constraint = UniqueConstraint("query_id", "target_id",
                                      "query_hsp_start", "query_hsp_end",
                                      "target_hsp_start", "target_hsp_end")
    hsp_evalue = Column(Float)
    hsp_bits = Column(Float)
    hsp_identity = Column(Float)
    hsp_length = Column(Integer)

    query_object = relationship(Query, uselist=False)
    target_object = relationship(Target, uselist=False)

    __table_args__ = (pk_constraint, query_index, target_index, combined_index)

    def __init__(self, hsp, counter, query_id, target_id):

        """

        :param hsp: an hsp object from the serialized XML

        :param counter: integer which indicates the HSP position in the HSP
         list for the hit
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
        line = [self.query, self.target, self.query_hsp_start,
                self.query_hsp_end, self.target_hsp_start,
                self.target_hsp_end, self.hsp_evalue]
        return "\t".join([str(x) for x in line])

    # @profile

    @classmethod
    def as_dict_static(cls, state_obj):
        """Method to return a dict representation of the object. Necessary for storing.
        This method returns a dictionary *without any attribute that requires joined data*.
        As a static method, it is useful to be used outside of the class.
        :param state_obj: an instance of the HSP class or a
        collections.namedtuple object from a direct query

        :rtype : dict
        """

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
            state[key] = getattr(state_obj, key)
        return state

    def as_dict(self):
        """Method to return a dict representation of the object. Necessary for storing.
        This method returns a dictionary *without any attribute that requires joined data*.
        It is meant to be used only by the method as_dict of the Hit class."""

        return self.as_dict_static(self)

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
        val = (self.query_hsp_end - self.query_hsp_start + 1)
        val /= self.query_object.query_length
        return val

    @hybrid_property
    def target_hsp_cov(self):
        """This property returns the percentage of the target which is covered by the HSP."""
        val = (self.target_hsp_end - self.target_hsp_start + 1)
        val /= self.target_object.target_length
        return val


class Hit(DBBASE):
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
    qt_index = Index("qt_index", "query_id", "target_id", unique=True)
    query_index = Index("hit_query_idx", "query_id", unique=False)
    target_index = Index("hit_target_idx", "target_id", unique=False)
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

    query_object = relationship(Query, uselist=False,
                                lazy="immediate",
                                backref=backref("hits"))
    target_object = relationship(Target,
                                 uselist=False,
                                 lazy="immediate",
                                 backref=backref("hits"))

    join_condition = "and_(Hit.query_id==Hsp.query_id, Hit.target_id==Hsp.target_id)"
    hsps = relationship(Hsp, uselist=True,
                        # lazy="immediate",
                        backref=backref("hit_object", uselist=False),
                        foreign_keys=[query_id, target_id],
                        primaryjoin=join_condition)

    __table_args__ = (qt_constraint, qt_index, query_index, target_index)

    def __init__(self, query_id, target_id, alignment,
                 evalue, bits, hit_number=1, query_multiplier=1,
                 target_multiplier=1):
        """This function takes as input the id of a target, the id of the query,
        and a hit-object from the XML. The multiplier keyword is used to calculate
        the ratio between the query and the target.

        :param query_id: reference key for the Query table
        :type query_id: int

        :param target_id: reference key for the Target table
        :type target_id: int

        :param alignment: The BLAST Hit object

        :param evalue: Evalue of the hit (recovered from the "description" object)
        :type evalue: float

        :param bits: BitScore of the hit (recovered from the "description" object)
        :type bits: float

        :param hit_number: progressive index that indicates the priority
        of the hit in the database
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

        self.global_identity = mean(
            [hsp.identities / hsp.align_length * 100 for hsp in alignment.hsps])
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
        line = [self.query, self.target, self.evalue,
                self.bits, self.query_start, self.query_end,
                self.target_start, self.target_end,
                self.query_len, self.target_len]

        return "\t".join(str(x) for x in line)

    # @profile

    @classmethod
    def as_dict_static(cls, state_obj):
        """Method to return a dict representation of the object.
        Static method to be called from outside the class.

        :param state_obj: a namedtuple or an instance of this class

        :rtype: dict
        """

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

        state = dict().fromkeys(keys)

        for key in keys:
            state[key] = getattr(state_obj, key)

        return state

    @classmethod
    def as_full_dict_static(cls, hit_tuple,
                            hsp_list,
                            query_tuple,
                            target_tuple):
        """
        :param hit_tuple: Hit namedtuple (from direct query to the DB)
        :type hit_tuple: collections.namedtuple
        :param hsp_list: list of hsp dictionaries from Hsp.as_dict_static
        :type hsp_list: list(collections.namedtuple)

        :param query_tuple: Query namedtuple
        :type query_tuple: collections.namedtuple
        :param target_tuple: Target namedtuple
        :type target_tuple: collections.namedtuple
        :rtype: dict
        """

        state = cls.as_dict_static(hit_tuple)
        hsps = [Hsp.as_dict_static(h) for h in hsp_list]

        state["query"] = query_tuple.query_name
        state["target"] = target_tuple.target_name
        state["query_len"] = query_tuple.query_length
        state["target_len"] = target_tuple.target_length
        state["query_hit_ratio"] = (query_tuple.query_length * hit_tuple.query_multiplier) / \
                                   (target_tuple.target_length * hit_tuple.target_multiplier)
        state["hit_query_ratio"] = (target_tuple.target_length * hit_tuple.target_multiplier) / \
                                   (query_tuple.query_length * hit_tuple.query_multiplier)
        state["hsps"] = []
        for hsp in hsps:
            hsp["query_hsp_cov"] = (hsp["query_hsp_end"] - hsp["query_hsp_start"] + 1)
            hsp["query_hsp_cov"] /= (state["query_len"])
            hsp["target_hsp_cov"] = (hsp["target_hsp_end"] - hsp["target_hsp_start"] + 1)
            hsp["target_hsp_cov"] /= (state["target_len"])
            state["hsps"].append(hsp)

        return state

    def as_dict(self):
        """Method to return a dict representation of the object.
        Necessary for storing.

        :rtype: dict
        """

        state = self.as_dict_static(self)

        # Retrieving the values ONCE
        query_object = self.query_object.as_tuple()
        target_object = self.target_object.as_tuple()

        state["query"] = query_object.query_name
        state["target"] = target_object.target_name
        state["query_len"] = query_object.query_length
        state["target_len"] = target_object.target_length
        state["query_hit_ratio"] = state["query_len"] * state["query_multiplier"] /\
            (state["target_len"] * state["target_multiplier"])

        state["hit_query_ratio"] = state["target_len"] * state["target_multiplier"] /\
            (state["query_len"] * state["query_multiplier"])

        state["hsps"] = []
        for hsp in self.hsps:
            dict_hsp = hsp.as_dict()
            dict_hsp["query_hsp_cov"] = (dict_hsp["query_hsp_end"] -
                                         dict_hsp["query_hsp_start"] + 1)
            dict_hsp["query_hsp_cov"] /= state["query_len"]
            dict_hsp["target_hsp_cov"] = (dict_hsp["target_hsp_end"] -
                                          dict_hsp["target_hsp_start"] + 1)
            dict_hsp["target_hsp_cov"] /= state["target_len"]
            state["hsps"].append(dict_hsp)

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
    """This class has the role of taking in input a blast XML file and (partially)
    serialise it into a database. We are using SQLalchemy, so the database type
    could be any of SQLite, MySQL, PSQL, etc."""

    __name__ = "XMLSerializer"
    logger = create_null_logger(__name__)

    def __init__(self, xml,
                 max_target_seqs=float("Inf"),
                 logger=None,
                 target_seqs=None,
                 query_seqs=None,
                 discard_definition=True, maxobjects=10000,
                 json_conf=None):
        """Initializing method. Arguments:

        :param xml: The XML to parse.

        Optional arguments:

        :param target_seqs: either a BioPython index of a FASTA file or the file itself.

        :param query_seqs: either a BioPython index of a FASTA file or the file itself.

        :param discard_definition: flag. If set to True, the "id" field in the XML
        instead of "definition" will be used for serializing.
        :type discard_definition: bool

        :param maxobjects: maximum number of objects to keep in memory
        before bulk loading. Default: 10^4
        :type maxobjects: int

        :param json_conf: a configuration dictionary.
        :type json_conf: dict


        """

        if logger is not None:
            self.logger = check_logger(logger)

        if xml is None:
            self.logger.warning("No BLAST XML provided. Exiting.")
            return

        self.engine = connect(json_conf)

        session = sessionmaker()
        session.configure(bind=self.engine)
        DBBASE.metadata.create_all(self.engine)  # @UndefinedVariable
        self.session = session()
        self.logger.debug("Created the session")
        if isinstance(xml, str):
            self.xml_parser = xparser(create_opener(xml))
        else:
            assert isinstance(xml, (list, set))
            if len(xml) < 1:
                raise ValueError("No input file provided!")
            elif len(xml) == 1:
                self.xml_parser = xparser(create_opener(list(xml)[0]))
            else:
                self.xml_parser = xparser(XMLMerger(xml))  # Merge in memory

        # Runtime arguments
        self.discard_definition = discard_definition
        self.__max_target_seqs = max_target_seqs
        self.maxobjects = maxobjects

        # Load sequences if necessary
        self.__determine_sequences(query_seqs, target_seqs)

    def __determine_sequences(self, query_seqs, target_seqs):

        """Private method to assign the sequence file variables
        if necessary.
        :param query_seqs:
        :param target_seqs:
        :return:
        """

        if isinstance(query_seqs, str):
            assert os.path.exists(query_seqs)
            self.query_seqs = SeqIO.index(query_seqs, "fasta")
        elif query_seqs is None:
            self.query_seqs = None
        else:
            assert "SeqIO.index" in repr(query_seqs)
            self.query_seqs = query_seqs

        if isinstance(target_seqs, str):
            assert os.path.exists(target_seqs)
            self.target_seqs = SeqIO.index(target_seqs, "fasta")
        elif target_seqs is None:
            self.target_seqs = None
        else:
            assert "SeqIO.index" in repr(target_seqs)
            self.target_seqs = target_seqs
        return

    def __serialize_queries(self, queries):

        """Private method used to serialise the queries.
        Additional argument: a set containing all queries already present
        in the database."""

        counter = 0
        self.logger.info("Started to serialise the queries")
        objects = []
        for record in self.query_seqs:
            if record in queries and queries[record][1] is not None:
                continue
            elif record in queries:
                self.session.query(Query).filter(Query.query_name == record).update(
                    {"query_length": record.query_length})
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
                self.engine.execute(Query.__table__.insert(), objects)

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
        self.engine.execute(Query.__table__.insert(), objects)
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
        for record in self.target_seqs:
            if record in targets and targets[record][1] is True:
                continue
            elif record in targets:
                self.session.query(Target).filter(Target.target_name == record).update(
                    {"target_length": record.query_length})
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
                counter += len(objects)
                self.logger.info("Loading %d objects into the \"target\" table",
                                 counter)
                # self.session.bulk_insert_mappings(Target, objects)
                self.session.commit()
                objects = []
                # pylint: disable=no-member
                self.engine.execute(Target.__table__.insert(), objects)
                # self.engine.execute(Target.__table__.insert(),
                #                     [{"target_name": obj.target_name,
                #                       "target_length": obj.target_length} for obj in objects])
                # # pylint: enable=no-member
                # self.logger.info("Loaded %d objects into the \"target\" table",
                #                  len(objects))
                # objects = []
        self.logger.info("Loading %d objects into the \"target\" table, (total %d)",
                         len(objects), counter)
        # pylint: disable=no-member
        self.engine.execute(Target.__table__.insert(), objects)
        # self.engine.execute(Target.__table__.insert(),
        #                     [{"target_name": obj.target_name,
        #                       "target_length": obj.target_length} for obj in objects])
        # self.session.bulk_insert_mappings(Target, objects)
        self.session.commit()
        # pylint: enable=no-member
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

    def serialize(self):

        """Method to serialize the BLAST XML file into a database
        provided with the __init__ method """

        # Load sequences in DB, precache IDs
        queries, targets = self.__serialise_sequences()

        query_counter = 0
        get_query = functools.partial(self.__get_query_for_blast,
                                      **{"queries": queries})

        hits = []
        hsps = []
        hit_counter = 0
        record_counter = 0

        for record in self.xml_parser:
            record_counter += 1
            if record_counter > 0 and record_counter % 10000 == 0:
                self.logger.info("Parsed %d queries", record_counter)
            if len(record.descriptions) == 0:
                continue

            q_mult, h_mult = self.__get_multipliers(record)
            query_counter += 1
            current_query, name = get_query(record)

            for ccc, alignment in enumerate(record.alignments):
                if ccc > self.__max_target_seqs:
                    break
                hit_counter += 1
                if hit_counter > 0 and hit_counter % 10000 == 0:
                    self.logger.info("Serialized %d alignments", hit_counter)

                self.logger.debug("Started the hit %s-%s",
                                  name, record.alignments[ccc].accession)
                evalue = record.descriptions[ccc].e
                bits = record.descriptions[ccc].bits
                hit_num = ccc + 1
                try:
                    current_target, targets = self.__get_target_for_blast(alignment,
                                                                          targets)
                except sqlalchemy.exc.IntegrityError as exc:
                    self.session.rollback()
                    self.logger.exception(exc)
                    continue

                hit_dict_params = dict()
                hit_dict_params["query_multiplier"] = q_mult
                hit_dict_params["target_multiplier"] = h_mult
                hit_dict_params["hit_number"] = hit_num
                hit_dict_params["evalue"] = evalue
                hit_dict_params["bits"] = bits

                # Prepare for bulk load
                hit, hit_hsps = self.__prepare_hit(alignment,
                                               current_query, current_target,
                                               **hit_dict_params)
                hits.append(hit)
                hsps.extend(hit_hsps)

                tot_objects = len(hits) + len(hsps)
                if tot_objects >= self.maxobjects:
                    # Bulk load
                    self.logger.debug("Loading %d BLAST objects into database", tot_objects)
                    self.engine.execute(Hit.__table__.insert(), hits)
                    self.engine.execute(Hsp.__table__.insert(), hsps)
                    # self.session.bulk_insert_mappings(Hit, hits)
                    # self.session.bulk_insert_mappings(Hsp, hsps)
                    self.session.commit()
                    hits = []
                    hsps = []

        self.engine.execute(Hit.__table__.insert(), hits)
        self.engine.execute(Hsp.__table__.insert(), hsps)
        # self.session.bulk_insert_mappings(Hit, hits)
        # self.session.bulk_insert_mappings(Hsp, hsps)
        self.session.commit()

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
            self.session.warn("%s not found among queries, adding to the DB now",
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
            self.session.warn("%s not found among targets, adding to the DB now",
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

    @staticmethod
    def __prepare_hit(hit, query_id, target_id, **kwargs):
        """Prepare the dictionary for fast loading of Hit and Hsp objects"""

        hit_dict = dict()
        hsp_dict_list = []
        hit_dict["global_identity"] = []
        q_intervals = []
        t_intervals = []

        for counter, hsp in enumerate(hit.hsps):
            hsp_dict = dict()
            hsp_dict["counter"] = counter
            hsp_dict["query_hsp_start"] = hsp.query_start
            hsp_dict["query_hsp_end"] = hsp.query_end
            # Prepare the list for later calculation
            q_intervals.append((hsp.query_start, hsp.query_end))

            hsp_dict["target_hsp_start"] = hsp.sbjct_start
            hsp_dict["target_hsp_end"] = hsp.sbjct_end
            # Prepare the list for later calculation
            t_intervals.append((hsp.sbjct_start, hsp.sbjct_end))

            hsp_dict["hsp_identity"] = float(hsp.identities) / hsp.align_length * 100
            # Prepare the list for later calculation
            hit_dict["global_identity"].append(hsp_dict["hsp_identity"])

            hsp_dict["hsp_length"] = hsp.align_length
            hsp_dict["hsp_bits"] = hsp.bits
            hsp_dict["hsp_evalue"] = hsp.expect
            hsp_dict["query_id"] = query_id
            hsp_dict["target_id"] = target_id
            hsp_dict_list.append(hsp_dict)

        hit_dict.update(kwargs)
        hit_dict["query_id"] = query_id
        hit_dict["target_id"] = target_id

        hit_dict["global_identity"] = mean(hit_dict["global_identity"])

        q_merged_intervals = sorted(merge(q_intervals), key=operator.itemgetter(0, 1))
        q_aligned = sum([tup[1] - tup[0] + 1 for tup in q_merged_intervals])
        hit_dict["query_aligned_length"] = q_aligned
        hit_dict["query_start"] = q_merged_intervals[0][0]
        hit_dict["query_end"] = q_merged_intervals[-1][1]

        t_merged_intervals = sorted(merge(t_intervals), key=operator.itemgetter(0, 1))
        t_aligned = sum([tup[1] - tup[0] + 1 for tup in t_merged_intervals])
        hit_dict["target_aligned_length"] = t_aligned
        hit_dict["target_start"] = t_merged_intervals[0][0]
        hit_dict["target_end"] = t_merged_intervals[-1][1]

        return hit_dict, hsp_dict_list
