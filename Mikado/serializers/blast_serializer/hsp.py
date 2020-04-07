"""
This module implements the HSP serialisation class.
"""

from sqlalchemy import Column, String, Integer, Float, ForeignKey, Index
from sqlalchemy.sql.schema import PrimaryKeyConstraint, UniqueConstraint
from sqlalchemy.orm import relationship, column_property
from sqlalchemy import select
from sqlalchemy.ext.hybrid import hybrid_property  # hybrid_method
from ...utilities.dbutils import DBBASE
from . import Query, Target
from .aln_string_parser import prepare_aln_strings


__author__ = 'Luca Venturini'


class Hsp(DBBASE):

    r"""
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

    :return match: the match line between query and target, with the following specs:
        - If the position is a match/positive, keep the original value
        - If the position is a gap *for the query*, insert a - (dash)
        - If the position is a gap *for the target*, insert a _ (underscore)
        - If the position is a gap *for both*, insert a \ (backslash)

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
    hsp_evalue_index = Index('hsp_evalue_idx', 'hsp_evalue', unique=False)
    combined_index = Index("hsp_combined_idx",
                           "query_id",
                           "target_id", unique=False)
    full_index = Index("hsp_full_idx",
                       "counter", "query_id", "target_id", unique=True)
    query_hsp_start = Column(Integer)
    query_hsp_end = Column(Integer)
    query_frame = Column(Integer)
    target_hsp_start = Column(Integer)
    target_hsp_end = Column(Integer)
    target_frame = Column(Integer)
    uni_constraint = UniqueConstraint("query_id", "target_id",
                                      "query_hsp_start", "query_hsp_end",
                                      "target_hsp_start", "target_hsp_end")
    match = Column(String(10000))
    hsp_evalue = Column(Float)
    hsp_bits = Column(Float)
    hsp_identity = Column(Float)
    hsp_positives = Column(Float)
    hsp_length = Column(Integer)

    query_object = relationship(Query, uselist=False)
    target_object = relationship(Target, uselist=False)

    query = column_property(select([Query.query_name]).where(
        Query.query_id == query_id))
    query_length = column_property(select([Query.query_length]).where(
        Query.query_id == query_id))
    target = select([Target.target_name]).where(
        Target.target_id == target_id)
    target_length = select([Target.target_length]).where(
        Target.target_id == target_id)

    __table_args__ = (pk_constraint, query_index, target_index, combined_index, hsp_evalue_index)

    def __init__(self, hsp, counter, query_id, target_id, qmultiplier=1, tmultiplier=1):

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
        hsp_dict, _, _ = prepare_hsp(hsp, counter, qmultiplier=qmultiplier, tmultiplier=tmultiplier)

        for key in hsp_dict:
            setattr(self, key, hsp_dict[key])

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
            "hsp_bits",
            "match",
            "query_frame",
            "target_frame"
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
        state["match"] = self.match

        return state

    @hybrid_property
    def query_hsp_cov(self):
        """This property returns the percentage of the query which is covered by the HSP."""
        val = (self.query_hsp_end - self.query_hsp_start + 1)
        val /= self.query_length
        return val

    @hybrid_property
    def target_hsp_cov(self):
        """This property returns the percentage of the target which is covered by the HSP."""
        val = (self.target_hsp_end - self.target_hsp_start + 1)
        val /= self.target_length
        return val


def prepare_hsp(hsp, counter, off_by_one=False, qmultiplier=1, tmultiplier=1):

    r"""
    Prepare a HSP for loading into the DB.
    The match line will be reworked in the following way:

    - If the position is a match/positive, keep the original value
    - If the position is a gap *for the query*, insert a - (dash)
    - If the position is a gap *for the target*, insert a _ (underscore)
    - If the position is a gap *for both*, insert a \ (backslash)

    :param hsp: An HSP object from Bio.Blast.NCBIXML
    :type hsp: Bio.SearchIO.HSP
    :param counter: a digit that indicates the priority of the HSP in the hit
    :return: hsp_dict, identical_positions, positives
    :rtype: (dict, set, set)
    """

    hsp_dict = dict()
    # We must start from 1, otherwise MySQL crashes as its indices start from 1 not 0
    hsp_dict["counter"] = counter + 1
    hsp_dict["query_hsp_start"] = hsp.query_start
    match, identical_positions, positives = prepare_aln_strings(hsp, off_by_one=off_by_one,
                                                                qmultiplier=qmultiplier)
    hsp_dict["query_hsp_end"] = hsp.query_end - off_by_one
    hsp_dict["query_frame"] = hsp.query_frame
    hsp_dict["target_hsp_start"] = hsp.hit_start
    hsp_dict["target_hsp_end"] = hsp.hit_end
    hsp_dict["target_frame"] = hsp.hit_frame
    hsp_dict["hsp_identity"] = hsp.ident_num / hsp.aln_span * 100
    hsp_dict["hsp_positives"] = (match.count("+") + match.count("|")) / hsp.aln_span * 100
    hsp_dict["match"] = match
    hsp_dict["hsp_length"] = hsp.aln_span
    hsp_dict["hsp_bits"] = hsp.bitscore
    hsp_dict["hsp_evalue"] = hsp.evalue
    return hsp_dict, identical_positions, positives
