"""
This module implements the Hit serialisation class.
"""

from sqlalchemy import Column, Integer, Float, ForeignKey, Index
from sqlalchemy.sql.schema import PrimaryKeyConstraint
from sqlalchemy.orm import relationship, column_property, backref
from sqlalchemy import select
from sqlalchemy.ext.hybrid import hybrid_property  # hybrid_method
from ...utilities.dbutils import DBBASE
from . import Query, Target, prepare_hit, Hsp

__author__ = 'Luca Venturini'


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
    evalue_index = Index('hit_evalue_idx', 'evalue', unique=False)
    evalue = Column(Float)
    bits = Column(Float)
    global_identity = Column(Float)
    global_positives = Column(Float)
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
                                backref=backref("hits", cascade="all, delete-orphan"))
    target_object = relationship(Target,
                                 uselist=False,
                                 lazy="immediate",
                                 backref=backref("hits", cascade="all, delete-orphan"))

    join_condition = "and_(Hit.query_id==Hsp.query_id, Hit.target_id==Hsp.target_id)"
    hsps = relationship(Hsp, uselist=True,
                        # lazy="immediate",
                        lazy="subquery",
                        backref=backref("hit_object", uselist=False),
                        cascade="all, delete-orphan",
                        single_parent=True,
                        foreign_keys=[query_id, target_id],
                        primaryjoin=join_condition)

    query = column_property(select([Query.query_name]).where(
        Query.query_id == query_id))
    query_length = column_property(select([Query.query_length]).where(
        Query.query_id == query_id))
    target = select([Target.target_name]).where(
        Target.target_id == target_id)
    target_length = select([Target.target_length]).where(
        Target.target_id == target_id)

    __table_args__ = (qt_constraint, qt_index, query_index, target_index, evalue_index)

    # All arguments are necessary and it is more convenient to have them here
    # rather than in a struct/dict/whatever
    # pylint: disable=too-many-arguments
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

        prepared_hit, _ = prepare_hit(alignment, query_id, target_id)
        for key in ["global_identity", "global_positives", "query_aligned_length",
                    "query_start", "query_end", "target_aligned_length",
                    "target_start", "target_end"]:
            setattr(self, key, prepared_hit[key])
    # pylint: enable=too-many-arguments

    def __str__(self):
        line = [self.query, self.target, self.evalue,
                self.bits, self.query_start, self.query_end,
                self.target_start, self.target_end,
                self.query_length, self.target_length]

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
            "global_positives",
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
        state["query_length"] = query_tuple.query_length
        state["target_length"] = target_tuple.target_length
        state["query_hit_ratio"] = (query_tuple.query_length * hit_tuple.query_multiplier) / \
                                   (target_tuple.target_length * hit_tuple.target_multiplier)
        state["hit_query_ratio"] = (target_tuple.target_length * hit_tuple.target_multiplier) / \
                                   (query_tuple.query_length * hit_tuple.query_multiplier)
        state["hsps"] = []
        for hsp in hsps:
            hsp["query_hsp_cov"] = (hsp["query_hsp_end"] - hsp["query_hsp_start"] + 1)
            hsp["query_hsp_cov"] /= (state["query_length"])
            hsp["target_hsp_cov"] = (hsp["target_hsp_end"] - hsp["target_hsp_start"] + 1)
            hsp["target_hsp_cov"] /= (state["target_length"])
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
        state["query_length"] = query_object.query_length
        state["target_length"] = target_object.target_length
        state["query_hit_ratio"] = state["query_length"] * state["query_multiplier"] /\
            (state["target_length"] * state["target_multiplier"])

        state["hit_query_ratio"] = state["target_length"] * state["target_multiplier"] /\
            (state["query_length"] * state["query_multiplier"])

        state["hsps"] = []
        for hsp in self.hsps:
            dict_hsp = hsp.as_dict()
            dict_hsp["query_hsp_cov"] = (dict_hsp["query_hsp_end"] -
                                         dict_hsp["query_hsp_start"] + 1)
            dict_hsp["query_hsp_cov"] /= state["query_length"]
            dict_hsp["target_hsp_cov"] = (dict_hsp["target_hsp_end"] -
                                          dict_hsp["target_hsp_start"] + 1)
            dict_hsp["target_hsp_cov"] /= state["target_length"]
            state["hsps"].append(dict_hsp)

        return state

    @hybrid_property
    def query_hit_ratio(self):
        """
        This property returns the quotient (Query Length)/(Target Length)
        """

        ratio = self.query_length * self.query_multiplier
        ratio /= self.target_length * self.target_multiplier

        return ratio

    @hybrid_property
    def hit_query_ratio(self):
        """
        This property returns the quotient (Target Length)/(Query Length)
        """

        ratio = self.target_length * self.target_multiplier
        ratio /= (self.query_length * self.query_multiplier)

        return ratio
