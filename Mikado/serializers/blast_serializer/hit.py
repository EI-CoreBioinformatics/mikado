"""
This module implements the Hit serialisation class.
"""

from ...exceptions import InvalidHit
from sqlalchemy import Column, Integer, Float, ForeignKey, Index
from sqlalchemy.sql.schema import PrimaryKeyConstraint
from sqlalchemy.orm import relationship, column_property, backref
from sqlalchemy import select
from sqlalchemy.ext.hybrid import hybrid_property  # hybrid_method
from ...utilities.dbutils import DBBASE
from . import Query, Target, Hsp, prepare_hsp
import numpy as np
from ...parsers.blast_utils import merge


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
    def __init__(self, query_id, target_id, query_length, alignment,
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

        #qmulti = kwargs["query_multiplier"]
        # tmulti = kwargs["target_multiplier"]
        # qlength = kwargs["query_length"]
        prepared_hit, _ = prepare_hit(alignment, query_id, target_id,
                                      query_multiplier=query_multiplier,
                                      target_multiplier=target_multiplier,
                                      query_length=query_length)
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
        state["query_hit_ratio"] = (query_tuple.query_length / hit_tuple.query_multiplier) / \
                                   (target_tuple.target_length / hit_tuple.target_multiplier)
        state["hit_query_ratio"] = (target_tuple.target_length / hit_tuple.target_multiplier) / \
                                   (query_tuple.query_length / hit_tuple.query_multiplier)
        state["query_cov"] = state["query_aligned_length"] / query_tuple.query_length
        assert state["query_cov"] <= 1, (state,)
        state["target_cov"] = state["target_aligned_length"] / target_tuple.target_length
        assert state["target_cov"] <= 1, (state,)

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
        state["query_hit_ratio"] = state["query_length"] / state["query_multiplier"] /\
            (state["target_length"] / state["target_multiplier"])

        state["hit_query_ratio"] = state["target_length"] / state["target_multiplier"] /\
            (state["query_length"] / state["query_multiplier"])

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

        ratio = self.query_length / self.query_multiplier
        ratio /= self.target_length / self.target_multiplier

        return ratio

    @hybrid_property
    def hit_query_ratio(self):
        """
        This property returns the quotient (Target Length)/(Query Length)
        """

        ratio = self.target_length / self.target_multiplier
        ratio /= (self.query_length / self.query_multiplier)

        return ratio


hit_cols = [col.name for col in Hit.__table__.columns]
hsp_cols = [col.name for col in Hsp.__table__.columns]


def prepare_hit(hit, query_id, target_id, off_by_one=False, as_list=False, **kwargs):
    """Prepare the dictionary for fast loading of Hit and Hsp objects.
    global_positives: the similarity rate for the global hit *using the query perspective*
    global_identity: the identity rate for the global hit *using the query perspective*

    :param hit: the hit to parse.
    :type hit: Bio.SearchIO.Hit

    :param query_id: the numeric ID of the query in the database. Necessary for serialisation.
    :type query_id: int

    :param target_id: the numeric ID of the target in the database. Necessary for serialisation.
    :type target_id: int

    :param kwargs: additional properties to give to the hit_dict. Retrieved e.g. from descriptions.
    :type kwargs: dict
    """

    hit_dict = dict()
    hsp_dict_list = []
    q_intervals = []
    t_intervals = []

    qmulti = kwargs["query_multiplier"]
    tmulti = kwargs["target_multiplier"]
    qlength = kwargs["query_length"]
    assert isinstance(qmulti, (int, float)), type(qmulti)
    assert isinstance(tmulti, (int, float)), type(tmulti)
    hit_dict.update(kwargs)
    hit_dict["query_id"] = query_id
    hit_dict["target_id"] = target_id

    query_array = np.zeros([2, int(qlength)], dtype=np.int)

    for counter, hsp in enumerate(hit.hsps):
        if hsp.query_start + off_by_one - 1 > qlength:
            raise ValueError("Invalid length: {}, {}", hsp.query_start + off_by_one - 1, qlength)

        hsp_dict, ident, posit = prepare_hsp(hsp, counter, qmultiplier=qmulti, tmultiplier=tmulti,
                                             off_by_one=off_by_one)
        if ident.max() > query_array.shape[1] or ident.min() < 0:
            raise IndexError("Invalid indexing values (max {}; frame {}; hsp: {})!"\
"Too low: {}\nToo high: {}".format(
                query_array.shape[1],
                hsp.query_frame, hsp.__dict__,
                list(ident[(ident < 0)]),
                list(ident[ident > query_array.shape[1]])))
        try:
            query_array[0, ident - off_by_one] = 1
        except IndexError as exc:
            raise IndexError("{}, off by one: {}; min, max: {}, {}; hsp {}".format(
                exc, off_by_one, ident.min(), ident.max(), hsp._items[0].__dict__))
        try:
            query_array[1, posit - off_by_one] = 1
        except IndexError as exc:
            raise IndexError("{}, off by one: {}; min, max: {}, {}; frame: {}".format(
                exc, off_by_one, posit.min(), posit.max(), hsp.query_frame))
        hsp_dict["query_id"] = query_id
        hsp_dict["target_id"] = target_id
        hsp_dict_list.append(hsp_dict)
        q_intervals.append((hsp.query_start, hsp.query_end - 1))
        t_intervals.append((hsp.hit_start, hsp.hit_end - 1))

    q_merged_intervals, q_aligned = merge(q_intervals)
    hit_dict["query_aligned_length"] = min(qlength, q_aligned)
    qstart, qend = q_merged_intervals[0][0], q_merged_intervals[-1][1]
    hit_dict["query_start"], hit_dict["query_end"] = qstart, qend
    identical = np.where(query_array[0] == 1)[0]
    positives = np.where(query_array[1] == 1)[0]

    if identical.shape[0] > q_aligned:
        raise ValueError(
            "Number of identical positions ({}) greater than number of aligned positions ({})!\n{}\n{}".format(
            identical.shape[0], q_aligned, q_intervals, q_merged_intervals))

    if positives.shape[0] > q_aligned:
        raise ValueError("Number of identical positions ({}) greater than number of aligned positions ({})!".format(
            positives.shape[0], q_aligned))

    t_merged_intervals, t_aligned = merge(t_intervals)
    hit_dict["target_aligned_length"] = min(t_aligned, hit.seq_len)
    hit_dict["target_start"] = t_merged_intervals[0][0]
    hit_dict["target_end"] = t_merged_intervals[-1][1] + 1
    hit_dict["global_identity"] = identical.shape[0] * 100 / q_aligned
    hit_dict["global_positives"] = positives.shape[0] * 100 / q_aligned

    if as_list is True:
        hit_list = [hit_dict[col] for col in hit_cols]
        hsp_dict_list = [[hsp[col] for col in hsp_cols] for hsp in hsp_dict_list]
        return hit_list, hsp_dict_list
    else:
        return hit_dict, hsp_dict_list
