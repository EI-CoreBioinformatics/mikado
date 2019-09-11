"""
Basic module with the Query serialiser.
"""

import collections
from sqlalchemy import Column, String, Integer
from ...utilities.dbutils import DBBASE
__author__ = 'Luca Venturini'


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
        if not isinstance(name, str):
            raise TypeError("Invalid name: {0}".format(name))
        if not isinstance(length, int) or length <= 0:
            raise TypeError("Invalid length value: {0}".format(length))

        self.query_name = name
        self.query_length = length

    def as_tuple(self):
        """Quick function to convert the SQLalchemy object into
        a named tuple with the same fields"""
        return self.named_tup(self.query_id, self.query_name, self.query_length)

    def as_dict(self):
        return {
            "query_name": self.query_name,
            "query_length": self.query_length
        }
