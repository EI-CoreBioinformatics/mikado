"""
This module implements a very simple class for serialising the BLAST targets into
a table containing numeric ID, name and length.
"""


import collections
from sqlalchemy import Column, String, Integer
from ...utilities.dbutils import DBBASE
__author__ = 'Luca Venturini'


# These two classes are OK like this, they do not need more public methods!
# pylint: disable=too-few-public-methods
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

        if not isinstance(target_name, str):
            raise TypeError("Invalid name: {0}".format(target_name))
        if not isinstance(target_length, int) or target_length <= 0:
            raise TypeError("Invalid length value: {0}".format(target_length))

        self.target_name = target_name
        self.target_length = target_length

    def as_tuple(self):
        """Quick function to convert the SQLalchemy object
        into a named tuple with the same fields"""
        return self.named_tup(
            self.target_id, self.target_name, self.target_length)

    def as_dict(self):
        return {
            "target_name": self.target_name,
            "target_length": self.target_length
        }

# pylint: enable=too-few-public-methods
