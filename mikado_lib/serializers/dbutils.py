# coding: utf-8

"""This initializer contains the base declaration for all the DB classes of the module."""

from sqlalchemy.engine.reflection import Inspector
from sqlalchemy.ext.declarative import declarative_base


Inspector = Inspector
dbBase = declarative_base()
