from sqlalchemy.engine.reflection import Inspector
from sqlalchemy.ext.declarative import declarative_base

"""This initializer contains the base declaration for all the DB classes of the module."""

Inspector=Inspector
dbBase = declarative_base()
