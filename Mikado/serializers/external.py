# coding: utf-8

"""
This module is necessary to serialise external data.
The initial model is very simple - each external data
will have a tag, and internally the files will be TAB-delimited
2-column files, one for the TID and one for the score.
"""

import os
import sqlite3
import pyfaidx
from sqlalchemy import Column, String, Integer, ForeignKey, CHAR, Index, Float, Boolean
from sqlalchemy.sql.schema import PrimaryKeyConstraint
from sqlalchemy.ext.hybrid import hybrid_property
import sqlalchemy.exc
from sqlalchemy.orm import relationship, backref, column_property
from sqlalchemy.orm.session import Session  # sessionmaker
from sqlalchemy import select
from ..utilities.dbutils import DBBASE, Inspector, connect
from .blast_serializer import Query
from ..utilities.log_utils import create_null_logger, check_logger


class ExternalSources(DBBASE):

    __tablename__ = "external_sources"

    source_id = Column(Integer, primary_key=True)
    source = Column(String)

    def __init__(self, source):

        self.source = source


class External(DBBASE):

    """This class serialises transcript ids, source of the score, and score itself."""

    __tablename__ = "external"

    query_id = Column(Integer, ForeignKey(Query.query_id), unique=False)
    source_id = Column(Integer, ForeignKey(ExternalSources.source_id), unique=False)
    ext_constraint = PrimaryKeyConstraint("query_id", "source_id", name="source_key")
    source = column_property(select([ExternalSources.source]).where(
        ExternalSources.source_id == source_id))
    score = Column(Float)

    query = column_property(select([Query.query_name]).where(
        Query.query_id == query_id))

    __table_args__ = (ext_constraint, )

    def __init__(self, source_id, score):

        self.source_id = source_id
        self.score = score
