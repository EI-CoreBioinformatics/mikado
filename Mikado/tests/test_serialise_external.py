import logging
import os
import unittest
import Mikado
import tempfile
import sqlalchemy.orm
from sqlalchemy import and_  # , or_
from pkg_resources import resource_stream
import gzip

__author__ = 'Luca Venturini'


class TestLoadExternal(unittest.TestCase):

    logger = Mikado.utilities.log_utils.create_null_logger("test_junction")
    dbfile = tempfile.mktemp(suffix=".db")

