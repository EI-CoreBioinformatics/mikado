# coding: utf-8

"""This initializer contains the base declaration for all the DB classes of the module."""

from sqlalchemy.engine.reflection import Inspector
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.engine import create_engine

Inspector = Inspector
dbBase = declarative_base()


def connect(json_conf):
    """
    Wrapper to create a working connection to the database.

    :param json_conf: configuration dictionary
    :type json_conf: dict

    :rtype: sqlalchemy.engine.base.Engine
    """
    if json_conf["dbtype"] == "sqlite":
        engine = create_engine("sqlite:///{0}".format(json_conf["db"]))
    else:
        engine = create_engine("{dbtype}://{dbuser}:{dbpasswd}@{dbhost}/{db}".format(
            dbtype=json_conf["dbtype"],
            dbuser=json_conf["dbuser"],
            dbpasswd=json_conf["dbpasswd"] if json_conf["dbpasswd"] is not None else "",
            dbhost=json_conf["dbhost"],
            db=json_conf["db"]))
    return engine
