# coding: utf-8

"""This initializer contains the base declaration for all the DB classes of the module."""

from sqlalchemy.engine.reflection import Inspector
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.engine import create_engine
import sqlite3
import logging
import functools

Inspector = Inspector
DBBASE = declarative_base()


def create_connector(json_conf, logger=None):
    """Creator function for the database connection. It necessitates the following information from
    the json_conf dictionary:

    - dbtype (one of sqlite, mysql, postgresql)
    - db (name of the database file, for sqlite, otherwise name of the database)

    If the database is MySQL/PostGreSQL, the method also requires:

    - dbuser
    - dbhost
    - dbpasswd
    - dbport

    These are controlled and added automatically by the json_utils functions.

    :param json_conf: configuration dictionary
    :type json_conf: dict

    :param logger: a logger instance
    :type logger: logging.Logger

    :rtype : MySQLdb.connect | sqlite3.connect | psycopg2.connect

    """

    if logger is None:
        # Create a default null handler
        logger = logging.Logger("null")
        logger.addHandler(logging.NullHandler())

    if json_conf["dbtype"] == "sqlite":
        if json_conf['run_options']['shm'] is False:
            logger.debug("Connecting to %s", json_conf["db"])
            func = sqlite3.connect(database=json_conf["db"], check_same_thread=False)
        else:
            logger.debug("Connecting to %s", json_conf["run_options"]["shm_db"])
            func = sqlite3.connect(database=json_conf["run_options"]["shm_db"],
                                   check_same_thread=False)
    elif json_conf["dbtype"] == "mysql":
        import MySQLdb
        logger.debug("Connecting to MySQL %s", json_conf["run_options"]["db"])
        func = MySQLdb.connect(host=json_conf["dbhost"],
                               user=json_conf["dbuser"],
                               passwd=json_conf["dbpasswd"],
                               db=json_conf["db"],
                               port=json_conf["dbport"])
    elif json_conf["dbtype"] == "postgresql":
        import psycopg2
        logger.debug("Connecting to PSQL %s", json_conf["run_options"]["db"])
        func = psycopg2.connect(
            host=json_conf["dbhost"],
            user=json_conf["dbuser"],
            password=json_conf["dbpasswd"],
            database=json_conf["db"],
            port=json_conf["dbport"]
        )
    else:
        raise ValueError("DB type not supported! {0}".format(json_conf["dbtype"]))
    return func


def connect(json_conf, logger=None):

    """
    Function to create an engine to connect to a DB with, using the
    configuration inside the provided json_conf.
    :param json_conf:
    :param logger:
    :return:
    """

    if json_conf is None:
        return create_engine("sqlite://:memory:")

    db_connection = functools.partial(create_connector, json_conf, logger=logger)
    engine = create_engine("{0}://".format(json_conf["dbtype"]),
                           creator=db_connection)
    return engine
