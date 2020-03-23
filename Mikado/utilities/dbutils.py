# coding: utf-8

"""This initializer contains the base declaration for all the DB classes of the module."""

from sqlalchemy.engine.reflection import Inspector
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.engine import create_engine, Engine
from sqlalchemy import event
from sqlalchemy_utils import database_exists, create_database
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

    db_settings = json_conf["db_settings"]

    func = None
    if db_settings["dbtype"] == "sqlite":
        if not database_exists("sqlite:///{}".format(db_settings["db"])):
            logger.debug("No database found, creating a mock one")
            create_database("sqlite:///{}".format(db_settings["db"]))
        logger.debug("Connecting to %s", db_settings["db"])
        func = sqlite3.connect(database=db_settings["db"], check_same_thread=False)
    elif db_settings["dbtype"] in ("mysql", "postgresql"):
        if db_settings["dbpasswd"] != '':
            passwd = ":{0}".format(db_settings["dbpasswd"])
        else:
            passwd = ''
        url = "{dialect}://{user}{passwd}@{host}:{port}/{db}".format(
            dialect=db_settings["dbtype"],
            host=db_settings["dbhost"],
            user=db_settings["dbuser"],
            passwd=passwd,
            db=db_settings["db"],
            port=db_settings["dbport"]
        )
        if database_exists(url) is False:
            create_database(url)

        if db_settings["dbtype"] == "mysql":
            import MySQLdb
            logger.debug("Connecting to MySQL %s", db_settings["db"])
            func = MySQLdb.connect(host=db_settings["dbhost"],
                                   user=db_settings["dbuser"],
                                   passwd=db_settings["dbpasswd"],
                                   db=db_settings["db"],
                                   port=db_settings["dbport"])
        elif db_settings["dbtype"] == "postgresql":
            import psycopg2
            logger.debug("Connecting to PSQL %s", db_settings["db"])
            func = psycopg2.connect(
                host=db_settings["dbhost"],
                user=db_settings["dbuser"],
                password=db_settings["dbpasswd"],
                database=db_settings["db"],
                port=db_settings["dbport"]
            )
    else:
        raise ValueError("DB type not supported! {0}".format(db_settings["dbtype"]))
    return func


def connect(json_conf, logger=None, **kwargs):

    """
    Function to create an engine to connect to a DB with, using the
    configuration inside the provided json_conf.
    :param json_conf:
    :param logger:
    :return: sqlalchemy.engine.base.Engine
    """

    @event.listens_for(Engine, "connect")
    def set_sqlite_pragma(dbapi_connection, connection_record):
        cursor = dbapi_connection.cursor()
        try:
            cursor.execute("PRAGMA foreign_keys=ON")
            cursor.execute("PRAGMA synchronous=OFF")
            cursor.execute("PRAGMA temp_store=MEMORY")
            cursor.execute("PRAGMA journal_mode=MEMORY")
            cursor.execute("PRAGMA count_changes=OFF")
        except sqlite3.OperationalError:
            pass
        finally:
            cursor.close()

    if json_conf is None:
        return create_engine("sqlite:///:memory:", **kwargs)

    db_connection = functools.partial(create_connector, json_conf, logger=logger)
    engine = create_engine("{0}://".format(json_conf["db_settings"]["dbtype"]),
                           creator=db_connection, **kwargs)
    DBBASE.metadata.create_all(engine, checkfirst=True)

    return engine
