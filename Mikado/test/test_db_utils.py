"""
TestCase to test the DButils module
"""

import os
import unittest
import Mikado.utilities.dbutils
import Mikado.configuration.configurator
import Mikado.serializers
import sqlalchemy
import sqlalchemy.orm.session
import sqlalchemy.engine.base
import sqlite3
import Mikado.serializers.junction
import Mikado.serializers.blast_serializer
import Mikado.serializers.orf
import shutil


__author__ = 'Luca Venturini'


class TestDbConnect(unittest.TestCase):

    def setUp(self):
        self.json = Mikado.configuration.configurator.to_json(
            os.path.join(os.path.dirname(__file__), "configuration.yaml"))
        self.assertEqual(self.json["db_settings"]["db"],
                         os.path.join(
                             os.path.dirname(__file__),
                             self.json["db_settings"]["db"],))

    def test_connector(self):
        connector = Mikado.utilities.dbutils.create_connector(self.json)
        self.assertIsInstance(connector, sqlite3.Connection)

    def test_engine(self):

        engine = Mikado.utilities.dbutils.connect(self.json)
        self.assertIsInstance(engine, sqlalchemy.engine.base.Engine)
        table_names = ['chrom', 'hit', 'hsp', 'junctions', 'orf', 'query', 'target']
        self.assertEqual(engine.table_names(), table_names)

    def test_content(self):

        engine = Mikado.utilities.dbutils.connect(self.json)
        sessionmaker = sqlalchemy.orm.sessionmaker(bind=engine)
        session = sessionmaker()
        # Simple tests based on the static content of the dictionary
        self.assertEqual(session.query(Mikado.serializers.junction.Junction).count(), 371)
        self.assertEqual(session.query(Mikado.serializers.orf.Orf).count(), 80)
        self.assertEqual(session.query(Mikado.serializers.blast_serializer.Target).count(), 38909)
        self.assertEqual(session.query(Mikado.serializers.blast_serializer.Query).count(), 94)
        self.assertEqual(session.query(Mikado.serializers.blast_serializer.Hit).count(), 568)
        self.assertEqual(session.query(Mikado.serializers.blast_serializer.Hsp).count(), 676)

        first_query = session.query(Mikado.serializers.blast_serializer.Query).limit(1).one()
        astup = first_query.as_tuple()
        self.assertTrue(astup._fields, ("query_id", "query_name", "query_length"))
        self.assertIsInstance(astup.query_id, int)
        self.assertIsInstance(astup.query_length, int)
        self.assertIsInstance(astup.query_name, str)
        
        first_target = session.query(
            Mikado.serializers.blast_serializer.Target).limit(1).one()
        astup = first_target.as_tuple()
        self.assertTrue(astup._fields, ("target_id", "target_name", "target_length"))
        self.assertIsInstance(astup.target_id, int)
        self.assertIsInstance(astup.target_length, int)
        self.assertIsInstance(astup.target_name, str)        

    def test_query_init(self):

        _ = Mikado.serializers.blast_serializer.Query("foo", 1000)
        with self.assertRaises(TypeError):
            _ = Mikado.serializers.blast_serializer.Query(100, 1000)
        with self.assertRaises(TypeError):
            _ = Mikado.serializers.blast_serializer.Query("foo", 0)
        with self.assertRaises(TypeError):
            _ = Mikado.serializers.blast_serializer.Query("foo", -10)
        with self.assertRaises(TypeError):
            _ = Mikado.serializers.blast_serializer.Query("foo", 1000.0)

    def test_target_init(self):

        _ = Mikado.serializers.blast_serializer.Target("foo", 1000)
        with self.assertRaises(TypeError):
            _ = Mikado.serializers.blast_serializer.Target(100, 1000)
        with self.assertRaises(TypeError):
            _ = Mikado.serializers.blast_serializer.Target("foo", 0)
        with self.assertRaises(TypeError):
            _ = Mikado.serializers.blast_serializer.Target("foo", -10)
        with self.assertRaises(TypeError):
            _ = Mikado.serializers.blast_serializer.Target("foo", 1000.0)

    def test_wrong_db(self):
        self.json["db_settings"]["dbtype"] = "sqlite_foo"
        with self.assertRaises(ValueError):
            _ = Mikado.utilities.dbutils.create_connector(self.json)

    @unittest.skipUnless(os.path.exists("/dev/shm"),
                         "/dev/shm is not available on this system.")
    def test_connect_to_shm(self):
        self.json["pick"]["run_options"]['shm'] = True
        shutil.copy(self.json["db_settings"]["db"], "/dev/shm/")
        self.json["pick"]["run_options"]['shm_db'] = os.path.join(
            "/dev/shm/",
            self.json["db_settings"]["db"])
        connector = Mikado.utilities.dbutils.connect(self.json)
        self.assertEqual(str(connector.url), "sqlite://")
        engine = Mikado.utilities.dbutils.connect(self.json)
        sessionmaker = sqlalchemy.orm.sessionmaker(bind=engine)
        session = sessionmaker()
        first_target = session.query(
            Mikado.serializers.blast_serializer.Target).limit(1).one()
        astup = first_target.as_tuple()
        self.assertTrue(astup._fields, ("target_id", "target_name", "target_length"))

    def test_to_memory(self):
        connector = Mikado.utilities.dbutils.connect(None)
        self.assertEqual(str(connector.url), "sqlite:///:memory:")


if __name__ == "__main__":
    unittest.main()
