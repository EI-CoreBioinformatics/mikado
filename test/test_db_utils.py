"""
TestCase to test the DButils module
"""

import os
import unittest
import mikado.utilities.dbutils
import mikado.configuration.configurator
import mikado.serializers
import sqlalchemy
import sqlalchemy.orm.session
import sqlalchemy.engine.base
import sqlite3
import mikado.serializers.junction
import mikado.serializers.blast_serializer
import mikado.serializers.orf


__author__ = 'Luca Venturini'


class TestDbConnect(unittest.TestCase):

    def setUp(self):
        self.json = mikado.configuration.configurator.to_json(
            os.path.join(os.path.dirname(__file__), "configuration.yaml"))
        self.assertEqual(self.json["db_settings"]["db"], "mikado.db")
        self.json["db_settings"]["db"] = os.path.join(os.path.dirname(__file__),
                                                      self.json["db_settings"]["db"],)


    def test_connector(self):
        connector = mikado.utilities.dbutils.create_connector(self.json)
        self.assertIsInstance(connector, sqlite3.Connection)

    def test_engine(self):

        engine = mikado.utilities.dbutils.connect(self.json)
        self.assertIsInstance(engine, sqlalchemy.engine.base.Engine)
        table_names = ['chrom', 'hit', 'hsp', 'junctions', 'orf', 'query', 'target']
        self.assertEqual(engine.table_names(), table_names)

    def test_content(self):

        engine = mikado.utilities.dbutils.connect(self.json)
        sessionmaker = sqlalchemy.orm.sessionmaker(bind=engine)
        session = sessionmaker()
        # Simple tests based on the static content of the dictionary
        self.assertEqual(session.query(mikado.serializers.junction.Junction).count(), 371)
        self.assertEqual(session.query(mikado.serializers.orf.Orf).count(), 82)
        self.assertEqual(session.query(mikado.serializers.blast_serializer.Target).count(), 38909)
        self.assertEqual(session.query(mikado.serializers.blast_serializer.Query).count(), 95)
        self.assertEqual(session.query(mikado.serializers.blast_serializer.Hit).count(), 576)
        self.assertEqual(session.query(mikado.serializers.blast_serializer.Hsp).count(), 684)

        first_query = session.query(mikado.serializers.blast_serializer.Query).limit(1).one()
        astup = first_query.as_tuple()
        self.assertTrue(astup._fields, ("query_id", "query_name", "query_length"))
        self.assertIsInstance(astup.query_id, int)
        self.assertIsInstance(astup.query_length, int)
        self.assertIsInstance(astup.query_name, str)
        
        first_target = session.query(
            mikado.serializers.blast_serializer.Target).limit(1).one()
        astup = first_target.as_tuple()
        self.assertTrue(astup._fields, ("target_id", "target_name", "target_length"))
        self.assertIsInstance(astup.target_id, int)
        self.assertIsInstance(astup.target_length, int)
        self.assertIsInstance(astup.target_name, str)        

    def test_query_init(self):

        _ = mikado.serializers.blast_serializer.Query("foo", 1000)
        with self.assertRaises(TypeError):
            _ = mikado.serializers.blast_serializer.Query(100, 1000)
        with self.assertRaises(TypeError):
            _ = mikado.serializers.blast_serializer.Query("foo", 0)
        with self.assertRaises(TypeError):
            _ = mikado.serializers.blast_serializer.Query("foo", -10)
        with self.assertRaises(TypeError):
            _ = mikado.serializers.blast_serializer.Query("foo", 1000.0)

    def test_target_init(self):

        _ = mikado.serializers.blast_serializer.Target("foo", 1000)
        with self.assertRaises(TypeError):
            _ = mikado.serializers.blast_serializer.Target(100, 1000)
        with self.assertRaises(TypeError):
            _ = mikado.serializers.blast_serializer.Target("foo", 0)
        with self.assertRaises(TypeError):
            _ = mikado.serializers.blast_serializer.Target("foo", -10)
        with self.assertRaises(TypeError):
            _ = mikado.serializers.blast_serializer.Target("foo", 1000.0)

if __name__ == "__main__":
    unittest.main()
