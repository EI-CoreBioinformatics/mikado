"""
TestCase to test the DButils module
"""

import os
import unittest
import mikado_lib.utilities.dbutils
import mikado_lib.configuration.configurator
import mikado_lib.serializers
import sqlalchemy
import sqlalchemy.orm.session
import sqlalchemy.engine.base
import sqlite3
import mikado_lib.serializers.junction
import mikado_lib.serializers.blast_serializer
import mikado_lib.serializers.orf


__author__ = 'Luca Venturini'


class TestDbConnect(unittest.TestCase):

    def setUp(self):
        self.json = mikado_lib.configuration.configurator.to_json(
            os.path.join(os.path.dirname(__file__), "configuration.yaml"))
        self.assertEqual(self.json["db_settings"]["db"], "mikado.db")
        self.json["db_settings"]["db"] = os.path.join(os.path.dirname(__file__),
                                                      self.json["db_settings"]["db"],)


    def test_connector(self):
        connector = mikado_lib.utilities.dbutils.create_connector(self.json)
        self.assertIsInstance(connector, sqlite3.Connection)

    def test_engine(self):

        engine = mikado_lib.utilities.dbutils.connect(self.json)
        self.assertIsInstance(engine, sqlalchemy.engine.base.Engine)
        table_names = ['chrom', 'hit', 'hsp', 'junctions', 'orf', 'query', 'target']
        self.assertEqual(engine.table_names(), table_names)

    def test_content(self):

        engine = mikado_lib.utilities.dbutils.connect(self.json)
        sessionmaker = sqlalchemy.orm.sessionmaker(bind=engine)
        session = sessionmaker()
        # Simple tests based on the static content of the dictionary
        self.assertEqual(session.query(mikado_lib.serializers.junction.Junction).count(), 371)
        self.assertEqual(session.query(mikado_lib.serializers.orf.Orf).count(), 82)
        self.assertEqual(session.query(mikado_lib.serializers.blast_serializer.Target).count(), 38909)
        self.assertEqual(session.query(mikado_lib.serializers.blast_serializer.Query).count(), 95)
        self.assertEqual(session.query(mikado_lib.serializers.blast_serializer.Hit).count(), 576)
        self.assertEqual(session.query(mikado_lib.serializers.blast_serializer.Hsp).count(), 684)


if __name__ == "__main__":
    unittest.main()
