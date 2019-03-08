import unittest
# import Mikado
from .. import utilities, configuration, serializers
import tempfile
import sqlalchemy.orm
import pandas as pd
import numpy as np
from ..utilities.dbutils import DBBASE
from ..serializers.external import External, ExternalSource, ExternalSerializer
from ..serializers.orf import Query
import sqlalchemy.exc


__author__ = 'Luca Venturini'


class TestExternal(unittest.TestCase):

    logger = utilities.log_utils.create_null_logger("test_junction")

    def setUp(self):
        self.dbfile = tempfile.mktemp(suffix=".db")
        self.json_conf = configuration.configurator.to_json(None)
        self.json_conf["db_settings"]["dbtype"] = "sqlite"
        self.json_conf["db_settings"]["db"] = self.dbfile
        self.__create_session()

    def __create_session(self):
        self.engine = utilities.dbutils.connect(
            self.json_conf, self.logger)
        self.sessionmaker = sqlalchemy.orm.sessionmaker(bind=self.engine)
        self.session = self.sessionmaker()
        DBBASE.metadata.create_all(self.engine)

    def test_serialize_source(self):

        source = ExternalSource("foo", np.dtype("float"), False)
        self.session.add(source)
        self.session.commit()

        self.assertEqual(self.session.query(serializers.external.ExternalSource).count(), 1)

    def test_add_score(self):

        source = ExternalSource("cdna_length", np.dtype("int"), False)
        query = Query("foo.1", 200)
        self.session.add(source)
        self.session.add(query)
        self.session.commit()

        score = External(query.query_id, source.source_id, 200)
        self.session.add(score)

        # Before commiting, everything is detached
        self.assertIs(score.source, None)
        self.session.commit()

        self.assertEqual(score.source, source.source)
        self.assertEqual(score.query, query.query_name)
        self.assertEqual(self.session.query(serializers.external.External).count(), 1)

    def test_wrong_score(self):
        source = ExternalSource("cdna_length", np.dtype("int"), False)
        query = Query("foo.1", 200)
        self.session.add(source)
        self.session.add(query)
        self.session.commit()

        # This fails because the query_id
        score = External(query.query_id + 1, source.source_id, 200)
        self.session.add(score)
        with self.assertRaises(sqlalchemy.exc.IntegrityError):
            self.session.commit()
        self.session.rollback()

        score = External(query.query_id, source.source_id + 1, 200)
        self.session.add(score)
        with self.assertRaises(sqlalchemy.exc.IntegrityError):
            self.session.commit()
        self.session.rollback()

        with self.assertRaises(sqlalchemy.exc.ArgumentError):
            _ = External(query_id=query.query_id, source_id=source.source_id, score="200")

    def test_load_table(self):

        queries = ["T{}".format(_) for _ in range(1, 20)]
        [self.session.add(Query(_, 100)) for _ in queries]
        self.session.commit()

        df = pd.DataFrame({'tid': queries,
                           'foo': [100 for _ in queries],
                           'bar': [10 for _ in queries]})
        df.set_index("tid", inplace=True)
        with tempfile.NamedTemporaryFile() as temp_df:  #, tempfile.NamedTemporaryFile() as fasta_tmp:
            # for query in queries:
            #     print(">{}".format(query), file=fasta_tmp)
            #     print("A" * 100, file=fasta_tmp)
            # fasta_tmp.flush()

            df.to_csv(temp_df.name, sep="\t", index_label="tid", mode="wt", header=True)
            temp_df.flush()
            self.json_conf["serialise"]["files"]["transcripts"] = None
            serializer = ExternalSerializer(temp_df.name, self.json_conf)
            serializer.serialize()
