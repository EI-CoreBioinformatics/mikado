import asyncio
from Mikado._transcripts.scoring_configuration import MinMaxScore
from Mikado.configuration.configurator import load_and_validate_config
from Mikado.loci import Superlocus
from Mikado.serializers.external import External, ExternalSource
from Mikado.serializers.blast_serializer.query import Query
from Mikado.transcripts import Transcript
from Mikado.utilities.dbutils import DBBASE
import unittest
from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine
import tempfile


class AsyncExternalTest(unittest.TestCase):

    def test_get_external(self):
        checked_conf = load_and_validate_config(None).copy()
        checked_conf.pick.output_format.report_all_external_metrics = True
        transcript = Transcript()
        transcript.chrom = "15"
        transcript.source = "protein_coding"
        transcript.start = 47631264
        transcript.end = 48051999

        exons = [(47631264, 47631416),
                 (47704590, 47704669),
                 (47762671, 47762742),
                 (47893062, 47893093),
                 (47895572, 47895655),
                 (48051942, 48051999)]

        transcript.strand = "+"
        transcript.add_exons(exons)
        transcript.id = "ENST00000560636"
        transcript.parent = "ENSG00000137872"
        transcript2 = transcript.copy()
        transcript2.id = "ENST00000560637"
        checked_conf.scoring.scoring["attributes.tpm"] = MinMaxScore.Schema().load(
            {"rescaling": "max", "default": 0, "rtype": "float", 'multiplier': 4,
             'use_raw': True, 'percentage': True})
        transcript.attributes["tpm"] = 10

        db = DBBASE
        int_source = ExternalSource('int', 'int', 0)
        float_source = ExternalSource('float', 'float', 0)
        bool_source = ExternalSource('bool', 'bool', 0)

        raw_int_source = ExternalSource('raw_int', 'int', 1)
        raw_float_source = ExternalSource('raw_float', 'float', 1)
        raw_bool_source = ExternalSource('raw_bool', 'bool', 1)

        int_score = External(1, 1, 10)
        float_score = External(1, 2, 10.0)
        bool_score = External(1, 3, int(False))  # We cast as int here following external.py serialize function

        raw_int_score = External(1, 4, 8)
        raw_float_score = External(1, 5, 8.0)
        raw_bool_score = External(1, 6, int(True))  # We cast as int here following external.py serialize function

        query = Query(transcript.id, transcript.cdna_length)
        query2 = Query(transcript2.id, transcript2.cdna_length)

        dbname = tempfile.NamedTemporaryFile(suffix=".db", delete=False)
        engine = create_engine("sqlite:///:memory:")
        db.metadata.create_all(engine)
        print(engine.url, dbname.name)
        SessionMaker = sessionmaker(bind=engine)
        session = SessionMaker()
        session.add_all([int_source, float_source, bool_source, raw_int_source, raw_float_source, raw_bool_source])
        session.add_all([query, query2])
        session.add_all([int_score, float_score, bool_score, raw_int_score, raw_float_score, raw_bool_score])
        session.commit()
        sup = Superlocus(transcript, configuration=checked_conf)
        sup.session = session
        tid = transcript.id
        self.assertIn(tid, sup.transcripts)
        from collections import namedtuple
        qobj = {1: namedtuple('t', field_names=('query_name'))}
        qobj[1].query_name = 'ENST00000560636'
        external = asyncio.run(sup.get_external(qobj, [1]))

        self.assertEqual(external, {
            'ENST00000560636': {
                   'int': (10, False),
                   'float': (10.0, False),
                   'bool': (False, False),
                   'raw_int': (8, True),
                   'raw_float': (8.0, True),
                   'raw_bool': (True, True)
            }
        })
