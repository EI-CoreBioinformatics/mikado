import asyncio
import itertools
from Mikado._transcripts.scoring_configuration import MinMaxScore, SizeFilter
from Mikado.configuration.configurator import load_and_validate_config
from Mikado.loci import Superlocus
from Mikado.parsers.bed12 import BED12
from Mikado.serializers.blast_serializer import Target, Hit, Hsp
from Mikado.serializers.external import External, ExternalSource
from Mikado.serializers.blast_serializer.query import Query
from Mikado.serializers.orf import Orf
from Mikado.serializers.junction import Chrom, Junction
from Mikado.transcripts import Transcript
from Mikado.utilities.dbutils import DBBASE as db
import unittest
from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine


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

        engine = create_engine("sqlite:///:memory:")
        db.metadata.create_all(engine)
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

        sup.configuration.pick.output_format.report_all_external_metrics = False
        external = asyncio.run(sup.get_external(qobj, [1]))
        self.assertEqual(len(external), 0)
        # These are meaningless it's just to verify we are loading *only* these metrics.
        # We should *NOT* have 'float' as it is not present in any section.
        sup.configuration.scoring.scoring["external.int"] = MinMaxScore(rescaling="max", filter=None)
        sup.configuration.scoring.requirements.parameters["external.raw_float"] = SizeFilter(operator="gt",
                                                                                             value=100)
        sup.configuration.scoring.cds_requirements.parameters["external.raw_int"] = SizeFilter(operator="lt",
                                                                                                 value=1)
        sup.configuration.scoring.as_requirements.parameters["external.raw_bool"] = SizeFilter(operator="lt",
                                                                                                 value=1)
        sup.configuration.scoring.not_fragmentary.parameters["external.bool"] = SizeFilter(operator="ne",
                                                                                           value=False)
        external = asyncio.run(sup.get_external(qobj, [1]))
        self.assertEqual(external, {
            'ENST00000560636': {
                'int': (10, False),
                'raw_float': (8.0, True),
                'bool': (False, False),
                'raw_int': (8, True),
                'raw_bool': (True, True)
            }
        })


class AsyncJunctionTest(unittest.TestCase):

    def test_retrieval(self):
        engine = create_engine("sqlite:///:memory:")
        db.metadata.create_all(engine)
        SessionMaker = sessionmaker(bind=engine)
        session = SessionMaker()

        transcript = Transcript(accept_undefined_multi=True)
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

        chrom_one = Chrom("1", 10**8)
        chrom_fifteen = Chrom("15", 5 * 10 ** 8)
        session.add_all([chrom_one, chrom_fifteen])
        session.commit()
        # junction_start, junction_end, name, strand, score, chrom_id)
        # This junction is on a different chrom
        junction_chrom_one = Junction(47704669 + 1, 47762671 - 1, "chrom_one", "+", 10, chrom_one.chrom_id)
        # This junction is too far away
        outside_chrom_15 = Junction(47704669 - 10 ** 6 + 1, 47762671 - 10 ** 6 - 1, "chrom_15_outside", "+", 10,
                                    chrom_fifteen.chrom_id)
        # This junction is in the right place but wrong strand
        wrong_strand_chrom_15 = Junction(47704669 + 1, 47762671 - 1, "chrom_15_wrong_strand", "-", 10,
                                         chrom_fifteen.chrom_id)
        # This one is correct
        chrom_15_junction = Junction(47704669 + 1, 47762671 - 1, "chrom_15", "+", 10, chrom_fifteen.chrom_id)
        session.add_all([junction_chrom_one, outside_chrom_15, wrong_strand_chrom_15, chrom_15_junction])
        session.commit()

        self.assertEqual(junction_chrom_one.chrom, "1")
        for junc in [outside_chrom_15, wrong_strand_chrom_15, chrom_15_junction]:
            self.assertEqual(junc.chrom, "15")

        for strand, stranded in itertools.product(("+", "-", None), (True, False)):
            transcript.unfinalize()
            transcript.strand = strand
            transcript.finalize()
            sup = Superlocus(transcript, stranded=stranded)
            self.assertTrue((chrom_15_junction.junction_start, chrom_15_junction.end) in
                            sup.introns, (chrom_15_junction, sup.introns))
            sup.session = session
            asyncio.run(sup._load_introns())
            if stranded is True and strand is not None:
                self.assertEqual(sup.locus_verified_introns, {(chrom_15_junction.junction_start,
                                                               chrom_15_junction.junction_end,
                                                               strand)},
                                 (stranded, strand))
            elif stranded is False:
                self.assertEqual(sup.locus_verified_introns, {(chrom_15_junction.junction_start,
                                                               chrom_15_junction.junction_end,
                                                               chrom_15_junction.strand),
                                                              (wrong_strand_chrom_15.junction_start,
                                                               wrong_strand_chrom_15.junction_end,
                                                               wrong_strand_chrom_15.strand)},
                                 (stranded, strand))
            elif stranded is True and strand is None:
                self.assertEqual(sup.locus_verified_introns, set())


class AsyncOrfLoading(unittest.TestCase):

    def test_load_orfs(self):

        transcript_line = 'Chr1\t100\t2000\tID=foo;coding=True;phase=0'\
                          '\t0\t+\t300\t1850\t0\t4\t400,400,400,200\t0,500,1100,1700'
        transcript = Transcript(transcript_line)
        orf = transcript.orfs[0].to_transcriptomic()
        transcript2 = transcript.copy()
        transcript2.unfinalize()
        transcript2.chrom = "Chr2"
        transcript2.id = "foo.2"
        transcript2.finalize()
        other_orf = transcript2.orfs[0].to_transcriptomic()
        engine = create_engine("sqlite:///:memory:")
        db.metadata.create_all(engine)
        SessionMaker = sessionmaker(bind=engine)
        session = SessionMaker()
        query = Query(transcript.id, transcript.cdna_length)
        query2 = Query(transcript2.id, transcript2.cdna_length)
        session.add_all([query, query2])
        session.commit()
        serialized_orf = Orf(orf, query.query_id)
        self.assertEqual(serialized_orf.thick_end, orf.thick_end)
        self.assertEqual(serialized_orf.cds_len, orf.cds_len)
        serialized_other_orf = Orf(other_orf, query2.query_id)
        session.add_all([serialized_orf, serialized_other_orf])
        session.commit()
        sup = Superlocus(transcript)
        sup.session = session
        sup_orfs = asyncio.run(sup.get_orfs([query.query_id]))
        self.assertEqual(len(sup_orfs), 1)
        self.assertIn(transcript.id, sup_orfs)
        self.assertEqual(len(sup_orfs[transcript.id]), 1)
        self.assertIsInstance(sup_orfs[transcript.id][0], BED12, type(sup_orfs[transcript.id][0]))
        self.assertTrue(sup_orfs[transcript.id][0] == orf, "\n" + "\n".join(
            [str(orf), str(sup_orfs[transcript.id][0])]))


# TODO: Create a test for the BLAST hits/hsps

class AsyncBlastTest(unittest.TestCase):
    """Test for the functionality of loading a BLAST hit from a Superlocus object."""
