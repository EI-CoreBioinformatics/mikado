# import Mikado
from .. import loci, configuration, transcripts
from ..transcripts import Transcript
from ..parsers.bed12 import BED12
from ..utilities.log_utils import create_default_logger
import unittest
import rapidjson as json
import io
from pkg_resources import resource_stream


with io.TextIOWrapper(resource_stream("Mikado.configuration",
                                      "requirements_blueprint.json")) as rs_blueprint:
    require_schema = json.loads(rs_blueprint.read())


class ScoreTester(unittest.TestCase):

    def setUp(self):

        b1 = BED12("1\t100\t900\tID=t1;coding=True\t10\t+\t200\t800\t.\t1\t800\t0", coding=True, transcriptomic=False)

        self.t1 = Transcript(b1)
        self.t1.finalize()
        self.assertTrue(self.t1.is_coding)

        b2 = BED12("1\t100\t1100\tID=t2;coding=True\t20\t+\t200\t900\t.\t2\t300,300\t0,700", coding=True)
        self.t2 = Transcript(b2)
        self.t2.finalize()
        self.assertTrue(self.t2.is_coding)

        # 200, 400  - 200
        # 650, 850  - 300
        # 1200, 1360    - 160

        b3 = BED12("1\t100\t1500\tID=t3;coding=True\t30\t+\t200\t1360\t.\t3\t300,300,300\t0,550,1100", coding=True)
        self.assertFalse(b3.header)
        self.assertEqual(b3.thick_start, 201)
        self.assertFalse(b3.invalid, b3.invalid_reason)
        self.assertEqual(b3.block_count, 3)
        self.assertTrue((b3.block_sizes == [300, 300, 300]).all())
        self.t3 = Transcript(b3)
        self.t3.finalize()
        self.assertTrue(self.t3.is_coding)
        self.json_conf = loci.abstractlocus.json_conf
        reqs = {"requirements":
                    {"expression": ["cdna_length"],
                     "parameters": {
                         "cdna_length": {"operator": "gt", "value": 0}
                     }
                     }
                }

        reqs = configuration.configurator.check_requirements(reqs, require_schema, "requirements")
        self.json_conf["requirements"] = reqs["requirements"]
        self.locus = loci.Superlocus(self.t1, json_conf=self.json_conf)
        self.locus.add_transcript_to_locus(self.t2)
        self.locus.add_transcript_to_locus(self.t3)

    def test_exon_num_max(self):

        for multiplier in (1, 2, 3):
            with self.subTest(multiplier=multiplier):
                scoring = {"exon_num": {"rescaling": "max", "use_raw": False, "multiplier": multiplier}}
                logger = create_default_logger("test_exon_num_max", level="WARNING")
                self.locus.json_conf["scoring"] = scoring
                self.assertEqual(self.locus.json_conf["scoring"]["exon_num"]["multiplier"], multiplier)
                self.assertIn("t3", self.locus.transcripts)
                self.locus.logger = logger
                self.locus.filter_and_calculate_scores()
                self.assertEqual(self.locus.scores["t1"]["exon_num"], 0)
                self.assertIn("t3", self.locus.transcripts)
                self.assertIn("t3", self.locus.scores, self.locus.scores)
                self.assertEqual(self.locus.scores["t3"]["exon_num"], multiplier,
                                 self.locus.json_conf["scoring"])
                self.assertEqual(self.locus.scores["t2"]["exon_num"], 0.5 * multiplier, self.locus.scores)
                self.locus.scores_calculated = False

    def test_exon_min(self):

        for multiplier in (1, 2, 3):
            with self.subTest(multiplier=multiplier):
                scoring = {"exon_num": {"rescaling": "min", "use_raw": False, "multiplier": multiplier}}

                logger = create_default_logger("test_exon_num_max", level="WARNING")

                self.locus.json_conf["scoring"] = scoring
                self.assertIn("t3", self.locus.transcripts)
                self.locus.logger = logger
                self.locus.filter_and_calculate_scores()
                self.assertEqual(self.locus.scores["t1"]["exon_num"], multiplier)
                self.assertIn("t3", self.locus.transcripts)
                self.assertIn("t3", self.locus.scores, self.locus.scores)
                self.assertEqual(self.locus.scores["t3"]["exon_num"], 0, self.locus.scores)
                self.assertEqual(self.locus.scores["t2"]["exon_num"], 0.5 * multiplier, self.locus.scores)
                self.locus.scores_calculated = False

    def test_exon_target(self):

        for multiplier in (1, 2, 3):
            with self.subTest(multiplier=multiplier):
                scoring = {"exon_num": {"rescaling": "target", "value": 2,
                           "use_raw": False, "multiplier": multiplier}}

                logger = create_default_logger("test_exon_num_max", level="WARNING")

                self.locus.json_conf["scoring"] = scoring
                self.assertIn("t3", self.locus.transcripts)
                self.locus.logger = logger
                self.locus.filter_and_calculate_scores()
                self.assertEqual(self.locus.scores["t1"]["exon_num"], 0)
                self.assertIn("t3", self.locus.transcripts)
                self.assertIn("t3", self.locus.scores, self.locus.scores)
                self.assertEqual(self.locus.scores["t3"]["exon_num"], 0, self.locus.scores)
                self.assertEqual(self.locus.scores["t2"]["exon_num"], multiplier, self.locus.scores)
                self.locus.scores_calculated = False

    def test_exon_target_min(self):

        for multiplier in (1, 2, 3):
            with self.subTest(multiplier=multiplier):
                scoring = {"exon_num": {"rescaling": "target", "value": 1,
                                        "use_raw": False, "multiplier": multiplier}}

                logger = create_default_logger("test_exon_num_max", level="WARNING")

                self.locus.json_conf["scoring"] = scoring
                self.assertIn("t3", self.locus.transcripts)
                self.locus.logger = logger
                self.locus.filter_and_calculate_scores()
                self.assertEqual(self.locus.scores["t1"]["exon_num"], multiplier)
                self.assertIn("t3", self.locus.transcripts)
                self.assertIn("t3", self.locus.scores, self.locus.scores)
                self.assertEqual(self.locus.scores["t3"]["exon_num"], 0, self.locus.scores)
                self.assertEqual(self.locus.scores["t2"]["exon_num"], 0.5 * multiplier, self.locus.scores)
                self.locus.scores_calculated = False

    def test_exon_max_filter(self):

        for multiplier in (1, 2, 3):
            with self.subTest(multiplier=multiplier):
                scoring = {"exon_num": {"rescaling": "max", "use_raw": False, "multiplier": multiplier,
                                        "filter": {"operator": "gt", "value": 2}
                                        }}

                logger = create_default_logger("test_exon_num_max", level="WARNING")

                self.locus.json_conf["scoring"] = scoring
                self.assertIn("t3", self.locus.transcripts)
                self.locus.logger = logger
                self.locus.filter_and_calculate_scores()
                self.assertEqual(self.locus.scores["t1"]["exon_num"], 0)
                self.assertIn("t3", self.locus.transcripts)
                self.assertIn("t3", self.locus.scores, self.locus.scores)
                self.assertEqual(self.locus.scores["t3"]["exon_num"], multiplier, self.locus.scores)
                self.assertEqual(self.locus.scores["t2"]["exon_num"], 0, self.locus.scores)
                self.locus.scores_calculated = False

    def test_combined_cds_length_exon_filter(self):

        for multiplier in (1, 2, 3):
            with self.subTest(multiplier=multiplier):
                scoring = {"combined_cds_length": {"rescaling": "max", "use_raw": False, "multiplier": multiplier,
                                                   "filter": {"operator": "gt", "value": 2, "metric": "exon_num"}
                                                   }}

                logger = create_default_logger("test_exon_num_max", level="WARNING")

                self.locus.json_conf["scoring"] = scoring
                self.assertIn("t3", self.locus.transcripts)
                self.locus.logger = logger
                self.locus.filter_and_calculate_scores()
                self.assertEqual(self.locus.scores["t1"]["combined_cds_length"], 0)
                self.assertIn("t3", self.locus.transcripts)
                self.assertIn("t3", self.locus.scores, self.locus.scores)
                self.assertEqual(self.locus.scores["t3"]["combined_cds_length"], multiplier, self.locus.scores)
                self.assertEqual(self.locus.scores["t2"]["combined_cds_length"], 0, self.locus.scores)
                self.locus.scores_calculated = False

    def test_combined_cds_length_exon_filter_lt2(self):

        for multiplier in (1, 2, 3):
            with self.subTest(multiplier=multiplier):
                scoring = {"combined_cds_length": {"rescaling": "max", "use_raw": False, "multiplier": multiplier,
                                                   "filter": {"operator": "lt", "value": 3, "metric": "exon_num"}
                                                   }}

                logger = create_default_logger("test_exon_num_max", level="WARNING")

                self.locus.json_conf["scoring"] = scoring
                self.assertIn("t3", self.locus.transcripts)
                self.locus.logger = logger
                self.locus.filter_and_calculate_scores()

                self.assertIn("t3", self.locus.transcripts)
                self.assertIn("t3", self.locus.scores, self.locus.scores)
                self.assertEqual(self.locus.scores["t3"]["combined_cds_length"], 0, self.locus.scores)
                self.assertEqual(self.locus.scores["t2"]["combined_cds_length"], 0, self.locus.scores)
                self.assertEqual(self.locus.scores["t1"]["combined_cds_length"], multiplier)
                self.locus.scores_calculated = False

    def test_default_scores(self):

        json_conf = loci.abstractlocus.json_conf
        locus = loci.Superlocus(self.t1, use_transcript_scores=True, json_conf=json_conf)
        locus.add_transcript_to_locus(self.t2)
        locus.add_transcript_to_locus(self.t3)
        self.assertTrue(locus._use_transcript_scores)
        self.assertTrue(locus.scores_calculated)
        for transcript in (self.t1, self.t2, self.t3):
            self.assertIn(transcript.id, locus)
            self.assertEqual(locus.scores[transcript.id], transcript.score)
        locus.filter_and_calculate_scores()
        for transcript in (self.t1, self.t2, self.t3):
            self.assertIn(transcript.id, locus)
            self.assertEqual(locus.scores[transcript.id], transcript.score)


class LocusMissedTester(unittest.TestCase):

    def setUp(self):
        self.json_conf = loci.abstractlocus.json_conf
        reqs = {"requirements":
                    {"expression": ["cdna_length"],
                     "parameters": {
                         "cdna_length": {"operator": "gt", "value": 0}
                     }
                     }
                }

        reqs = configuration.configurator.check_requirements(reqs, require_schema, "requirements")
        self.json_conf["requirements"] = reqs["requirements"]
        self.json_conf["pick"]["alternative_splicing"]["pad"] = False

    def test_transcript_not_missed(self):

        b1 = BED12("1\t100\t500\tID=t1;coding=True\t20\t+\t200\t459\t.\t2\t200,100,\t0,300", coding=True)
        self.assertFalse(b1.invalid, b1.invalid_reason)
        self.assertTrue(b1.coding)
        self.assertTrue(b1.is_transcript)
        self.assertEqual(b1.thick_start, 201)
        self.assertEqual(b1.thick_end, 459)
        logger = create_default_logger("test_transcript_missed", level="ERROR")
        t1 = transcripts.Transcript(b1, logger=logger)
        t1.finalize()
        self.assertEqual(sorted(t1.exons), [(101, 300), (401,500)])

        b2 = BED12("1\t100\t1500\tID=t2;coding=True\t25\t+\t200\t901\t.\t4\t200,150,400,200\t0,300,700,1200",
                   coding=True)
        self.assertFalse(b2.invalid, b2.invalid_reason)
        self.assertTrue(b2.coding)
        self.assertTrue(b2.is_transcript)
        self.assertEqual(b2.thick_start, 201)
        self.assertEqual(b2.thick_end, 901)
        # logger.setLevel("DEBUG")
        t2 = transcripts.Transcript(b2, logger=logger)
        t2.finalize()
        self.assertTrue(t2.is_coding)
        self.assertEqual(sorted(t2.exons), [(101, 300), (401, 550), (801, 1200), (1301, 1500)])
        self.assertEqual([t2.combined_cds_start, t2.combined_cds_end], [201, 901])

        b3 = BED12("1\t800\t1500\tID=t3;coding=True\t40\t+\t820\t1370\t.\t4\t50,100,100,200\t0,150,300,500",
                   coding=True)
        self.assertFalse(b3.invalid, b3.invalid_reason)
        self.assertTrue(b3.coding)
        self.assertTrue(b3.is_transcript)
        self.assertEqual(b3.thick_start, 821)
        self.assertEqual(b3.thick_end, 1370)
        # logger.setLevel("DEBUG")
        t3 = transcripts.Transcript(b3, logger=logger)
        t3.finalize()
        self.assertTrue(t3.is_coding)
        self.assertEqual(sorted(t3.exons), [(801, 850), (951, 1050), (1101, 1200), (1301, 1500)])
        self.assertEqual([t3.combined_cds_start, t3.combined_cds_end], [821, 1370])

        locus = loci.Superlocus(t1, use_transcript_scores=True, json_conf=self.json_conf, logger=logger)
        locus.add_transcript_to_locus(t2)
        locus.add_transcript_to_locus(t3)
        locus.define_loci()
        self.assertEqual(len(locus.loci), 2)
        primaries = set([locus.loci[_].primary_transcript_id for _ in locus.loci])
        self.assertEqual(primaries, {t3.id, t1.id})

    def test_transcript_missed(self):

        """This unit-test describes a situation which currently is pathological in Mikado (issue #131).
        Namely, given three transcripts (A, B, C), with A sharing an intron with B and having *lower* score;
        and B sharing a *splice site* with C but having lower score than the latter; Mikado will end up selecting
        B and C at the sublocus stage, then discard B - without recovering A. This will lead to potential loss of good
        loci.
        """

        b1 = BED12("1\t100\t500\tID=t1;coding=True\t20\t+\t200\t459\t.\t2\t200,100,\t0,300", coding=True)
        self.assertFalse(b1.invalid, b1.invalid_reason)
        self.assertTrue(b1.coding)
        self.assertTrue(b1.is_transcript)
        self.assertEqual(b1.thick_start, 201)
        self.assertEqual(b1.thick_end, 459)
        logger = create_default_logger("test_transcript_missed", level="ERROR")
        t1 = transcripts.Transcript(b1, logger=logger)
        t1.finalize()
        self.assertEqual(sorted(t1.exons), [(101, 300), (401, 500)])

        b2 = BED12("1\t100\t1500\tID=t2;coding=True\t25\t+\t200\t901\t.\t4\t200,150,400,170\t0,300,700,1230",
                   coding=True)
        self.assertFalse(b2.invalid, b2.invalid_reason)
        self.assertTrue(b2.coding)
        self.assertTrue(b2.is_transcript)
        self.assertEqual(b2.thick_start, 201)
        self.assertEqual(b2.thick_end, 901)
        # logger.setLevel("DEBUG")
        t2 = transcripts.Transcript(b2, logger=logger)
        t2.finalize()
        self.assertTrue(t2.is_coding)
        self.assertEqual(sorted(t2.exons), [(101, 300), (401, 550), (801, 1200), (1331, 1500)])
        self.assertEqual([t2.combined_cds_start, t2.combined_cds_end], [201, 901])

        b3 = BED12("1\t800\t1500\tID=t3;coding=True\t40\t+\t820\t1370\t.\t4\t50,100,100,200\t0,150,300,500",
                   coding=True)
        self.assertFalse(b3.invalid, b3.invalid_reason)
        self.assertTrue(b3.coding)
        self.assertTrue(b3.is_transcript)
        self.assertEqual(b3.thick_start, 821)
        self.assertEqual(b3.thick_end, 1370)
        # logger.setLevel("DEBUG")
        t3 = transcripts.Transcript(b3, logger=logger)
        t3.finalize()
        self.assertTrue(t3.is_coding)
        self.assertEqual(sorted(t3.exons), [(801, 850), (951, 1050), (1101, 1200), (1301, 1500)])
        self.assertEqual([t3.combined_cds_start, t3.combined_cds_end], [821, 1370])

        locus = loci.Superlocus(t1, use_transcript_scores=True, json_conf=self.json_conf, logger=logger)
        locus.add_transcript_to_locus(t2)
        locus.add_transcript_to_locus(t3)
        locus.define_loci()
        self.assertEqual(len(locus.loci), 2)
        primaries = set([locus.loci[_].primary_transcript_id for _ in locus.loci])
        self.assertEqual(primaries, {t3.id, t1.id})
        # self.assertEqual(set(locus.lost_transcripts.keys()), {t1.id})
