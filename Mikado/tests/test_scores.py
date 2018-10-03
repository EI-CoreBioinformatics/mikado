import Mikado
from Mikado.transcripts import Transcript
from Mikado.parsers.bed12 import BED12
from Mikado.utilities.log_utils import create_default_logger
import unittest
import ujson as json
import io
from pkg_resources import resource_stream


with io.TextIOWrapper(resource_stream("Mikado.configuration",
                                      "requirements_blueprint.json")) as rs_blueprint:
    require_schema = json.load(rs_blueprint)


class ScoreTester(unittest.TestCase):

    def setUp(self):

        b1 = BED12("1\t100\t900\tID=t1;coding=True\t0\t+\t200\t800\t.\t1\t800\t0", coding=True, transcriptomic=False)

        self.t1 = Transcript(b1)
        self.t1.finalize()
        self.assertTrue(self.t1.is_coding)

        b2 = BED12("1\t100\t1100\tID=t2;coding=True\t0\t+\t200\t900\t.\t2\t300,300\t0,700", coding=True)
        self.t2 = Transcript(b2)
        self.t2.finalize()
        self.assertTrue(self.t2.is_coding)

        # 200, 400  - 200
        # 650, 850  - 300
        # 1200, 1360    - 160

        b3 = BED12("1\t100\t1500\tID=t3;coding=True\t0\t+\t200\t1360\t.\t3\t300,300,300\t0,550,1100", coding=True)
        self.assertFalse(b3.header)
        self.assertEqual(b3.thick_start, 201)
        self.assertFalse(b3.invalid, b3.invalid_reason)
        self.assertEqual(b3.block_count, 3)
        self.assertEqual(b3.block_sizes, [300, 300, 300])
        self.t3 = Transcript(b3)
        self.t3.finalize()
        self.assertTrue(self.t3.is_coding)
        json_conf = Mikado.loci.abstractlocus.json_conf

        self.locus = Mikado.loci.Superlocus(self.t1, json_conf=json_conf)
        self.locus.add_transcript_to_locus(self.t2)
        self.locus.add_transcript_to_locus(self.t3)
        print(self.locus.json_conf["requirements"])
        self.locus.json_conf["requirements"] = {}

        reqs = {"requirements":
             {"expression": ["cdna_length"],
              "parameters": {
                  "cdna_length": {"operator": "gt", "value": 0}
                }
              }
         }

        reqs = Mikado.configuration.configurator.check_requirements(reqs, require_schema, "requirements")

        self.locus.json_conf["requirements"] = reqs["requirements"]

    def test_exon_num_max(self):

        for multiplier in (1, 2, 3):
            with self.subTest(multiplier=multiplier):
                scoring = {"exon_num": {"rescaling": "max", "use_raw": False, "multiplier": multiplier}}
                logger = create_default_logger("test_exon_num_max", level="WARNING")
                self.locus.json_conf["scoring"] = scoring
                self.assertEqual(self.locus.json_conf["scoring"]["exon_num"]["multiplier"], multiplier)
                self.assertIn("t3", self.locus.transcripts)
                self.locus.logger = logger
                self.locus.calculate_scores()
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
                self.locus.calculate_scores()
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
                self.locus.calculate_scores()
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
                self.locus.calculate_scores()
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
                self.locus.calculate_scores()
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
                self.locus.calculate_scores()
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
                self.locus.calculate_scores()

                self.assertIn("t3", self.locus.transcripts)
                self.assertIn("t3", self.locus.scores, self.locus.scores)
                self.assertEqual(self.locus.scores["t3"]["combined_cds_length"], 0, self.locus.scores)
                self.assertEqual(self.locus.scores["t2"]["combined_cds_length"], 0, self.locus.scores)
                self.assertEqual(self.locus.scores["t1"]["combined_cds_length"], multiplier)
                self.locus.scores_calculated = False