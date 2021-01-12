import unittest
from ..loci import Transcript, Superlocus
import pkg_resources
import os
from ..configuration.configurator import to_json


class ExternalTester(unittest.TestCase):
    
    def setUp(self):
        
        self.transcript = Transcript()
        self.transcript.chrom = "15"
        self.transcript.source = "protein_coding"
        self.transcript.start = 47631264
        self.transcript.end = 48051999

        exons = [(47631264, 47631416),
                 (47704590, 47704669),
                 (47762671, 47762742),
                 (47893062, 47893093),
                 (47895572, 47895655),
                 (48051942, 48051999)]

        self.transcript.strand = "+"
        self.transcript.add_exons(exons)
        self.transcript.id = "ENST00000560636"
        self.transcript.parent = "ENSG00000137872"
        
    def test_copying(self):
        
        self.transcript.external_scores.update({"test": 0, "test1": 1})
        self.assertEqual(self.transcript.external_scores.test, 0)
        self.assertEqual(self.transcript.external_scores.test1, 1)
        transcript = self.transcript.deepcopy()
        self.assertEqual(transcript.external_scores.test, 0)
        self.assertEqual(transcript.external_scores.test1, 1)

    def test_real(self):
        conf = to_json(None)
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

        self.assertIn("scoring", conf)
        transcript.attributes["tpm"] = 10
        conf["scoring"]["attributes.tpm"] = {"rescaling": "max", "default": 0, "rtype": "float"}
        sup = Superlocus(transcript, json_conf=conf)
        sup.get_metrics()
        self.assertIn("attributes.tpm", sup._metrics[self.transcript.id])
        self.assertEqual(sup._metrics[self.transcript.id]["attributes.tpm"], 10)
