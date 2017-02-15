import unittest

from Mikado.loci import Transcript


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
