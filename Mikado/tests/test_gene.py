from Mikado.loci import Transcript, Gene
import unittest


class GeneTester(unittest.TestCase):

    def test_null_creation(self):

        null_gene = Gene(None)
        self.assertIsNone(null_gene.chrom)
        self.assertIsNone(null_gene.strand)
        self.assertIsNone(null_gene.source)
        self.assertIsNone(null_gene.start)
        self.assertEqual(null_gene.feature, "gene")

    def test_creation_from_transcript(self):

        t = Transcript()
        t.chrom = "Chr1"
        t.start = 100
        t.end = 200
        t.strand = "+"
        t.id = "test"
        t.parent = "parent"
        gene = Gene(t)
        self.assertIn(t.id, gene, gene.keys())
        self.assertIn(t, gene, gene.keys())
        self.assertEqual(t.chrom, gene.chrom)
        self.assertEqual(t.start, gene.start)
        self.assertEqual(t.end, gene.end)
        self.assertEqual(t.strand, gene.strand)
        self.assertIs(t, gene[t.id])

if __name__ == '__main__':
    unittest.main()