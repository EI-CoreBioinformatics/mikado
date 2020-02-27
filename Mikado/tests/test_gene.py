from ..loci import Transcript, Gene
import unittest
from .. import parsers
from pkg_resources import resource_filename


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

    def test_pesky_gtf(self):
        fname = resource_filename("Mikado.tests", "gtf_check.gtf.gz")
        with parsers.GTF.GTF(fname) as gtf:
            gene = None
            for gline in gtf:
                self.assertFalse(gline.header)
                self.assertEqual(gline.gene, "ENSG00000142621")
                if gline.is_gene:
                    gene = Gene(gline)
                elif gline.is_transcript:
                    gene.add(Transcript(gline))
                else:
                    gene.add_exon(gline)

        gene.finalize()
        self.assertTrue(gene.is_coding)
        self.assertTrue(sorted(list(gene.transcripts.keys())),
                    sorted("""ENST00000314668
ENST00000358897
ENST00000375995
ENST00000375996
ENST00000375997
ENST00000375998
ENST00000433640
ENST00000444385
ENST00000459961
ENST00000462113
ENST00000471347
ENST00000472086
ENST00000472131
ENST00000472449
ENST00000477846
ENST00000481324
ENST00000483120
ENST00000483694
ENST00000495195
ENST00000524761
ENST00000529606
ENST00000532408
""".split("\n")))


if __name__ == '__main__':
    unittest.main()
