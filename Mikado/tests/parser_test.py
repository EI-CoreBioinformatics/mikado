import unittest
from .. import parsers
import tempfile
import os

__author__ = 'Luca Venturini'


class TestParser(unittest.TestCase):

    def test_invalid_file(self):

        with self.assertRaises(FileNotFoundError):
            _ = parsers.GFF.GFF3("foo")

        with self.assertRaises(FileNotFoundError):
            _ = parsers.GTF.GTF("foo")

    def test_with_construct(self):

        gff_line = "Chr1\tmikado\ttranscript\t1000\t2000\t.\t+\t.\tID=foo.1.1; Parent=foo.1;"
        gtf_line = "Chr1\tmikado\ttranscript\t1000\t2000\t.\t+\t.\tgene_id \"foo.1\"; transcript_id \"foo.1.1\";"

        gtf_temp = tempfile.NamedTemporaryFile("wt", suffix=".gtf", delete=False)
        print(gtf_line, file=gtf_temp)
        gtf_temp.flush()

        gtf_temp_reader = open(gtf_temp.name, "rt")
        gtf_temp_reader.close()
        with self.assertRaises(ValueError):
            with parsers.GTF.GTF(gtf_temp_reader) as gtf_reader:
                _ = next(gtf_reader)

        with self.assertRaises(ValueError):
            with parsers.GFF.GFF3(gtf_temp_reader) as gtf_reader:
                _ = next(gtf_reader)

        os.remove(gtf_temp.name)

    def test_name(self):
        gff_line = "Chr1\tmikado\ttranscript\t1000\t2000\t.\t+\t.\tID=foo.1.1; Parent=foo.1;"
        gtf_line = "Chr1\tmikado\ttranscript\t1000\t2000\t.\t+\t.\tgene_id \"foo.1\"; transcript_id \"foo.1.1\";"

        gtf_temp = tempfile.NamedTemporaryFile("wt", suffix=".gtf")
        gff_temp = tempfile.NamedTemporaryFile("wt", suffix=".gff3")

        print(gff_line, file=gff_temp)
        print(gtf_line, file=gtf_temp)

        gff_temp.flush()
        gtf_temp.flush()

        with parsers.GTF.GTF(open(gtf_temp.name)) as gtf_reader:
            self.assertEqual(gtf_temp.name, gtf_reader.name)
            self.assertEqual(next(gtf_reader)._line, gtf_line)
        self.assertTrue(gtf_reader.closed)

        with parsers.GFF.GFF3(open(gff_temp.name)) as gff_reader:
            self.assertEqual(gff_temp.name, gff_reader.name)
            self.assertEqual(next(gff_reader)._line, gff_line)
        self.assertTrue(gff_reader.closed)
        gtf_temp.close()
        gff_temp.close()

        with self.assertRaises(TypeError):
            gtf_reader.closed = "foo"  # not a boolean

        with self.assertRaises(TypeError):
            gff_reader.closed = "foo"  # not a boolean

    def test_phase_frame(self):
        gtf_line = "Chr1\tmikado\texon\t1000\t2000\t.\t+\t0\tgene_id \"foo.1\"; transcript_id \"foo.1.1\";"
        gtf_line = parsers.GTF.GtfLine(gtf_line)
        self.assertEqual(gtf_line.frame, 0)
        gff_line = parsers.GFF.GffLine(
            "Chr1\tmikado\texon\t1000\t2000\t.\t+\t0\tParent=foo.1.1;"
        )
        self.assertEqual(gtf_line.transcript, gff_line.transcript[0])
        self.assertEqual(gtf_line.phase, gtf_line.frame)
        self.assertEqual(gtf_line.phase, gff_line.phase)

        for val in ["foo", 3, 10, dict(), [], [10]]:
            with self.subTest(val=val, msg="Testing {}".format(val)):
                with self.assertRaises(AttributeError):
                    gtf_line.frame = val
                with self.assertRaises(ValueError):
                    gff_line.phase = val
                with self.assertRaises(ValueError):
                    gtf_line.phase = val

        gff_line.phase = -1
        self.assertEqual(gff_line.phase, 0)
        with self.assertRaises(ValueError):
            gtf_line.phase = -1

        for phase, frame in {0: 0, 1: 2, 2: 1, None: None}.items():
            gtf_line.phase = phase
            self.assertEqual(gtf_line.frame, frame)
            gff_line.phase = phase

    def test_gene_property(self):

        gff_line = "Chr1\tmikado\tgene\t1000\t2000\t.\t+\t0\tID=\"transcript:foo.1\";"
        gff_line = parsers.GFF.GffLine(gff_line)
        self.assertFalse(gff_line.is_gene)
        gff_line.id = "gene:foo.1"
        self.assertTrue(gff_line.is_gene)
        gff_line.feature = "gene_ens"
        self.assertTrue(gff_line.is_gene)
        gff_line.id = "foo.1"
        self.assertFalse(gff_line.is_gene)

    def test_cds_property(self):

        gff_line = "Chr1\tmikado\tgene\t1000\t2000\t.\t+\t0\tID=\"transcript:foo.1\";"
        gff_line = parsers.GFF.GffLine(gff_line)
        self.assertFalse(gff_line.is_cds)
        self.assertFalse(gff_line.is_exon)

        gff_line.header = True
        self.assertFalse(gff_line.is_cds)
        self.assertFalse(gff_line.is_exon)

        gff_line.header = False
        gff_line.feature = "CDS"
        gff_line.parent = gff_line.id
        self.assertTrue(gff_line.is_cds)

    def test_gff_order(self):

        gff_line = "Chr1\tmikado\texon\t1000\t2000\t.\t+\t0\tParent=foo.1.1;"
        gff_line = parsers.GFF.GffLine(gff_line)

        lines = [gff_line]

        for feature in ["five_prime_utr", "three_prime_utr", "CDS", "start_codon", "stop_codon"]:
            l = lines[0].copy()
            l.feature = feature
            lines.append(l)

        __positive_order = ["five_prime_utr",
                            "exon",
                            "start_codon",
                            "CDS",
                            "stop_codon",
                            "three_prime_utr"]

        for pos, feat in enumerate(__positive_order):
            self.assertEqual(gff_line._sort_feature(feat), pos, feat)

        self.assertEqual([_.feature for _ in sorted(lines)],
                         __positive_order)

        for _ in lines:
            _.strand = "-"

        self.assertEqual([_.feature for _ in sorted(lines)],
                         ["three_prime_utr", "exon", "stop_codon", "CDS", "start_codon", "five_prime_utr"])

        for _ in lines:
            _.strand = None

        self.assertEqual([_.feature for _ in sorted(lines)],
                         __positive_order)

        for _ in lines:
            _.strand = "+"

        l = lines[0].copy()
        l.start += 10000
        l.end += 10000
        lines.append(l)

        self.assertEqual([_.feature for _ in sorted(lines)],
                         __positive_order + ["exon"])

        for _ in lines:
            _.strand = "-"

        self.assertEqual([_.feature for _ in sorted(lines)],
                         ["three_prime_utr", "exon", "stop_codon", "CDS", "start_codon", "five_prime_utr", "exon"])

        self.assertFalse(lines[0] == lines[1])
        self.assertEqual(lines[1], lines[1])

        l = lines[0].copy()
        l.chrom = "Chr2"
        self.assertFalse(lines[0] == l)
        self.assertGreater(l, lines[0])
        l2 = lines[0]
        l2.start += 100000
        l2.end += 100000
        self.assertGreater(l, lines[0])

    def test_gtf_order(self):

        gtf_line = "Chr1\tmikado\texon\t1000\t2000\t.\t+\t0\tgene_id \"foo.1\"; transcript_id \"foo.1.1\";"
        gtf_line = parsers.GTF.GtfLine(gtf_line)

        lines = [gtf_line]
        for feature in ["5UTR", "3UTR", "CDS", "start_codon", "stop_codon"]:
            l = lines[0].copy()
            l.feature = feature
            lines.append(l)

        __positive_order = ["5UTR",
                            "exon",
                            "start_codon",
                            "CDS",
                            "stop_codon",
                            "3UTR"]

        for pos, feat in enumerate(__positive_order):
            self.assertEqual(gtf_line._sort_feature(feat), pos, feat)

        self.assertEqual([_.feature for _ in sorted(lines)],
                         ["5UTR", "exon", "start_codon", "CDS", "stop_codon", "3UTR"]
                         )

        for _ in lines:
            _.strand = "-"

        self.assertEqual([_.feature for _ in sorted(lines)],
                         ["3UTR", "exon", "stop_codon", "CDS", "start_codon", "5UTR"])

        for _ in lines:
            _.strand = None

        self.assertEqual([_.feature for _ in sorted(lines)],
                         ["5UTR", "exon", "start_codon", "CDS", "stop_codon", "3UTR"]
                         )

        for _ in lines:
            _.strand = "+"

        l = lines[0].copy()
        l.start += 10000
        l.end += 10000
        lines.append(l)

        self.assertEqual([_.feature for _ in sorted(lines)],
                         ["5UTR", "exon", "start_codon", "CDS", "stop_codon", "3UTR", "exon"])

        for _ in lines:
            _.strand = "-"

        self.assertEqual([_.feature for _ in sorted(lines)],
                         ["3UTR", "exon", "stop_codon", "CDS", "start_codon", "5UTR", "exon"])

        self.assertFalse(lines[0] == lines[1])
        self.assertEqual(lines[1], lines[1])

        l = lines[0].copy()
        l.chrom = "Chr2"
        self.assertFalse(lines[0] == l)
        self.assertGreater(l, lines[0])
        l2 = lines[0].copy()
        l2.start += 100000
        l2.end += 100000
        self.assertGreater(l, lines[0])

        for _ in lines:
            _.strand = "+"
        lines = sorted(lines)
        self.assertEqual(lines[0].start, 1000)
        lines = [_ for _ in lines if _.start == 1000]
        l = lines[0].copy()
        l.feature = "transcript"
        lines.append(l)
        self.assertEqual([_.feature for _ in sorted(lines)],
                         ["5UTR", "exon", "start_codon", "CDS", "stop_codon", "3UTR", "transcript"])

        for _ in lines:
            _.strand = "-"

        self.assertEqual([_.feature for _ in sorted(lines)],
                         ["3UTR", "exon", "stop_codon", "CDS", "start_codon", "5UTR", "transcript"])

    def test_length(self):

        gtf_line = "Chr1\tmikado\texon\t1000\t2000\t.\t+\t0\tgene_id \"foo.1\"; transcript_id \"foo.1.1\";"
        gtf_line = parsers.GTF.GtfLine(gtf_line)

        self.assertEqual(len(gtf_line), 2000-1000+1)

        gff_line = "Chr1\tmikado\texon\t1000\t2000\t.\t+\t0\tParent=foo.1.1;"
        gff_line = parsers.GFF.GffLine(gff_line)

        self.assertEqual(len(gtf_line), 2000 - 1000 + 1)

        gff_line.header, gff_line.start, gff_line.end = True, None, gff_line.end
        self.assertEqual(len(gff_line), 0)
        gtf_line.header, gtf_line.start, gtf_line.end = True, None, gtf_line.end
        self.assertEqual(len(gff_line), 0)


if __name__ == '__main__':
    unittest.main()
