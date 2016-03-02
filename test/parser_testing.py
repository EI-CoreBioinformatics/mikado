import unittest
import Mikado.parsers
import tempfile
import os

__author__ = 'Luca Venturini'


class TestParser(unittest.TestCase):

    def test_invalid_file(self):

        with self.assertRaises(FileNotFoundError):
            _ = Mikado.parsers.GFF.GFF3("foo")

        with self.assertRaises(FileNotFoundError):
            _ = Mikado.parsers.GTF.GTF("foo")

    def test_with_construct(self):

        gff_line = "Chr1\tmikado\ttranscript\t1000\t2000\t.\t+\t.\tID=foo.1.1; Parent=foo.1;"
        gtf_line = "Chr1\tmikado\ttranscript\t1000\t2000\t.\t+\t.\tgene_id \"foo.1\"; transcript_id \"foo.1.1\";"

        gtf_temp = tempfile.NamedTemporaryFile("wt", suffix=".gtf", delete=False)
        print(gtf_line, file=gtf_temp)
        gtf_temp.flush()

        gtf_temp_reader = open(gtf_temp.name, "rt")
        gtf_temp_reader.close()
        with self.assertRaises(ValueError):
            with Mikado.parsers.GTF.GTF(gtf_temp_reader) as gtf_reader:
                _ = next(gtf_reader)

        with self.assertRaises(ValueError):
            with Mikado.parsers.GFF.GFF3(gtf_temp_reader) as gtf_reader:
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

        with Mikado.parsers.GTF.GTF(open(gtf_temp.name)) as gtf_reader:
            self.assertEqual(gtf_temp.name, gtf_reader.name)
            self.assertEqual(next(gtf_reader)._line, gtf_line)
        self.assertTrue(gtf_reader.closed)

        with Mikado.parsers.GFF.GFF3(open(gff_temp.name)) as gff_reader:
            self.assertEqual(gff_temp.name, gff_reader.name)
            self.assertEqual(next(gff_reader)._line, gff_line)
        self.assertTrue(gff_reader.closed)
        gtf_temp.close()
        gff_temp.close()

        with self.assertRaises(TypeError):
            gtf_reader.closed = "foo"  # not a boolean

        with self.assertRaises(TypeError):
            gff_reader.closed = "foo"  # not a boolean


if __name__ == '__main__':
    unittest.main()
