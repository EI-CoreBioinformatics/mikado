import Mikado.loci.transcript
import Mikado.utilities.log_utils
import Mikado.exceptions
import Mikado.parsers
import re
import unittest

__author__ = 'Luca Venturini'


class MultOrfTester(unittest.TestCase):

    logger = Mikado.utilities.log_utils.create_default_logger("mult_orf")

    tr_gff = """
    Chr1    TAIR10    mRNA    5928    8737    .    -    .    ID=AT1G01020.1;Parent=AT1G01020;Name=AT1G01020.1;Index=1
    Chr1    TAIR10    five_prime_UTR    8667    8737    .    -    .    Parent=AT1G01020.1
    Chr1    TAIR10    CDS    8571    8666    .    -    0    Parent=AT1G01020.1;
    Chr1    TAIR10    exon    8571    8737    .    -    .    Parent=AT1G01020.1
    Chr1    TAIR10    CDS    8417    8464    .    -    0    Parent=AT1G01020.1;
    Chr1    TAIR10    exon    8417    8464    .    -    .    Parent=AT1G01020.1
    Chr1    TAIR10    CDS    8236    8325    .    -    0    Parent=AT1G01020.1;
    Chr1    TAIR10    exon    8236    8325    .    -    .    Parent=AT1G01020.1
    Chr1    TAIR10    CDS    7942    7987    .    -    0    Parent=AT1G01020.1;
    Chr1    TAIR10    exon    7942    7987    .    -    .    Parent=AT1G01020.1
    Chr1    TAIR10    exon    7762    7835    .    -    .    Parent=AT1G01020.1
    Chr1    TAIR10    CDS    7564    7649    .    -    0    Parent=AT1G01020.1;
    Chr1    TAIR10    exon    7564    7649    .    -    .    Parent=AT1G01020.1
    Chr1    TAIR10    CDS    7384    7450    .    -    1    Parent=AT1G01020.1;
    Chr1    TAIR10    exon    7384    7450    .    -    .    Parent=AT1G01020.1
    Chr1    TAIR10    CDS    7157    7232    .    -    0    Parent=AT1G01020.1;
    Chr1    TAIR10    exon    7157    7232    .    -    .    Parent=AT1G01020.1
    Chr1    TAIR10    CDS    6915    7069    .    -    2    Parent=AT1G01020.1;
    Chr1    TAIR10    three_prime_UTR    6437    6914    .    -    .    Parent=AT1G01020.1
    Chr1    TAIR10    exon    6437    7069    .    -    .    Parent=AT1G01020.1
    Chr1    TAIR10    three_prime_UTR    5928    6263    .    -    .    Parent=AT1G01020.1
    Chr1    TAIR10    exon    5928    6263    .    -    .    Parent=AT1G01020.1"""

    tr_lines = [line.strip() for line in tr_gff.split("\n") if line.strip() != '']
    for pos, line in enumerate(tr_lines):
        tr_lines[pos] = re.sub("\s+", "\t", line)
        assert len(tr_lines[pos].split("\t")) == 9, line.split("\t")

    tr_gff_lines = [Mikado.parsers.GFF.GffLine(line) for line in tr_lines]

    for l in tr_gff_lines:
        assert l.header is False

    def setUp(self):
        """Basic creation test."""

        self.tr = Mikado.loci.transcript.Transcript(self.tr_gff_lines[0], logger=self.logger)
        for line in self.tr_gff_lines[1:]:
            self.tr.add_exon(line)

    def test_missing_cds(self):

        with self.assertLogs("mult_orf", level="WARNING") as cm:
            self.tr.finalize()
        self.assertFalse(self.tr.is_coding)

    def test_split_cds_exon_one(self):

        lines = """
        Chr1    TAIR10    CDS    7762    7800    .    -    .    Parent=AT1G01020.1
        Chr1    TAIR10    CDS    7810    7835    .    -    .    Parent=AT1G01020.1
        """

        lines = ["\t".join(line.strip().split()) for line in lines.split("\n") if line.strip() != '']
        for line in lines:
            line = Mikado.parsers.GFF.GffLine(line)
            self.assertFalse(line.header)
            self.assertEqual(line.parent[0], self.tr.id)
            self.tr.add_exon(line)
        with self.assertLogs("mult_orf", level="WARNING") as cm:
            self.tr.finalize()
        self.assertFalse(self.tr.is_coding)

    def test_split_cds_exon_two(self):

        lines = """
        Chr1    TAIR10    CDS    7762    7800    .    -    .    Parent=AT1G01020.1
        Chr1    TAIR10    CDS    7800    7835    .    -    .    Parent=AT1G01020.1
        """

        lines = ["\t".join(line.strip().split()) for line in lines.split("\n") if line.strip() != '']
        for line in lines:
            line = Mikado.parsers.GFF.GffLine(line)
            self.assertFalse(line.header)
            self.tr.add_exon(line)
        with self.assertLogs("mult_orf", level="WARNING") as cm:
            self.tr.finalize()
        self.assertFalse(self.tr.is_coding)

    def test_restore(self):
        lines = """
        Chr1    TAIR10    CDS    7762    7835    .    -    0    Parent=AT1G01020.1
        """

        lines = ["\t".join(line.strip().split()) for line in lines.split("\n") if line.strip() != '']
        for line in lines:
            line = Mikado.parsers.GFF.GffLine(line)
            self.tr.add_exon(line)
        self.assertTrue(self.tr.is_coding)


if __name__ == '__main__':
    unittest.main()
