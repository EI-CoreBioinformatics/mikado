import Mikado.utilities
import Mikado.exceptions
from Mikado.parsers.GFF import GffLine
from Mikado.loci import Transcript
from Mikado.subprograms.util.trim import trim_coding, trim_noncoding
import unittest

__author__ = 'Luca Venturini'


class TestTrimming(unittest.TestCase):

    def test_wrong_cds(self):

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
        cds_line = "\t".join(["15",
                              "protein_coding",
                              "CDS",
                              "48051996",
                              "48051996",
                              ".",
                              "+",
                              "0",
                              "ID=ENST00000560636.cds1;Parent=ENST00000560636"])
        cds_line = GffLine(cds_line)
        transcript.add_exon(cds_line)
        logger = Mikado.utilities.log_utils.create_null_logger()
        transcript.logger = logger
        with self.assertLogs("null", level="WARNING"):
            transcript.finalize()

        trimmed = trim_coding(transcript, logger, max_length=50)
        self.assertEqual(trimmed.start, 47631366)
        self.assertEqual(trimmed.end, 48051992)

    # @unittest.SkipTest
    def test_correct_cds(self):

        transcript = Transcript()
        transcript.chrom = "Chr1"
        transcript.source = "test"
        transcript.start = 10000
        transcript.end = 20000

        exons = [(10000, 11500),
                 (12000, 13000),
                 (15000, 18000),
                 (19000, 20000)]

        cds = [(11400, 11500),  # 101
               (12000, 13000),  # 1001 ==> 1102
               (15000, 17998)]  # 2998 == > 3090 (y)

        transcript.add_exons(exons)
        transcript.add_exons(cds, features="CDS")

        transcript.strand = "+"
        transcript.finalize()

        logger = Mikado.utilities.log_utils.create_null_logger("correct_cds")

        copied = transcript.deepcopy()

        trimmed = trim_coding(copied, logger, max_length=50)
        self.assertEqual(trimmed.start, 11400)
        self.assertEqual(trimmed.end, 19050)

        copied = transcript.deepcopy()
        self.assertEqual(copied.start, 10000)
        trimmed = trim_coding(copied, logger, max_length=200)
        self.assertEqual(trimmed.start, 11300)
        self.assertEqual(trimmed.end, 19200)

    def test_noncoding(self):

        logger = Mikado.utilities.log_utils.create_null_logger("correct_cds2")
        transcript = Transcript(logger=logger)
        transcript.chrom = "Chr1"
        transcript.source = "test"
        transcript.start = 10000
        transcript.end = 20000

        exons = [(10000, 11500),
                 (12000, 13000),
                 (15000, 18000),
                 (19000, 20000)]

        transcript.add_exons(exons)
        transcript.strand = "+"
        transcript.finalize()

        copied = transcript.deepcopy()

        trimmed = trim_noncoding(copied, max_length=50)
        self.assertEqual(trimmed.start, 11450)
        self.assertEqual(trimmed.end, 19050)

        copied = transcript.deepcopy()

        trimmed = trim_noncoding(copied, max_length=200)
        self.assertEqual(trimmed.start, 11300)
        self.assertEqual(trimmed.end, 19200)

if __name__ == "__main__":
    unittest.main()
