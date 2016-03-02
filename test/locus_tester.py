# coding: utf-8

"""
Very basic, all too basic test for some functionalities of locus-like classes.
"""

import unittest
import os.path

from Mikado.configuration import configurator
from Mikado import exceptions
from Mikado.parsers import GFF  # ,GTF, bed12
from Mikado.loci_objects import transcript, superlocus,  abstractlocus
from Mikado.utilities.log_utils import create_null_logger


class OverlapTester(unittest.TestCase):

    def test_overlap(self):
        """
        Test for overlap function
        :return:
        """

        self.assertEqual(abstractlocus.Abstractlocus.overlap((100, 200), (100, 200)),
                         100)
        self.assertEqual(abstractlocus.Abstractlocus.overlap((100, 200), (1000, 2001)),
                         -800)
        self.assertEqual(abstractlocus.Abstractlocus.overlap((100, 200), (200,300)),
                         0)


class LocusTester(unittest.TestCase):

    def setUp(self):

        gff_transcript1 = """Chr1\tfoo\ttranscript\t101\t300\t.\t+\t.\tID=t0
Chr1\tfoo\texon\t101\t300\t.\t+\t.\tID=t0:exon1;Parent=t0""".split("\n")
        gff_transcript1 = [GFF.GffLine(x) for x in gff_transcript1]
        self.assertEqual(gff_transcript1[0].chrom, "Chr1", gff_transcript1[0])
        self.transcript1 = transcript.Transcript(gff_transcript1[0])
        for exon in gff_transcript1[1:]:
            self.transcript1.add_exon(exon)
        self.transcript1.finalize()
        self.assertTrue(self.transcript1.monoexonic)
        self.assertEqual(self.transcript1.chrom, gff_transcript1[0].chrom)

        gff_transcript2 = """Chr1\tfoo\ttranscript\t101\t600\t.\t+\t.\tID=t1
Chr1\tfoo\texon\t101\t200\t.\t+\t.\tID=t1:exon1;Parent=t1
Chr1\tfoo\texon\t301\t400\t.\t+\t.\tID=t1:exon2;Parent=t1
Chr1\tfoo\texon\t501\t600\t.\t+\t.\tID=t1:exon3;Parent=t1""".split("\n")
        gff_transcript2 = [GFF.GffLine(x) for x in gff_transcript2]
        self.transcript2 = transcript.Transcript(gff_transcript2[0])
        for exon in gff_transcript2[1:-1]:
            self.transcript2.add_exon(exon)
        # Test that a transcript cannot be finalized if the exons do not define the external boundaries
        with self.assertRaises(exceptions.InvalidTranscript):
            self.transcript2.finalize()
        self.transcript2.add_exon(gff_transcript2[-1])
        self.transcript2.finalize()
        self.assertFalse(self.transcript2.monoexonic)
        self.assertEqual(self.transcript2.exon_num, len(gff_transcript2) - 1)
        # Test that trying to modify a transcript after it has been finalized causes errors
        with self.assertRaises(exceptions.ModificationError):
            for exon in gff_transcript2[1:]:
                self.transcript2.add_exon(exon)
        # Test that creating a superlocus without configuration fails
        with self.assertRaises(exceptions.NoJsonConfigError):
            _ = superlocus.Superlocus(self.transcript1)

    def test_locus(self):
        """Basic testing of the Locus functionality."""

        my_json = os.path.join(os.path.dirname(__file__), "configuration.yaml")
        my_json = configurator.to_json(my_json)
        logger = create_null_logger("null")
        logger.setLevel("WARNING")
        logger.info("Started")
        slocus = superlocus.Superlocus(self.transcript1, json_conf=my_json,
                                       logger=logger)
        slocus.add_transcript_to_locus(self.transcript2)
        self.assertEqual(slocus.strand, self.transcript1.strand)
        self.assertEqual(slocus.start, min(self.transcript1.start, self.transcript2.start))
        self.assertEqual(slocus.end, max(self.transcript1.end, self.transcript2.end))
        logger.info(slocus.transcripts)
        slocus.define_subloci()
        logger.info(slocus.subloci)
        logger.info(slocus.transcripts)
        self.assertEqual(len(slocus.transcripts), 2)
        self.assertEqual(len(slocus.subloci), 2)
        slocus.define_monosubloci()
        self.assertEqual(len(slocus.monosubloci), 2)
        slocus.define_loci()
        self.assertEqual(len(slocus.loci), 1)
        self.assertEqual(list(slocus.loci[list(slocus.loci.keys())[0]].transcripts.keys())[0], "t1")
        gff_transcript3 = """Chr1\tfoo\ttranscript\t101\t200\t.\t-\t.\tID=tminus0
Chr1\tfoo\texon\t101\t200\t.\t-\t.\tID=tminus0:exon1;Parent=tminus0""".split("\n")
        gff_transcript3 = [GFF.GffLine(x) for x in gff_transcript3]
        transcript3 = transcript.Transcript(gff_transcript3[0])
        for exon in gff_transcript3[1:]:
                transcript3.add_exon(exon)
        transcript3.finalize()
        minusuperlocus = superlocus.Superlocus(transcript3, json_conf=my_json)
        minusuperlocus.define_loci()
        self.assertEqual(len(minusuperlocus.loci), 1)
        self.assertTrue(transcript3.strand != self.transcript1.strand)

if __name__ == '__main__':
    unittest.main()
