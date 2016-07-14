# coding: utf-8

"""
Very basic, all too basic test for some functionalities of locus-like classes.
"""

import unittest
import os.path
import logging
from Mikado.configuration import configurator
from Mikado import exceptions
from Mikado.parsers import GFF  # ,GTF, bed12
from Mikado.parsers.GTF import GtfLine
from Mikado.loci import Transcript, Superlocus, Abstractlocus, Locus, MonosublocusHolder
from Mikado.utilities.log_utils import create_null_logger, create_default_logger
from Mikado.utilities import overlap
import Mikado.loci


class OverlapTester(unittest.TestCase):

    def test_overlap(self):
        """
        Test for overlap function
        :return:
        """

        self.assertEqual(Abstractlocus.overlap((100, 200), (100, 200)),
                         100)
        self.assertEqual(Abstractlocus.overlap((100, 200), (100, 200)),
                         overlap((100, 200), (100, 200)))


class LocusTester(unittest.TestCase):

    logger = create_null_logger("locus_tester")

    def setUp(self):

        gff_transcript1 = """Chr1\tfoo\ttranscript\t101\t300\t.\t+\t.\tID=t0
Chr1\tfoo\texon\t101\t300\t.\t+\t.\tID=t0:exon1;Parent=t0
Chr1\tfoo\tCDS\t101\t250\t.\t+\t.\tID=t0:exon1;Parent=t0""".split("\n")
        gff_transcript1 = [GFF.GffLine(x) for x in gff_transcript1]
        self.assertEqual(gff_transcript1[0].chrom, "Chr1", gff_transcript1[0])
        self.transcript1 = Transcript(gff_transcript1[0])
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
        self.transcript2 = Transcript(gff_transcript2[0], logger=self.logger)

        for exon in gff_transcript2[1:-1]:
            self.transcript2.add_exon(exon)
        # Test that a transcript cannot be finalized if
        # the exons do not define the external boundaries
        with self.assertLogs("null", level="WARNING") as _:
            self.transcript2.finalize()
        with self.assertRaises(exceptions.ModificationError):
            self.transcript2.add_exon(gff_transcript2[-1])

        self.transcript2.finalized = False
        self.transcript2.start = 101
        self.transcript2.end = 600
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
            _ = Superlocus(self.transcript1)

    def test_locus(self):
        """Basic testing of the Locus functionality."""

        my_json = os.path.join(os.path.dirname(__file__), "configuration.yaml")
        my_json = configurator.to_json(my_json)
        logger = create_null_logger("null")
        logger.setLevel("WARNING")
        logger.info("Started")
        slocus = Superlocus(self.transcript1,
                            json_conf=my_json,
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
        self.assertEqual(list(slocus.loci[
                                  list(slocus.loci.keys())[0]].transcripts.keys())[0], "t0")
        gff_transcript3 = """Chr1\tfoo\ttranscript\t101\t200\t.\t-\t.\tID=tminus0
Chr1\tfoo\texon\t101\t200\t.\t-\t.\tID=tminus0:exon1;Parent=tminus0""".split("\n")
        gff_transcript3 = [GFF.GffLine(x) for x in gff_transcript3]
        transcript3 = Transcript(gff_transcript3[0])
        for exon in gff_transcript3[1:]:
                transcript3.add_exon(exon)
        transcript3.finalize()
        minusuperlocus = Superlocus(transcript3, json_conf=my_json)
        minusuperlocus.define_loci()
        self.assertEqual(len(minusuperlocus.loci), 1)
        self.assertTrue(transcript3.strand != self.transcript1.strand)


class ASeventsTester(unittest.TestCase):

    logger = create_null_logger("ASevents")

    def setUp(self):
        
        self.conf = dict()
        self.conf["pick"] = dict()
        self.conf["pick"]["alternative_splicing"] = dict()
        self.conf["pick"]["alternative_splicing"]["max_utr_length"] = 10000
        self.conf["pick"]["alternative_splicing"]["max_fiveutr_length"] = 10000
        self.conf["pick"]["alternative_splicing"]["max_threeutr_length"] = 10000
        self.conf["pick"]["alternative_splicing"]["valid_ccodes"] = ["j", "J", "O", "mo"]
        self.conf["pick"]["alternative_splicing"]["redundant_ccodes"] = ["c", "=", "_", "m"]
        self.conf["pick"]["alternative_splicing"]["only_confirmed_introns"] = False
        self.conf["pick"]["alternative_splicing"]["min_score_perc"] = 0.5
        self.conf["pick"]["alternative_splicing"]["keep_retained_introns"] = True
        self.conf["pick"]["alternative_splicing"]["min_cdna_overlap"] = 0.2
        self.conf["pick"]["alternative_splicing"]["min_cds_overlap"] = 0.2
        self.conf["pick"]["alternative_splicing"]["max_isoforms"] = 3        
    
        self.t1 = Transcript()
        self.t1.chrom = "Chr1"
        self.t1.strand = "+"
        self.t1.score = 20
        self.t1.id = "G1.1"
        self.t1.parent = "G1"
        self.t1.start = 101
        self.t1.end = 1500
        
        self.t1.add_exons([(101, 500), (601, 700), (1001, 1300), (1401, 1500)],
                          "exon")
        self.t1.add_exons([(401, 500), (601, 700), (1001, 1300), (1401, 1440)],
                          "CDS")
        self.t1.finalize()
        
        self.locus = Locus(self.t1)
        self.locus.logger = self.logger
        self.locus.json_conf = self.conf

    def test_not_intersecting(self):

        # This one is contained and should be rejected
        t2 = Transcript()
        t2.chrom = "Chr1"
        t2.strand = "+"
        t2.score = 20
        t2.id = "G1.1"
        t2.parent = "G1"
        t2.start = 601
        t2.end = 1420
        t2.add_exons([(601, 700), (1001, 1300), (1401, 1420)],
                     "exon")
        t2.add_exons([(601, 700), (1001, 1300), (1401, 1420)],
                     "CDS")
        t2.finalize()

        self.assertEqual(self.locus.is_alternative_splicing(t2)[:2], (False, "c"))

    def test_valid_as(self):

        t2 = Transcript()
        t2.chrom = "Chr1"
        t2.strand = "+"
        t2.score = 20
        t2.id = "G2.1"
        t2.parent = "G2"
        t2.start = 101
        t2.end = 1600
        
        t2.add_exons([(101, 500), (601, 700), (1001, 1300), (1401, 1460), (1501, 1600)],
                     "exon")
        t2.add_exons([(401, 500), (601, 700), (1001, 1300), (1401, 1440)],
                     "CDS")
        t2.finalize()        

        self.assertEqual(self.locus.is_alternative_splicing(t2)[:2], (True, "J"))

        self.locus.add_transcript_to_locus(t2)
        self.assertEqual(len(self.locus.transcripts), 2, self.locus.transcripts)
        
    def test_redundant_as(self):

        t2 = Transcript()
        t2.chrom = "Chr1"
        t2.strand = "+"
        t2.score = 20
        t2.id = "G2.1"
        t2.parent = "G2"
        t2.start = 101
        t2.end = 1600
        
        t2.add_exons([(101, 500), (601, 700), (1001, 1300), (1401, 1460), (1501, 1600)],
                     "exon")
        t2.add_exons([(401, 500), (601, 700), (1001, 1300), (1401, 1440)],
                     "CDS")

        t2.finalize()

        self.locus.add_transcript_to_locus(t2)
        self.assertEqual(len(self.locus.transcripts), 2, self.locus.transcripts)        
        
        t3 = Transcript()
        t3.chrom = "Chr1"
        t3.strand = "+"
        t3.score = 20
        t3.id = "G3.1"
        t3.parent = "G3"
        t3.start = 201
        t3.end = 1630
        
        t3.add_exons([(201, 500), (601, 700), (1001, 1300), (1401, 1460), (1501, 1630)],
                     "exon")
        t3.add_exons([(401, 500), (601, 700), (1001, 1300), (1401, 1440)],
                     "CDS")
        t3.finalize()

        self.assertEqual(self.locus.is_alternative_splicing(t3)[:2], (False, "J"))
        self.locus.add_transcript_to_locus(t3)
        self.assertEqual(len(self.locus.transcripts), 2, self.locus.transcripts)

    def test_non_redundant_as(self):

        t2 = Transcript()
        t2.chrom = "Chr1"
        t2.strand = "+"
        t2.score = 20
        t2.id = "G2.1"
        t2.parent = "G2"
        t2.start = 101
        t2.end = 1600

        t2.add_exons([(101, 500), (601, 700), (1001, 1300), (1401, 1460), (1501, 1600)],
                     "exon")
        t2.add_exons([(401, 500), (601, 700), (1001, 1300), (1401, 1440)],
                     "CDS")
        t2.finalize()

        self.locus.add_transcript_to_locus(t2)
        self.assertEqual(len(self.locus.transcripts), 2, self.locus.transcripts)

        t3 = Transcript()
        t3.chrom = "Chr1"
        t3.strand = "+"
        t3.score = 20
        t3.id = "G3.1"
        t3.parent = "G3"
        t3.start = 201
        t3.end = 1630

        t3.add_exons([(201, 500), (601, 670), (1031, 1300), (1401, 1460), (1501, 1630)],
                     "exon")
        t3.add_exons([(401, 500), (601, 670), (1031, 1300), (1401, 1440)],
                     "CDS")
        t3.logger = self.logger
        t3.finalize()

        self.assertEqual(self.locus.is_alternative_splicing(t3)[:2], (True, "j"))
        self.locus.add_transcript_to_locus(t3)
        self.assertEqual(len(self.locus.transcripts), 3, self.locus.transcripts)

    def test_lowscore(self):

        t2 = Transcript()
        t2.chrom = "Chr1"
        t2.strand = "+"
        t2.score = 1
        t2.id = "G2.1"
        t2.parent = "G2"
        t2.start = 101
        t2.end = 1600

        t2.add_exons([(101, 500), (601, 700), (1001, 1300), (1401, 1460), (1501, 1600)],
                     "exon")
        t2.add_exons([(401, 500), (601, 700), (1001, 1300), (1401, 1440)],
                     "CDS")
        t2.finalize()

        self.locus.add_transcript_to_locus(t2)
        self.assertEqual(len(self.locus.transcripts), 1, self.locus.transcripts)


class MonoHolderTester(unittest.TestCase):

    logger = create_default_logger("MonoHolderTester")

    def setUp(self):

        self.conf = dict()

        self.t1 = Transcript()
        self.t1.chrom = "Chr1"
        self.t1.strand = "+"
        self.t1.score = 20
        self.t1.id = "G1.1"
        self.t1.parent = "G1"
        self.t1.start = 101
        self.t1.end = 1500

        self.t1.add_exons([(101, 500), (601, 700), (1001, 1300), (1401, 1500)],
                          "exon")
        self.t1.add_exons([(401, 500), (601, 700), (1001, 1300), (1401, 1440)],
                          "CDS")
        self.t1.finalize()

    def testCdsOverlap(self):

        t2 = Transcript()
        t2.chrom = "Chr1"
        t2.strand = "+"
        t2.score = 1
        t2.id = "G2.1"
        t2.parent = "G2"
        t2.start = 101
        t2.end = 1600

        t2.add_exons([(101, 500), (601, 700), (1001, 1300), (1401, 1460), (1501, 1600)],
                     "exon")
        t2.add_exons([(401, 500), (601, 700), (1001, 1300), (1401, 1440)],
                     "CDS")
        t2.finalize()

        self.assertTrue(MonosublocusHolder.is_intersecting(self.t1, t2))

    def test_intronMatch(self):

        t2 = Transcript()
        t2.chrom = "Chr1"
        t2.strand = "+"
        t2.score = 1
        t2.id = "G2.1"
        t2.parent = "G2"
        t2.start = 101
        t2.end = 1600

        t2.add_exons([(101, 500), (601, 700), (1001, 1320), (1451, 1460), (1501, 1600)],
                     "exon")
        t2.add_exons([(401, 500), (601, 700), (1001, 1320), (1451, 1460), (1501, 1510)],
                     "CDS")
        t2.finalize()

        self.assertTrue(self.t1.is_coding)
        self.assertTrue(t2.is_coding)

        self.assertTrue(MonosublocusHolder.is_intersecting(self.t1, t2, logger=self.logger))
        self.assertTrue(MonosublocusHolder.is_intersecting(self.t1, t2, cds_only=True, logger=self.logger))

    def test_intronOverlap(self):

        self.t1.strip_cds()
        t2 = Transcript()
        t2.chrom = "Chr1"
        t2.strand = "+"
        t2.score = 1
        t2.id = "G2.1"
        t2.parent = "G2"
        t2.start = 101
        t2.end = 1470
        t2.add_exons([(101, 510), (601, 700), (960, 1350), (1420, 1470)])

        t2.finalize()
        self.assertTrue(MonosublocusHolder.is_intersecting(self.t1, t2))

    def test_noIntronOverlap(self):

        self.t1.strip_cds()
        t2 = Transcript()
        t2.chrom = "Chr1"
        t2.strand = "+"
        t2.score = 1
        t2.id = "G2.1"
        t2.parent = "G2"
        t2.start = 1250
        t2.end = 2000
        t2.add_exons([(1250, 1560), (1800, 2000)])
        t2.finalize()
        self.assertFalse(MonosublocusHolder.is_intersecting(self.t1, t2))

    def test_noCDSOverlap(self):

        self.t1.strip_cds()
        self.assertEqual(self.t1.combined_cds_introns, set())
        self.t1.finalized = False
        self.t1.add_exons([(401, 500), (601, 700), (1001, 1100)],
                          "CDS")
        self.t1.finalize()

        t2 = Transcript()
        t2.logger = self.logger
        t2.chrom = "Chr1"
        t2.strand = "+"
        t2.score = 1
        t2.id = "G2.1"
        t2.parent = "G2"
        t2.start = 101
        t2.end = 1470
        t2.add_exons([(101, 510), (601, 700), (960, 1350), (1421, 1470)])
        t2.add_exons([(1201, 1350), (1421, 1450)], "CDS")
        t2.finalize()

        self.assertTrue(self.t1.is_coding)
        self.assertTrue(t2.is_coding)

        self.assertGreaterEqual(0,
                                overlap(
                                    (self.t1.combined_cds_start, self.t1.combined_cds_end),
                                    (t2.combined_cds_start, t2.combined_cds_end)),
                                [(self.t1.combined_cds_start, self.t1.combined_cds_end),
                                 (t2.combined_cds_start, t2.combined_cds_end)])

        self.assertTrue(MonosublocusHolder.is_intersecting(self.t1, t2, logger=self.logger))
        self.assertFalse(MonosublocusHolder.is_intersecting(self.t1, t2, cds_only=True, logger=self.logger))

    def test_only_CDS_overlap(self):

        t2 = Transcript()
        t2.chrom = "Chr1"
        t2.strand = "+"
        t2.score = 1
        t2.id = "G2.1"
        t2.parent = "G2"
        t2.start = 1250
        t2.end = 2000
        t2.add_exons([(1250, 1560), (1801, 2000)])
        t2.add_exons([(1401, 1560), (1801, 1850)], "CDS")
        t2.finalize()
        self.assertFalse(MonosublocusHolder.is_intersecting(self.t1, t2))

        t2.strip_cds()
        t2.finalized = False
        t2.add_exons([(1461, 1560), (1801, 1850)], "CDS")
        # No CDS overlap this time
        self.assertFalse(MonosublocusHolder.is_intersecting(self.t1, t2))

    def test_no_overlap(self):

        t2 = Transcript()
        t2.chrom = "Chr1"
        t2.strand = "+"
        t2.score = 1
        t2.id = "G2.1"
        t2.parent = "G2"
        t2.start = 1600
        t2.end = 2000
        t2.add_exons([(1600, 1700), (1801, 2000)])
        t2.add_exons([(1661, 1700), (1801, 1850)], "CDS")
        t2.finalize()
        self.assertFalse(MonosublocusHolder.is_intersecting(self.t1, t2))

    def test_same_id(self):

        t2 = Transcript()
        t2.chrom = "Chr1"
        t2.strand = "+"
        t2.score = 1
        t2.id = "G1.1"
        t2.parent = "G1"
        t2.start = 1250
        t2.end = 2000
        t2.add_exons([(1250, 1560), (1801, 2000)])
        t2.add_exons([(1401, 1560), (1801, 1850)], "CDS")
        t2.finalize()
        # This fails because they have the same ID
        self.assertFalse(MonosublocusHolder.is_intersecting(self.t1, t2))


class TestLocus(unittest.TestCase):

    """
    This unit test is focussed on the locus definition and alternative splicings.
    """

    logger = Mikado.utilities.log_utils.create_default_logger("tester")

    def setUp(self):

        """Set up for the unit test."""

        # Mock dictionary to be used for the alternative splicing checks
        self.json_conf = dict()
        self.json_conf["pick"] = dict()
        self.json_conf["pick"]["alternative_splicing"] = dict()
        self.json_conf["pick"]["alternative_splicing"]["max_utr_length"] = 2000
        self.json_conf["pick"]["alternative_splicing"]["max_fiveutr_length"] = 1000
        self.json_conf["pick"]["alternative_splicing"]["max_threeutr_length"] = 1000
        self.json_conf["pick"]["alternative_splicing"]["max_isoforms"] = 3
        self.json_conf["pick"]["alternative_splicing"]["keep_retained_introns"] = False
        self.json_conf["pick"]["alternative_splicing"]["min_cds_overlap"] = 0
        self.json_conf["pick"]["alternative_splicing"]["min_cdna_overlap"] = 0
        self.json_conf["pick"]["alternative_splicing"]["min_score_perc"] = 0.1
        self.json_conf["pick"]["alternative_splicing"]["valid_ccodes"] = ["j", "n", "O", "mo"]
        self.json_conf["pick"]["alternative_splicing"]["redundant_ccodes"] = ["c", "=", "_", "m"]
        self.json_conf["pick"]["alternative_splicing"]["only_confirmed_introns"] = False

        t1 = """Chr1\tfoo\ttranscript\t1001\t3000\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.1";
        Chr1\tfoo\texon\t1001\t1300\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.1";
        Chr1\tfoo\tCDS\t1101\t1300\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.1";
        Chr1\tfoo\texon\t1701\t2000\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.1";
        Chr1\tfoo\tCDS\t1701\t2000\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.1";
        Chr1\tfoo\texon\t2101\t2500\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.1";
        Chr1\tfoo\tCDS\t2101\t2500\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.1";
        Chr1\tfoo\texon\t2801\t3000\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.1";
        Chr1\tfoo\tCDS\t2801\t2902\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.1";
        """

        t1lines = [GtfLine(line) for line in t1.split("\n") if line]
        self.t1 = Mikado.loci.Transcript(t1lines[0])
        for exon in t1lines[1:]:
            if exon.header:
                continue
            self.t1.add_exon(exon)
        self.t1.score = 20
        self.t1.finalize()

        # Just a fragment of the best
        t1_contained = """Chr1\tfoo\ttranscript\t1001\t2000\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.2";
        Chr1\tfoo\texon\t1001\t1300\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.2";
        Chr1\tfoo\tCDS\t1101\t1300\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.2";
        Chr1\tfoo\texon\t1701\t2000\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.2";
        Chr1\tfoo\tCDS\t1701\t1902\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.2";
        """

        t1_contained_lines = [GtfLine(line) for line in t1_contained.split("\n") if line]
        self.t1_contained = Mikado.loci.Transcript(t1_contained_lines[0])
        for exon in t1_contained_lines[1:]:
            if exon.header:
                continue
            self.t1_contained.add_exon(exon)
        self.t1_contained.score = 15
        self.t1_contained.finalize()

        # Valid AS event
        t1_as = """Chr1\tfoo\ttranscript\t1001\t3000\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.3";
        Chr1\tfoo\texon\t1001\t1300\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.3";
        Chr1\tfoo\tCDS\t1101\t1300\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.3";
        Chr1\tfoo\texon\t1701\t2000\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.3";
        Chr1\tfoo\tCDS\t1701\t2000\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.3";
        Chr1\tfoo\texon\t2101\t2400\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.3";
        Chr1\tfoo\tCDS\t2101\t2400\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.3";
        Chr1\tfoo\texon\t2801\t3000\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.3";
        Chr1\tfoo\tCDS\t2801\t2900\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.3";
        """

        t1_as_lines = [GtfLine(line) for line in t1_as.split("\n") if line]
        self.t1_as = Mikado.loci.Transcript(t1_as_lines[0])
        for exon in t1_as_lines[1:]:
            if exon.header:
                continue
            self.t1_as.add_exon(exon)
        self.t1_as.score = 19
        self.t1_as.finalize()

        # Retained intron AS event
        t1_retained = """Chr1\tfoo\ttranscript\t1001\t2900\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.4";
        Chr1\tfoo\texon\t1001\t1300\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.4";
        Chr1\tfoo\tCDS\t1101\t1300\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.4";
        Chr1\tfoo\texon\t1701\t2000\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.4";
        Chr1\tfoo\tCDS\t1701\t2000\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.4";
        Chr1\tfoo\texon\t2101\t2900\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.4";
        Chr1\tfoo\tCDS\t2101\t2472\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.4";
        """

        t1_retained_lines = [GtfLine(line) for line in t1_retained.split("\n") if line]
        self.t1_retained = Mikado.loci.Transcript(t1_retained_lines[0])
        for exon in t1_retained_lines[1:]:
            if exon.header:
                continue
            self.t1_retained.add_exon(exon)
        self.t1_retained.score = 10
        self.t1_retained.finalize()

        # self.logger = logging.getLogger("tester")
        # self.handler = logging.StreamHandler()
        self.logger.setLevel(logging.WARNING)
        # self.logger.addHandler(self.handler)

    def test_validity(self):
        """
        First test of validity to ensure the CCodes are as expected.
        :return:
        """

        # The fragment should have a c assigned
        result, _ = Mikado.scales.assigner.Assigner.compare(self.t1_contained, self.t1)
        self.assertEqual(result.ccode[0], "c")

        # The valid AS should have a j assigned
        result, _ = Mikado.scales.assigner.Assigner.compare(self.t1_as, self.t1)
        self.assertEqual(result.ccode[0], "j")

        # The retained intron AS should have a j assigned
        result, _ = Mikado.scales.assigner.Assigner.compare(self.t1_retained, self.t1)
        self.assertEqual(result.ccode[0], "j", result.ccode)

    def testCreate(self):

        """
        Test the creation of the locus
        :return:
        """

        locus = Mikado.loci.Locus(self.t1, logger=self.logger)
        locus.json_conf = self.json_conf
        self.assertEqual(len(locus.transcripts), 1)

    def test_exclude_contained(self):

        """Test that we exclude a transcript with a contained class code (c)"""

        locus = Mikado.loci.Locus(self.t1, logger=self.logger)
        locus.json_conf = self.json_conf
        self.assertEqual(len(locus.transcripts), 1)
        locus.add_transcript_to_locus(self.t1_contained)
        self.assertEqual(len(locus.transcripts), 1)

    def test_add_contained(self):

        """Test that we add a transcript with a contained class code (c) if
        we explicitly ask for it"""

        locus = Mikado.loci.Locus(self.t1, logger=self.logger)
        locus.json_conf = self.json_conf
        locus.json_conf["pick"]["alternative_splicing"]["valid_ccodes"].append("c")
        self.assertEqual(len(locus.transcripts), 1)
        locus.add_transcript_to_locus(self.t1_contained)
        self.assertEqual(len(locus.transcripts), 2)

    def test_addValid(self):

        """Test that we can successfully add a transcript to the locus if
        it passes the muster."""

        locus = Mikado.loci.Locus(self.t1, logger=self.logger)
        locus.json_conf = self.json_conf
        self.assertEqual(len(locus.transcripts), 1)
        locus.add_transcript_to_locus(self.t1_as)
        self.assertEqual(len(locus.transcripts), 2)

    def test_excludeValid(self):

        """Test that a usually valid AS is excluded when:
        - we ask for no more than one AS event
        - we exclude its ccode (j)
        - we ask for perfect (100%) CDS overlap
        """

        locus = Mikado.loci.Locus(self.t1, logger=self.logger)
        locus.json_conf = self.json_conf
        self.assertEqual(len(locus.transcripts), 1)

        locus.json_conf["pick"]["alternative_splicing"]["max_isoforms"] = 1
        locus.add_transcript_to_locus(self.t1_as)
        self.assertEqual(len(locus.transcripts), 1)

        locus.json_conf["pick"]["alternative_splicing"]["max_isoforms"] = 3
        locus.json_conf["pick"]["alternative_splicing"]["valid_ccodes"] = ["n", "O", "h"]
        locus.add_transcript_to_locus(self.t1_as)
        self.assertEqual(len(locus.transcripts), 1)

        locus.json_conf["pick"]["alternative_splicing"]["valid_ccodes"].append("j")
        locus.json_conf["pick"]["alternative_splicing"]["min_cds_overlap"] = 100
        locus.add_transcript_to_locus(self.t1_as)
        self.assertEqual(len(locus.transcripts), 1)

    def test_exclude_with_retained_intron(self):

        """Test that a transcript with a retained intron is chucked out"""

        locus = Mikado.loci.Locus(self.t1, logger=self.logger)
        locus.json_conf = self.json_conf
        self.assertEqual(len(locus.transcripts), 1)
        locus.add_transcript_to_locus(self.t1_retained)
        self.assertEqual(len(locus.transcripts), 1)

    def test_add_with_retained_intron(self):
        """Test that a transcript with a retained intron is kept as valid
        if we change the switch."""

        locus = Mikado.loci.Locus(self.t1, logger=self.logger)
        json_conf = self.json_conf.copy()
        json_conf["pick"]["alternative_splicing"]["keep_retained_introns"] = True
        locus.json_conf = json_conf
        self.assertEqual(len(locus.transcripts), 1)
        locus.add_transcript_to_locus(self.t1_retained)
        self.assertEqual(len(locus.transcripts), 2)

    def test_exclude_opposite_strand(self):

        candidate = self.t1_as
        candidate.reverse_strand()
        logger = self.logger
        # logger.setLevel(logging.DEBUG)
        locus = Mikado.loci.Locus(self.t1, logger=logger)
        locus.json_conf = self.json_conf
        self.assertEqual(len(locus.transcripts), 1)
        locus.add_transcript_to_locus(candidate)
        self.assertEqual(len(locus.transcripts), 1)


if __name__ == '__main__':
    unittest.main(verbosity=2)
