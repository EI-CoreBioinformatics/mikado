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
from Mikado.loci import Transcript, Superlocus, Abstractlocus, Locus, Monosublocus, MonosublocusHolder, Sublocus
from Mikado.utilities.log_utils import create_null_logger, create_default_logger
from Mikado.utilities import overlap
from Mikado.utilities.intervaltree import Interval
import Mikado.loci
import pickle
import inspect
from Mikado.parsers.bed12 import BED12
from Mikado.scales.contrast import compare as c_compare


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

    logger = create_null_logger(inspect.getframeinfo(inspect.currentframe())[2])

    def setUp(self):

        gff_transcript1 = """Chr1\tfoo\ttranscript\t101\t400\t.\t+\t.\tID=t0
Chr1\tfoo\texon\t101\t400\t.\t+\t.\tID=t0:exon1;Parent=t0
Chr1\tfoo\tCDS\t101\t350\t.\t+\t.\tID=t0:exon1;Parent=t0""".split("\n")
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
        self.my_json = os.path.join(os.path.dirname(__file__), "configuration.yaml")
        self.my_json = configurator.to_json(self.my_json)
        self.assertIn("scoring", self.my_json, self.my_json.keys())

    def test_locus(self):
        """Basic testing of the Locus functionality."""

        logger = create_default_logger(inspect.getframeinfo(inspect.currentframe())[2])
        logger.setLevel("CRITICAL")
        logger.info("Started")
        self.transcript1.logger = logger
        self.transcript2.logger = logger
        self.assertTrue(self.transcript1.monoexonic)
        slocus = Superlocus(self.transcript1,
                            json_conf=self.my_json, logger=logger)
        slocus.add_transcript_to_locus(self.transcript2)
        self.assertEqual(len(slocus.transcripts), 2)
        self.assertEqual(slocus.strand, self.transcript1.strand)
        self.assertEqual(slocus.start, min(self.transcript1.start, self.transcript2.start))
        self.assertEqual(slocus.end, max(self.transcript1.end, self.transcript2.end))
        logger.info(slocus.transcripts)
        slocus.define_subloci()
        self.assertEqual(len(slocus.transcripts), 2)
        logger.info(slocus.subloci)
        logger.info(slocus.transcripts)
        self.assertEqual(len(slocus.transcripts), 2)
        self.assertEqual(len(slocus.subloci), 2)
        slocus.define_monosubloci()
        self.assertEqual(len(slocus.monosubloci), 2)
        slocus.calculate_mono_metrics()
        self.assertEqual(len(slocus.monoholders), 1)
        slocus.define_loci()
        self.assertEqual(len(slocus.loci), 1)
        self.assertEqual(list(slocus.loci[
                                  list(slocus.loci.keys())[0]].transcripts.keys())[0], "t0")
        gff_transcript3 = """Chr1\tfoo\ttranscript\t101\t1000\t.\t-\t.\tID=tminus0
Chr1\tfoo\texon\t101\t600\t.\t-\t.\tID=tminus0:exon1;Parent=tminus0
Chr1\tfoo\tCDS\t201\t500\t.\t-\t.\tID=tminus0:exon1;Parent=tminus0
Chr1\tfoo\texon\t801\t1000\t.\t-\t.\tID=tminus0:exon1;Parent=tminus0""".split("\n")
        gff_transcript3 = [GFF.GffLine(x) for x in gff_transcript3]
        transcript3 = Transcript(gff_transcript3[0])
        for exon in gff_transcript3[1:]:
                transcript3.add_exon(exon)
        transcript3.finalize()
        self.assertGreater(transcript3.combined_cds_length, 0)
        self.my_json["pick"]["clustering"]["purge"] = True
        logger.setLevel("WARNING")
        minusuperlocus = Superlocus(transcript3, json_conf=self.my_json)
        minusuperlocus.logger = logger
        minusuperlocus.define_subloci()
        self.assertGreater(len(minusuperlocus.subloci), 0)
        minusuperlocus.calculate_mono_metrics()
        self.assertGreater(len(minusuperlocus.monoholders), 0)
        minusuperlocus.define_loci()
        self.assertEqual(len(minusuperlocus.loci), 1)
        self.assertTrue(transcript3.strand != self.transcript1.strand)

    def test_verified_introns(self):

        """This method will check that during run-time, the verified introns are considered at
        the Superlocus level, not at the Sublocus one."""

        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = "1", "+", "t1"
        t1.start, t1.end = 100, 1000
        t1.add_exons([(100, 200), (300, 500), (600, 1000)])
        t1.finalize()
        t1.verified_introns.add((201, 299))
        t1.verified_introns.add((501, 599))

        t2 = Transcript()
        t2.chrom, t2.strand, t2.id = "1", "+", "t2"
        t2.start, t2.end = 100, 1000
        t2.add_exons([(100, 200), (300, 1000)])
        t2.finalize()
        t2.verified_introns.add((201, 299))

        t3 = Transcript()
        t3.chrom, t3.strand, t3.id = "1", "+", "t3"
        t3.start, t3.end = 100, 1000
        t3.add_exons([(100, 250), (300, 505), (600, 1000)])
        t3.finalize()

        jconf = configurator.to_json(None)

        loc = Superlocus(t1, json_conf=jconf)
        loc.add_transcript_to_locus(t2)
        loc.add_transcript_to_locus(t3)

        loc.define_subloci()

        self.assertEqual(t1.proportion_verified_introns, 1)
        self.assertEqual(t1.proportion_verified_introns_inlocus, 1)

        self.assertEqual(t2.proportion_verified_introns, 1)
        self.assertEqual(t2.proportion_verified_introns_inlocus, 0.5)

        self.assertEqual(t3.proportion_verified_introns, 0)
        self.assertEqual(t3.proportion_verified_introns_inlocus, 0)

    def test_boolean_requirement(self):

        logger = create_null_logger(inspect.getframeinfo(inspect.currentframe())[2])
        logger.setLevel("DEBUG")
        logger.info("Started")

        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = "1", "+", "t1"
        t1.start, t1.end = 100, 1000
        t1.add_exons([(100, 200), (300, 500), (600, 1000)])
        t1.finalize()
        t1.verified_introns.add((201, 299))
        t1.verified_introns.add((501, 599))

        t2 = Transcript()
        t2.chrom, t2.strand, t2.id = "1", "+", "t2"
        t2.start, t2.end = 100, 1000
        t2.add_exons([(100, 200), (300, 1000)])
        t2.finalize()
        t2.verified_introns.add((201, 299))

        t3 = Transcript()
        t3.chrom, t3.strand, t3.id = "1", "+", "t3"
        t3.start, t3.end = 100, 1000
        t3.add_exons([(100, 250), (300, 505), (600, 1000)])
        t3.finalize()

        jconf = configurator.to_json(None)
        # print(jconf["requirements"])

        del jconf["requirements"]

        # del jconf["scoring_file"]

        jconf["requirements"] = dict()
        jconf["requirements"]["parameters"] = dict()
        jconf["requirements"]["expression"] = ["suspicious_splicing"]
        jconf["requirements"]["parameters"]["suspicious_splicing"] = dict()
        jconf["requirements"]["parameters"]["suspicious_splicing"]["operator"] = "ne"
        jconf["requirements"]["parameters"]["suspicious_splicing"]["name"] = "suspicious_splicing"
        jconf["requirements"]["parameters"]["suspicious_splicing"]["value"] = True

        jconf["pick"]["alternative_splicing"]["report"] = False
        # Necessary to make sure that the externally-specified requirements are taken in
        configurator.check_all_requirements(jconf)
        self.assertEqual(
            jconf["requirements"]["expression"],
            "evaluated[\"suspicious_splicing\"]")

        jconf = configurator.check_json(jconf)

        self.assertEqual(
            jconf["requirements"]["expression"],
            "evaluated[\"suspicious_splicing\"]")

        logger = create_default_logger(inspect.getframeinfo(inspect.currentframe())[2])
        for suspicious in (False, True):
            with self.subTest(suspicious=suspicious):
                loc = Superlocus(t1, json_conf=jconf, logger=logger)
                t2.attributes["canonical_on_reverse_strand"] = suspicious
                loc.add_transcript_to_locus(t2)
                loc.add_transcript_to_locus(t3)
                self.assertEqual(len(loc.transcripts), 3)
                loc.define_subloci()
                self.assertEqual(len(loc.transcripts), 3 if not suspicious else 2)


class ASeventsTester(unittest.TestCase):

    logger = create_null_logger("ASevents")

    def setUp(self):
        
        self.conf = configurator.to_json(None)
        self.conf["pick"]["alternative_splicing"] = dict()
        self.conf["pick"]["alternative_splicing"]["report"] = True
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

        # self.t1.add_exons([(101, 500), (601, 700), (1001, 1300), (1401, 1500)],
        #                   "exon")
        # self.t1.add_exons([(401, 500), (601, 700), (1001, 1300), (1401, 1440)],
        #                   "CDS")

        t2.add_exons([(101, 500), (601, 700), (1001, 1300), (1401, 1460), (1501, 1600)],
                     "exon")
        t2.add_exons([(401, 500), (601, 700), (1001, 1300), (1401, 1440)],
                     "CDS")
        t2.finalize()

        # self.locus.add_transcript_to_locus(t2)
        self.assertEqual(self.locus.is_alternative_splicing(t2)[:2], (True, "J"))
        self.locus.json_conf["pick"]["clustering"]["cds_only"] = True

        self.assertEqual(self.locus.is_alternative_splicing(t2)[:2], (False, "="))

    def test_redundant_cds_non_redundant_cdna(self):

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
        self.assertEqual(len(self.locus.transcripts), 2, self.locus.transcripts)


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

    def test_intron_contained_in_exon(self):

        """Here the intron is completely contained within an exon. Returns true."""

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
        self.assertTrue(MonosublocusHolder.is_intersecting(self.t1, t2))

    def test_intron_not_contained_in_exon(self):

        self.t1.strip_cds()
        t2 = Transcript()
        t2.chrom = "Chr1"
        t2.strand = "+"
        t2.score = 1
        t2.id = "G2.1"
        t2.parent = "G2"
        t2.start = 1400
        t2.end = 3000
        t2.add_exons([(1400, 1560), (2800, 3000)])
        t2.finalize()

        logger = create_default_logger("test_intron_not_contained_in_exon")

        for min_cdna_overlap in (0.01, 1):
            with self.subTest(min_cdna_overlap=min_cdna_overlap):
                self.assertIs(MonosublocusHolder.is_intersecting(
                    self.t1, t2,
                    logger=logger,
                    cds_only=False,
                    min_cdna_overlap=min_cdna_overlap,
                    min_cds_overlap=min_cdna_overlap), (min_cdna_overlap < 0.28))

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
        logger = create_default_logger(inspect.getframeinfo(inspect.currentframe())[2])

        for min_cds_overlap in [0.05, 0.1, 0.15, 0.2, 0.5]:
            with self.subTest(min_cds_overlap=min_cds_overlap):
                self.assertIs(MonosublocusHolder.is_intersecting(self.t1, t2,
                                                                 cds_only=True,
                                                                 logger=logger,
                                                                 min_cds_overlap=min_cds_overlap,
                                                                 min_cdna_overlap=0.01),
                              (min_cds_overlap <= 0.19))

        t2.strip_cds()
        t2.finalized = False
        t2.add_exons([(1461, 1560), (1801, 1850)], "CDS")
        t2.finalize()
        self.assertGreater(len(t2.introns), 0)
        self.assertGreater(len(t2.combined_cds_introns), 0)
        # No CDS overlap this time, but cDNA overlap.
        for cds_only in (True, False):
            with self.subTest(cds_only=cds_only):
                self.assertIs(MonosublocusHolder.is_intersecting(self.t1,
                                                            t2,
                                                            cds_only=cds_only,
                                                            logger=logger), not cds_only)

        t2 = Transcript()
        t2.chrom = "Chr1"
        t2.strand = "+"
        t2.score = 1
        t2.id = "G2.1"
        t2.parent = "G2"
        t2.start = 1350
        t2.end = 3850
        t2.add_exons([(1350, 1560), (2801, 3850)])
        t2.add_exons([(1401, 1560), (2801, 3850)], "CDS")
        # logger.setLevel("DEBUG")
        t2.logger = logger
        t2.finalize()
        self.assertTrue(t2.is_coding)
        for min_overlap in [0.1, 0.2, 0.3, 0.5]:
            with self.subTest(min_overlap=min_overlap):
                self.assertIs(MonosublocusHolder.is_intersecting(self.t1, t2,
                                                                 cds_only=False,
                                                                 min_cds_overlap=0.07,
                                                                 min_cdna_overlap=min_overlap,
                                                                 logger=logger), (min_overlap <= 0.12))

        self.assertTrue(t2.is_coding)

        for min_overlap in [0.01, 0.05, 0.1, 0.2]:
            with self.subTest(min_overlap=min_overlap):
                self.assertIs(MonosublocusHolder.is_intersecting(self.t1,
                                                                 t2,
                                                                 cds_only=True,
                                                                 min_cds_overlap=min_overlap,
                                                                 min_cdna_overlap=min_overlap,
                                                                 logger=logger), (min_overlap <= 0.07))

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

    def test_sameness(self):

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

    def test_holder_clustering(self):

        """This test has been added starting from the annotation of IWGSC.
        It verifies that in a complex locus we create the holders correctly."""

        chrom, strand = "chr7A", "+"
        transcripts = dict()
        transcripts["TA_PGSB_v1_dez2016_mRNA_662095"] = [
            [(711041145, 711041431), (711041641, 711042803), (711059154, 711059935)],
            36.17, "sublocus:chr3A+:711041145-711059935.multi"]
        transcripts["TA_PGSB_v1_dez2016_mRNA_662100"] = [
            [(711056723, 711056806), (711056870, 711057549), (711057994, 711059935)],
            49.8, "sublocus:chr3A+:711040605-711059935.multi"]
        transcripts["TA_PGSB_v1_dez2016_mRNA_662101"] = [
            [(711056723, 711057549), (711057991, 711059935)],
            48.02, "sublocus:chr3A+:711056723-711059935.multi"]
        transcripts["TA_PGSB_v1_dez2016_mRNA_662106"] = [
            [(711056723, 711057549), (711057995, 711058007)],
            39.51, "sublocus:chr3A+:711056723-711058007.multi"]
        transcripts["TA_PGSB_v1_dez2016_mRNA_662109"] = [
            [(711056723, 711057141), (711057213, 711057237)],
            35.85, "sublocus:chr3A+:711056723-711057237.multi"]
        transcripts["TA_PGSB_v1_dez2016_mRNA_662111"] = [
            [(711056723, 711057549), (711057994, 711059935)],
            49.97, "sublocus:chr3A+:711040605-711059935.multi"]
        transcripts["TA_PGSB_v1_dez2016_mRNA_662116"] = [
            [(711056723, 711057610)],
            30.7, "sublocus:chr3A+:711056723-711057610.mono"
        ]
        transcripts["TA_PGSB_v1_dez2016_mRNA_662121"] = [
            [(711058325, 711059913), (711060068, 711060089)],
            36.48, "sublocus:chr3A+:711058325-711060089.multi"
        ]

        superlocus = None
        json_conf = configurator.to_json(None)

        for transcript in transcripts:
            tr = Transcript()
            exons, score, sublocus = transcripts[transcript]
            tr.chrom, tr.strand, tr.id = chrom, strand, transcript
            tr.add_exons(exons, features="exon")
            tr.add_exons(exons, features="CDS")
            tr.finalize()
            subl = Sublocus(tr, json_conf=json_conf)
            if superlocus is None:
                superlocus = Superlocus(tr, json_conf=json_conf)
            else:
                superlocus.add_transcript_to_locus(tr, check_in_locus=False)
            superlocus.subloci.append(subl)
            superlocus.scores[transcript] = score

        superlocus.subloci_defined = True
        self.assertEqual(len(superlocus.subloci), len(transcripts))
        superlocus.calculate_mono_metrics()
        self.assertEqual(len(superlocus.monoholders), 1,
                         "\n".join([", ".join(list(_.transcripts.keys())) for _ in superlocus.monoholders]))

    def test_alternative_splicing_monoexonic_not_enough_overlap(self):

        """This test verifies that while we can cluster together the transcripts at the holder stage,
        if the overlap is not enough they will fail to be recognised as valid AS events."""

        jconf = configurator.to_json(None)

        t1, t2 = Transcript(), Transcript()
        t1.chrom, t2.chrom = "1", "1"
        t1.strand, t2.strand = "+", "+"
        t1.add_exons([(1260208, 1260482), (1262216, 1262412), (1262621, 1263851)])
        t1.add_exons([(1262291, 1262412), (1262621, 1263143)], features="CDS")
        t1.id = "cls-0-sta-combined-0_1.27.12"

        t2.add_exons([(1262486, 1264276)])
        t2.add_exons([(1263571, 1264236)], features="CDS")
        t2.id = "trn-0-sta-combined-0_1_TRINITY_GG_1373_c0_g2_i1.path1"

        self.assertTrue(MonosublocusHolder.is_intersecting(
            t1, t2, cds_only=False,
            min_cdna_overlap=jconf["pick"]["alternative_splicing"]["min_cdna_overlap"],
            min_cds_overlap=jconf["pick"]["alternative_splicing"]["min_cds_overlap"],
            simple_overlap_for_monoexonic=True))
        self.assertFalse(MonosublocusHolder.is_intersecting(
            t1, t2, cds_only=False,
            min_cdna_overlap=jconf["pick"]["clustering"]["min_cdna_overlap"],
            min_cds_overlap=jconf["pick"]["clustering"]["min_cds_overlap"],
            simple_overlap_for_monoexonic=False))

        for simple in (True, False):
            with self.subTest(simple=simple):
                jconf["pick"]["clustering"]["simple_overlap_for_monoexonic"] = simple

                slocus = Superlocus(t1, json_conf=jconf)
                slocus.add_transcript_to_locus(t2)
                locus = Locus(t1, json_conf=jconf)
                slocus.loci[locus.id] = locus
                slocus.define_alternative_splicing()
                self.assertEqual(len(slocus.loci[locus.id].transcripts), 1)


class TestLocus(unittest.TestCase):

    """
    This unit test is focused on the locus definition and alternative splicings.
    """

    logger = Mikado.utilities.log_utils.create_default_logger("tester")

    def setUp(self):

        """Set up for the unit test."""

        # Mock dictionary to be used for the alternative splicing checks
        self.json_conf = configurator.to_json(None)
        # self.json_conf["pick"] = dict()
        self.json_conf["pick"]["alternative_splicing"] = dict()
        self.json_conf["pick"]["alternative_splicing"]["report"] = True
        self.json_conf["pick"]["alternative_splicing"]["max_utr_length"] = 2000
        self.json_conf["pick"]["alternative_splicing"]["max_fiveutr_length"] = 1000
        self.json_conf["pick"]["alternative_splicing"]["max_threeutr_length"] = 1000
        self.json_conf["pick"]["alternative_splicing"]["max_isoforms"] = 3
        self.json_conf["pick"]["alternative_splicing"]["keep_retained_introns"] = False
        self.json_conf["pick"]["alternative_splicing"]["min_cds_overlap"] = 0
        self.json_conf["pick"]["alternative_splicing"]["min_cdna_overlap"] = 0
        self.json_conf["pick"]["alternative_splicing"]["min_score_perc"] = 0.1
        self.json_conf["pick"]["alternative_splicing"]["valid_ccodes"] = ["j", "G", "g"]
        self.json_conf["pick"]["alternative_splicing"]["redundant_ccodes"] = ["c", "=", "_", "m", "n"]
        self.json_conf["pick"]["alternative_splicing"]["only_confirmed_introns"] = False

        self.json_conf = configurator.check_json(self.json_conf)

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

        locus.json_conf["pick"]["alternative_splicing"]["max_isoforms"] = 3
        locus.json_conf["pick"]["alternative_splicing"]["valid_ccodes"] = ["n", "O", "h"]
        locus.add_transcript_to_locus(self.t1_as)
        self.assertEqual(len(locus.transcripts), 1)

        locus.json_conf["pick"]["alternative_splicing"]["valid_ccodes"].append("j")
        locus.json_conf["pick"]["alternative_splicing"]["min_cds_overlap"] = 1

        locus.add_transcript_to_locus(self.t1_as)
        self.assertEqual(len(locus.transcripts), 1)

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

    def test_serialisation(self):

        """Check that the main types can be serialised correctly."""

        candidate = self.t1
        pickle.dumps(candidate)

        json_conf = configurator.to_json(None)

        for obj in Superlocus, Sublocus, Locus:
            with self.subTest(obj=obj):
                locus = obj(candidate, json_conf=json_conf)
                pickle.dumps(locus)

    def test_double_orf(self):

        t = Transcript()
        t.add_exons([(101, 1000), (1101, 1200), (2001, 2900)])
        t.id = "t1"
        t.strand = "+"

        orf1 = BED12()
        orf1.transcriptomic = True
        orf1.chrom = t.id
        orf1.start = 1
        orf1.end = sum([_[1] - _[0] + 1 for _ in t.exons])
        orf1.strand = "+"
        orf1.name = "t1.orf1"
        orf1.block_sizes = (900,)
        orf1.thick_start = 1
        orf1.thick_end = 900
        orf1.block_starts = (1,)
        orf1.block_count = 1

        orf2 = BED12()
        orf2.transcriptomic = True
        orf2.strand = "+"
        orf2.chrom = t.id
        orf2.start = 1
        orf2.end = sum([_[1] - _[0] + 1 for _ in t.exons])
        orf2.name = "t1.orf2"
        orf2.block_sizes = (900,)
        orf2.thick_start = 1001
        orf2.thick_end = 1900
        orf2.block_starts = (1,)
        orf2.block_count = 1

        self.assertFalse(orf1.invalid)
        self.assertFalse(orf2.invalid)

        t.load_orfs([orf1, orf2])
        self.assertEqual(t.number_internal_orfs, 2)

        locus = Locus(t)
        locus.calculate_scores()
        self.assertTrue(list(locus.scores.keys()), [t.id])
        rows = list(locus.print_scores())
        self.assertEqual(len(rows), 1, rows)
        self.assertEqual(rows[0]["tid"], t.id, rows[0])

    def test_remove_AS_overlapping(self):

        logger = create_null_logger(inspect.getframeinfo(inspect.currentframe())[2],
                                    level="WARNING")

        t1, t2, t1_1, t2_1 = Transcript(), Transcript(), Transcript(), Transcript()
        t1.chrom = t2.chrom = t1_1.chrom = t2_1.chrom = "1"
        t1.id, t2.id, t1_1.id, t2_1.id = "t1", "t2", "t1_1", "t2_1"
        t1.strand = t2.strand = t1_1.strand = t2_1.strand = "+"
        t1.add_exons([(101, 500), (801, 1000)])
        t1.add_exons([(101, 500), (801, 1000)], features="CDS")
        t1_1.add_exons([(101, 500), (901, 1100), (1301, 1550)])
        t1_1.add_exons([(101, 500), (901, 1100), (1301, 1550)], features="CDS")
        t2.add_exons([(1601, 1800), (1901, 2000)])
        t2.add_exons([(1601, 1800), (1901, 2000)], features="CDS")
        t2_1.add_exons([(1351, 1550), (1651, 1850), (1901, 2000)])
        t2_1.add_exons([(1351, 1550), (1651, 1850), (1901, 2000)], features="CDS")

        for tr in [t1, t2, t1_1, t2_1]:
            with self.subTest(tr=tr):
                tr.finalize()
                self.assertGreater(tr.combined_cds_length, 0, tr.id)

        conf = configurator.to_json(None)
        conf["pick"]["alternative_splicing"]["valid_ccodes"] = ["j", "J", "g", "G"]
        conf["pick"]["alternative_splicing"]["only_confirmed_introns"] = False

        conf["as_requirements"] = {"_expression": "cdna_length",
                                   "expression": "evaluated['cdna_length']",
                                   "parameters": {
                                       "cdna_length": {"operator": "gt", "value": 0, "name": "cdna_length"}
                                   }}

        with self.subTest():
            superlocus_one = Superlocus(t1, json_conf=conf)
            superlocus_one.add_transcript_to_locus(t1_1)
            locus_one = Locus(t1, json_conf=conf)
            locus_one.logger = logger
            superlocus_one.loci[locus_one.id] = locus_one
            superlocus_one.loci_defined = True
            superlocus_one.logger = logger
            superlocus_one.define_alternative_splicing()
            self.assertEqual(len(superlocus_one.loci), 1)
            locus_id = [_ for _ in superlocus_one.loci.keys() if
                        t1.id in superlocus_one.loci[_].transcripts][0]
            self.assertEqual(len(superlocus_one.loci[locus_id].transcripts), 2)
        with self.subTest():
            superlocus_two = Superlocus(t2, json_conf=conf)
            superlocus_two.add_transcript_to_locus(t2_1)
            locus_two = Locus(t2, json_conf=conf)
            superlocus_two.loci[locus_two.id] = locus_two
            superlocus_two.loci_defined = True
            superlocus_two.logger = logger
            superlocus_two.define_alternative_splicing()
            self.assertEqual(len(superlocus_two.loci), 1)
            locus_id = [_ for _ in superlocus_two.loci.keys() if
                        t2.id in superlocus_two.loci[_].transcripts][0]
            self.assertEqual(len(superlocus_two.loci[locus_id].transcripts), 2)
        with self.subTest():
            superlocus = Superlocus(t1, json_conf=conf)
            superlocus.add_transcript_to_locus(t2, check_in_locus=False)
            superlocus.add_transcript_to_locus(t1_1)
            superlocus.add_transcript_to_locus(t2_1)
            locus_one = Locus(t1, json_conf=conf)
            locus_two = Locus(t2, json_conf=conf)
            superlocus.loci[locus_one.id] = locus_one
            superlocus.loci[locus_two.id] = locus_two
            superlocus.loci_defined = True
            superlocus.logger = logger
            self.assertEqual(len(superlocus.loci), 2)
            superlocus.define_alternative_splicing()
            locus_one_id = [_ for _ in superlocus.loci.keys() if
                            t1.id in superlocus.loci[_].transcripts][0]
            locus_two_id = [_ for _ in superlocus.loci.keys() if
                            t2.id in superlocus.loci[_].transcripts][0]
            self.assertEqual(len(superlocus.loci), 2)
            self.assertEqual(len(superlocus.loci[locus_one_id].transcripts), 1)
            self.assertEqual(len(superlocus.loci[locus_two_id].transcripts), 1)


class RetainedIntronTester(unittest.TestCase):

    def setUp(self):
        self.my_json = os.path.join(os.path.dirname(__file__), "configuration.yaml")
        self.my_json = configurator.to_json(self.my_json)

    def test_real_retained_pos(self):

        """Here we verify that a real retained intron is called as such"""

        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = 1, "+", "t1"
        t1.add_exons([(101, 500), (801, 1000), (1201, 1300), (1501, 1800)])
        t1.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1300),  #100
                      (1501, 1530)  # 30
                      ], features="CDS")
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.strand, t2.id = 1, "+", "t2"
        t2.add_exons([(101, 500), (801, 1000), (1201, 1600)])
        t2.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1420),  # 220
                      ], features="CDS")
        t2.finalize()

        t3 = Transcript()
        t3.chrom, t3.strand, t3.id = 1, "+", "t3"
        t3.add_exons([(101, 500), (801, 970), (1100, 1180)])
        t3.add_exons([(101, 500), (801, 970), (1100, 1130)], features="CDS")
        t3.finalize()

        for pred, retained in [(t2, True), (t3, False)]:
            with self.subTest(pred=pred, retained=retained):
                sup = Superlocus(t1, json_conf=self.my_json)
                sup.add_transcript_to_locus(pred)
                sup.json_conf["pick"]["run_options"]["consider_truncated_for_retained"] = True
                sup.find_retained_introns(pred)
                self.assertEqual((len(sup.transcripts[pred.id].retained_introns) > 0),
                                 retained)

    def test_retained_pos_truncated(self):
        """Here we verify that a real retained intron is called as such,
        even when the transcript is truncated."""

        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = 1, "+", "t1"
        t1.add_exons([(101, 500), (801, 1000), (1201, 1300), (1501, 1800)])
        t1.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1300),  #100
                      (1501, 1530)  # 30
                      ], features="CDS")
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.strand, t2.id = 1, "+", "t2"
        t2.add_exons([(101, 500), (801, 1000), (1201, 1420)])
        t2.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1420),  # 220
                      ], features="CDS")
        t2.finalize()
        self.assertEqual(t2.combined_cds_end, 1420)

        t3 = Transcript()
        t3.chrom, t3.strand, t3.id = 1, "+", "t3"
        t3.add_exons([(101, 500), (801, 970), (1100, 1130)])
        t3.add_exons([(101, 500), (801, 970), (1100, 1130)], features="CDS")
        t3.finalize()

        for pred, retained in [(t2, True), (t3, False)]:
            with self.subTest(pred=pred, retained=retained):
                sup = Superlocus(t1, json_conf=self.my_json)
                sup.add_transcript_to_locus(pred)
                sup.json_conf["pick"]["run_options"]["consider_truncated_for_retained"] = True
                sup.find_retained_introns(pred)
                self.assertEqual((len(sup.transcripts[pred.id].retained_introns) > 0),
                                 retained)
                # Now check that things function also after unpickling
                unpickled_t1 = pickle.loads(pickle.dumps(t1))
                unpickled_other = pickle.loads(pickle.dumps(pred))
                sup = Superlocus(unpickled_t1, json_conf=self.my_json)
                sup.add_transcript_to_locus(unpickled_other)
                sup.json_conf["pick"]["run_options"]["consider_truncated_for_retained"] = True
                sup.find_retained_introns(pred)
                self.assertEqual((len(sup.transcripts[pred.id].retained_introns) > 0),
                                 retained)

    def test_real_retained_pos_truncated_skip(self):
        """Here we verify that a real retained intron is *NOT* called as such when
        the transcript is truncated and we elect not to investigate the 3' end."""

        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = 1, "+", "t1"
        t1.add_exons([(101, 500), (801, 1000), (1201, 1300), (1501, 1800)])
        t1.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1300),  #100
                      (1501, 1530)  # 30
                      ], features="CDS")
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.strand, t2.id = 1, "+", "t2"
        t2.add_exons([(101, 500), (801, 1000), (1201, 1420)])
        t2.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1420),  # 220
                      ], features="CDS")
        t2.finalize()
        self.assertEqual(t2.combined_cds_end, 1420)

        sup = Superlocus(t1, json_conf=self.my_json)
        sup.add_transcript_to_locus(t2)
        sup.json_conf["pick"]["run_options"]["consider_truncated_for_retained"] = False

        sup.find_retained_introns(t2)

        self.assertEqual(sup.transcripts["t2"].retained_introns, (),)

    def test_real_retained_neg_truncated(self):
        """Here we verify that a real retained intron is called as such,
        even when the transcript is truncated."""

        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = 1, "-", "t1"
        t1.add_exons([(101, 500), (801, 1000), (1201, 1300), (1501, 1800)])
        t1.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1300),  #100
                      (1501, 1530)  # 30
                      ], features="CDS")
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.strand, t2.id = 1, "-", "t2"
        t2.add_exons([(601, 1000), (1201, 1300), (1501, 1800)])
        t2.add_exons([(601, 1000),  # 200
                      (1201, 1300),  #100
                      (1501, 1530)  # 30
                      ], features="CDS")
        t2.finalize()
        self.assertEqual(t2.combined_cds_end, 601)

        t3 = Transcript()
        t3.chrom, t3.strand, t3.id = 1, "-", "t3"
        t3.add_exons([(551, 580), (801, 1000), (1201, 1300), (1501, 1800)])
        t3.add_exons([(551, 580),
                      (801, 1000),  # 200
                      (1201, 1300),  #100
                      (1501, 1530)  # 30
                      ], features="CDS")
        t3.finalize()
        self.assertEqual(t3.combined_cds_end, 551)

        for pred, retained in [(t2, True), (t3, False)]:
            with self.subTest(pred=pred, retained=retained):
                sup = Superlocus(t1, json_conf=self.my_json)
                sup.add_transcript_to_locus(pred)
                sup.json_conf["pick"]["run_options"]["consider_truncated_for_retained"] = True
                sup.find_retained_introns(pred)
                self.assertEqual((len(sup.transcripts[pred.id].retained_introns) > 0),
                                 retained)

    def test_real_retained_neg_truncated_skip(self):
        """Here we verify that a real retained intron is *NOT* called as such when
        the transcript is truncated and we elect not to investigate the 3' end."""

        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = 1, "-", "t1"
        t1.add_exons([(101, 500), (801, 1000), (1201, 1300), (1501, 1800)])
        t1.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1300),  #100
                      (1501, 1530)  # 30
                      ], features="CDS")
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.strand, t2.id = 1, "-", "t2"
        t2.add_exons([(601, 1000), (1201, 1300), (1501, 1800)])
        t2.add_exons([(601, 1000),  # 200
                      (1201, 1300),  #100
                      (1501, 1530)  # 30
                      ], features="CDS")
        t2.finalize()
        self.assertEqual(t2.combined_cds_end, 601)

        sup = Superlocus(t1, json_conf=self.my_json)
        sup.add_transcript_to_locus(t2)
        sup.json_conf["pick"]["run_options"]["consider_truncated_for_retained"] = False

        sup.find_retained_introns(t2)

        self.assertEqual(sup.transcripts["t2"].retained_introns, ())

    def test_real_retained_pos_noCDS(self):
        """Here we verify that a real retained intron is called as such, even when the transcript lacks a CDS"""

        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = 1, "+", "t1"
        t1.add_exons([(101, 500), (801, 1000), (1201, 1300), (1501, 1800)])
        t1.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1300),  # 100
                      (1501, 1530)  # 30
                      ], features="CDS")
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.strand, t2.id = 1, "+", "t2"
        t2.add_exons([(101, 500), (801, 1000), (1201, 1600)])
        t2.finalize()

        sup = Superlocus(t1, json_conf=self.my_json)
        sup.add_transcript_to_locus(t2)

        sup.find_retained_introns(t2)

        self.assertEqual(sup.transcripts["t2"].retained_introns, ((1201, 1600),))

    def test_not_retained_pos(self):

        """Here we verify that a false retained intron is not called as such"""

        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = 1, "+", "t1"
        t1.add_exons([(101, 500), (801, 1000), (1201, 1300), (1501, 1800)])
        t1.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1300),  # 100
                      (1501, 1530)  # 30
                      ], features="CDS")
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.strand, t2.id = 1, "+", "t2"
        t2.add_exons([(101, 500), (801, 1000), (1201, 1600)])
        t2.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1540),  # 340
                      ], features="CDS")
        t2.finalize()

        t3 = Transcript()
        t3.chrom, t3.strand, t3.id = 1, "+", "t3"
        t3.add_exons([(101, 500), (801, 970), (1100, 1130)])
        t3.add_exons([(101, 500), (801, 970), (1100, 1130)], features="CDS")
        t3.finalize()

        for pred in [t2, t3]:
            with self.subTest(pred=pred):
                sup = Superlocus(t1, json_conf=self.my_json)
                sup.add_transcript_to_locus(pred)
                sup.find_retained_introns(pred)
                self.assertEqual(sup.transcripts[pred.id].retained_intron_num, 0)
                unpickled_t1 = pickle.loads(pickle.dumps(t1))
                unpickled_other = pickle.loads(pickle.dumps(pred))
                sup = Superlocus(unpickled_t1, json_conf=self.my_json)
                sup.add_transcript_to_locus(unpickled_other)
                sup.find_retained_introns(unpickled_other)
                self.assertEqual(sup.transcripts[unpickled_other.id].retained_intron_num, 0)

    def test_real_retained_neg(self):
        """Here we verify that a real retained intron is called as such"""

        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = 1, "-", "t1"
        t1.add_exons([(101, 500), (801, 1000), (1201, 1300), (1501, 1800)])
        t1.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1300),  # 100
                      (1501, 1530)  # 30
                      ], features="CDS")
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.strand, t2.id = 1, "-", "t2"
        t2.add_exons([(401, 1000), (1201, 1300), (1501, 1800)])
        t2.add_exons([(1501, 1530),  # 30
                      (1201, 1300),  # 100
                      (771, 1000)  # 230
                      ], features="CDS")
        t2.finalize()

        with self.subTest():
            sup = Superlocus(t1, json_conf=self.my_json)
            sup.add_transcript_to_locus(t2)

            sup.find_retained_introns(t2)
            self.assertEqual(sup.transcripts["t2"].retained_introns, ((401, 1000),))

        with self.subTest():
            unpickled_t1 = pickle.loads(pickle.dumps(t1))
            unpickled_other = pickle.loads(pickle.dumps(t2))
            sup = Superlocus(unpickled_t1, json_conf=self.my_json)
            sup.add_transcript_to_locus(unpickled_other)
            sup.find_retained_introns(unpickled_other)
            self.assertEqual(sup.transcripts["t2"].retained_introns, ((401, 1000),))

    def test_not_real_retained_neg(self):
        """Here we verify that a real retained intron is called as such"""

        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = 1, "-", "t1"
        t1.add_exons([(101, 500), (801, 1000), (1201, 1300), (1501, 1800)])
        t1.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1300),  # 100
                      (1501, 1530)  # 30
                      ], features="CDS")
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.strand, t2.id = 1, "-", "t2"
        t2.add_exons([(601, 1000), (1201, 1300), (1501, 1800)])
        t2.add_exons([(1501, 1530),  # 30
                      (1201, 1300),  # 100
                      (771, 1000)  # 230
                      ], features="CDS")
        t2.finalize()

        t3 = Transcript()
        t3.chrom, t3.strand, t3.id = 1, "-", "t3"
        t3.add_exons([(401, 1000), (1201, 1300), (1501, 1800)])
        t3.add_exons([(831, 1000),  # 200
                      (1201, 1300),
                      (1501, 1530)
                      ], features="CDS")
        t3.finalize()

        self.assertFalse(
            Abstractlocus._is_exon_retained_in_transcript((401, 1000), [Interval(401, 830)], t1))

        for alt in [t2, t3]:
            unpickled_t1 = pickle.loads(pickle.dumps(t1))
            unpickled_alt = pickle.loads(pickle.dumps(alt))
            with self.subTest(alt=alt):
                sup = Superlocus(t1, json_conf=self.my_json)
                sup.find_retained_introns(alt)
                self.assertEqual(alt.retained_intron_num, 0,
                                 alt.retained_introns)
            with self.subTest(alt=alt):
                sup = Superlocus(unpickled_t1, json_conf=self.my_json)
                sup.find_retained_introns(unpickled_alt)
                self.assertEqual(unpickled_alt.retained_intron_num, 0,
                                 unpickled_alt.retained_introns)

    def test_not_retained_neg(self):
        """Here we verify that a false retained intron is not called as such"""

        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = 1, "-", "t1"
        t1.add_exons([(101, 500), (801, 1000), (1201, 1300), (1501, 1800)])
        t1.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1300),  # 100
                      (1501, 1530)  # 30
                      ], features="CDS")
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.strand, t2.id = 1, "-", "t2"
        t2.add_exons([(301, 1000), (1201, 1300), (1501, 1800)])
        t2.add_exons([(1501, 1530),  # 30
                      (1201, 1300),  # 100
                      (471, 1000)  # 230
                      ], features="CDS")
        t2.finalize()

        sup = Superlocus(t1, json_conf=self.my_json)
        sup.add_transcript_to_locus(t2)

        self.assertEqual(t2.cds_tree.find(301, 1000),
                         [Interval(471, 1000)])

        self.assertEqual(Abstractlocus._exon_to_be_considered((301, 1000), t2),
                         (True, [(301, 470)]),
                         Abstractlocus._exon_to_be_considered((301, 1000), t2))

        self.assertFalse(Abstractlocus._is_exon_retained_in_transcript((301, 1000),
                                                                       [(301, 470)],
                                                                       t1))

        sup.find_retained_introns(t2)

        self.assertEqual(sup.transcripts["t2"].retained_intron_num, 0,
                         sup.transcripts["t2"].retained_introns)

    def test_exon_switching_pos(self):

        """Checking that an exon switching is treated correctly as a NON-retained intron. Positive strand case"""

        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = 1, "+", "t1"
        t1.add_exons([(101, 500), (801, 1000), (1201, 1300), (2501, 2800)])
        t1.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1300),  # 100
                      (2501, 2530)  # 30
                      ], features="CDS")
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.strand, t2.id = 1, "+", "t2"
        t2.add_exons([(101, 500), (801, 1000), (1201, 1300), (1501, 1800)])
        t2.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1300),  # 100
                      (1501, 1530)  # 30
                      ], features="CDS")
        t2.finalize()

        sup = Superlocus(t1, json_conf=self.my_json)
        sup.add_transcript_to_locus(t2)

        sup.find_retained_introns(t2)

        self.assertEqual(sup.transcripts["t2"].retained_intron_num, 0)

    def test_exon_switching_pos_noCDS(self):
        """Checking that an exon switching is treated correctly as a NON-retained intron even when the CDS is absent.
        Positive strand case"""

        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = 1, "+", "t1"
        t1.add_exons([(101, 500), (801, 1000), (1201, 1300), (2501, 2800)])
        t1.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1300),  # 100
                      (2501, 2530)  # 30
                      ], features="CDS")
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.strand, t2.id = 1, "+", "t2"
        t2.add_exons([(101, 500), (801, 1000), (1201, 1300), (1501, 1800)])
        # t2.add_exons([(201, 500),  # 300
        #               (801, 1000),  # 200
        #               (1201, 1300),  # 100
        #               (1501, 1530)  # 30
        #               ], features="CDS")
        t2.finalize()

        sup = Superlocus(t1, json_conf=self.my_json)
        sup.add_transcript_to_locus(t2)

        sup.find_retained_introns(t2)

        self.assertEqual(sup.transcripts["t2"].retained_intron_num, 0,
                         sup.transcripts["t2"].retained_introns)

    def test_exon_switching_neg(self):
        """Checking that an exon switching is treated correctly as a NON-retained intron. Positive strand case"""

        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = 1, "-", "t1"
        t1.add_exons([(101, 500), (801, 1000), (2201, 2300), (2501, 2800)])
        t1.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (2201, 2300),  # 100
                      (2501, 2530)  # 30
                      ], features="CDS")
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.strand, t2.id = 1, "-", "t2"
        t2.add_exons([(101, 500), (1701, 2000), (2201, 2300), (2501, 2800)])
        t2.add_exons([
                      (1801, 2000),  # 200
                      (2201, 2300),  # 100
                      (2501, 2530)  # 30
                      ], features="CDS")
        t2.finalize()
        self.assertEqual(len(t2.cds_tree), len(t2.combined_cds))
        self.assertEqual(len(t2.cds_tree), 3)

        sup = Superlocus(t1, json_conf=self.my_json)
        sup.add_transcript_to_locus(t2)

        sup.find_retained_introns(t2)

        self.assertEqual(sup.transcripts["t2"].retained_intron_num, 0)

    def test_exon_switching_neg_noCDS(self):
        """Checking that an exon switching is treated correctly as a NON-retained intron even when the CDS is absent.
        Positive strand case"""

        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = 1, "-", "t1"
        t1.add_exons([(101, 500), (801, 1000), (1201, 1300), (2501, 2800)])
        t1.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1300),  # 100
                      (2501, 2530)  # 30
                      ], features="CDS")
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.strand, t2.id = 1, "-", "t2"
        t2.add_exons([(101, 500), (801, 1000), (1201, 1300), (1501, 1800)])
        # t2.add_exons([(201, 500),  # 300
        #               (801, 1000),  # 200
        #               (1201, 1300),  # 100
        #               (1501, 1530)  # 30
        #               ], features="CDS")
        t2.finalize()

        sup = Superlocus(t1, json_conf=self.my_json)
        sup.add_transcript_to_locus(t2)

        sup.find_retained_introns(t2)

        self.assertEqual(sup.transcripts["t2"].retained_intron_num, 0)

    def test_mixed_strands(self):

        """Verify that no retained intron is called if the strands are mixed."""

        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = 1, "+", "t1"
        t1.add_exons([(101, 500), (801, 1000), (1201, 1300), (1501, 1800)])
        t1.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1300),  # 100
                      (1501, 1530)  # 30
                      ], features="CDS")
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.strand, t2.id = 1, "-", "t2"
        t2.add_exons([(601, 1000), (1201, 1300), (1501, 1800)])
        t2.add_exons([(1501, 1530),  # 30
                      (1201, 1300),  # 100
                      (771, 1000)  # 230
                      ], features="CDS")
        t2.finalize()

        sup = Superlocus(t1, json_conf=self.my_json, stranded=False)
        sup.add_transcript_to_locus(t2)

        sup.find_retained_introns(t2)

        self.assertEqual(sup.transcripts["t2"].retained_intron_num, 0)


class PicklingTest(unittest.TestCase):

    def setUp(self):
        t1 = Transcript()
        t1.chrom, t1.strand, t1.id = 1, "+", "t1"
        t1.add_exons([(101, 500), (801, 1000), (1201, 1300), (1501, 1800)])
        t1.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1300),  # 100
                      (1501, 1530)  # 30
                      ], features="CDS")
        t1.finalize()

        t2 = Transcript()
        t2.chrom, t2.strand, t2.id = 1, "+", "t2"
        t2.add_exons([(101, 500), (801, 1000), (1201, 1600)])
        t2.add_exons([(201, 500),  # 300
                      (801, 1000),  # 200
                      (1201, 1420),  # 220
                      ], features="CDS")
        t2.finalize()

        t3 = Transcript()
        t3.chrom, t3.strand, t3.id = 1, "+", "t3"
        t3.add_exons([(101, 500), (801, 970), (1100, 1180)])
        t3.add_exons([(101, 500), (801, 970), (1100, 1130)], features="CDS")
        t3.finalize()

        self.t1, self.t2, self.t3 = t1, t2, t3
        self.json_conf = configurator.to_json(None)

    def test_transcript_pickling(self):

        for transcript in [self.t1, self.t2, self.t3]:
            with self.subTest(transcript=transcript):
                pickled = pickle.dumps(transcript)
                unpickled = pickle.loads(pickled)
                self.assertEqual(transcript, unpickled)
                self.assertEqual(len(transcript.combined_cds), len(unpickled.cds_tree))
                self.assertEqual(len(transcript.cds_introntree), len(unpickled.cds_introntree))
                self.assertEqual(len(transcript.segmenttree), len(unpickled.segmenttree))

    def test_locus_unpickling(self):

        for transcript in [self.t1, self.t2, self.t3]:
            for (loc_type, loc_name) in [(_, _.__name__) for _ in (Superlocus, Sublocus, Monosublocus, Locus)]:
                with self.subTest(transcript=transcript, loc_type=loc_type, loc_name=loc_name):
                    loc = loc_type(transcript, json_conf=self.json_conf)
                    pickled = pickle.dumps(transcript)
                    unpickled = pickle.loads(pickled)
                    self.assertEqual(transcript, unpickled)


if __name__ == '__main__':
    unittest.main(verbosity=2)
