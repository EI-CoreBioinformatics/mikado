"""
This module tests the alternative splicing and retained intron functionalities.
"""

import unittest
import mikado_lib.loci_objects
from mikado_lib.parsers.GTF import GtfLine
import mikado_lib.scales
import mikado_lib.utilities
import logging

__author__ = 'Luca Venturini'


class TestLocus(unittest.TestCase):

    """
    This unit test is focussed on the locus definition and alternative splicings.
    """

    logger = mikado_lib.utilities.log_utils.create_default_logger("tester")

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
        self.json_conf["pick"]["alternative_splicing"]["valid_ccodes"] = ["j", "n", "O", "h"]

        t1 = """Chr1\tfoo\ttranscript\t1000\t3000\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.1";
        Chr1\tfoo\texon\t1000\t1300\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.1";
        Chr1\tfoo\tCDS\t1100\t1300\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.1";
        Chr1\tfoo\texon\t1700\t2000\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.1";
        Chr1\tfoo\tCDS\t1700\t2000\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.1";
        Chr1\tfoo\texon\t2100\t2500\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.1";
        Chr1\tfoo\tCDS\t2100\t2500\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.1";
        Chr1\tfoo\texon\t2800\t3000\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.1";
        Chr1\tfoo\tCDS\t2800\t2900\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.1";
        """

        t1lines = [GtfLine(line) for line in t1.split("\n") if line]
        self.t1 = mikado_lib.loci_objects.Transcript(t1lines[0])
        for exon in t1lines[1:]:
            if exon.header:
                continue
            self.t1.add_exon(exon)
        self.t1.finalize()

        # Just a fragment of the best
        t1_contained = """Chr1\tfoo\ttranscript\t1000\t2000\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.2";
        Chr1\tfoo\texon\t1000\t1300\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.2";
        Chr1\tfoo\tCDS\t1100\t1300\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.2";
        Chr1\tfoo\texon\t1700\t2000\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.2";
        Chr1\tfoo\tCDS\t1700\t1900\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.2";
        """

        t1_contained_lines = [GtfLine(line) for line in t1_contained.split("\n") if line]
        self.t1_contained = mikado_lib.loci_objects.Transcript(t1_contained_lines[0])
        for exon in t1_contained_lines[1:]:
            if exon.header:
                continue
            self.t1_contained.add_exon(exon)
        self.t1_contained.finalize()

        # Valid AS event
        t1_as = """Chr1\tfoo\ttranscript\t1000\t3000\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.3";
        Chr1\tfoo\texon\t1000\t1300\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.3";
        Chr1\tfoo\tCDS\t1100\t1300\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.3";
        Chr1\tfoo\texon\t1700\t2000\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.3";
        Chr1\tfoo\tCDS\t1700\t2000\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.3";
        Chr1\tfoo\texon\t2100\t2400\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.3";
        Chr1\tfoo\tCDS\t2100\t2400\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.3";
        Chr1\tfoo\texon\t2800\t3000\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.3";
        Chr1\tfoo\tCDS\t2800\t2900\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.3";
        """

        t1_as_lines = [GtfLine(line) for line in t1_as.split("\n") if line]
        self.t1_as = mikado_lib.loci_objects.Transcript(t1_as_lines[0])
        for exon in t1_as_lines[1:]:
            if exon.header:
                continue
            self.t1_as.add_exon(exon)
        self.t1_as.finalize()

        # Retained intron AS event
        t1_retained = """Chr1\tfoo\ttranscript\t1000\t2800\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.4";
        Chr1\tfoo\texon\t1000\t1300\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.4";
        Chr1\tfoo\tCDS\t1100\t1300\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.4";
        Chr1\tfoo\texon\t1700\t2000\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.4";
        Chr1\tfoo\tCDS\t1700\t2000\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.4";
        Chr1\tfoo\texon\t2100\t2800\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.4";
        Chr1\tfoo\tCDS\t2100\t2450\t.\t+\t.\tgene_id "Chr1.1"; transcript_id "Chr1.1.4";
        """

        t1_retained_lines = [GtfLine(line) for line in t1_retained.split("\n") if line]
        self.t1_retained = mikado_lib.loci_objects.Transcript(t1_retained_lines[0])
        for exon in t1_retained_lines[1:]:
            if exon.header:
                continue
            self.t1_retained.add_exon(exon)
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
        result, _ = mikado_lib.scales.assigner.Assigner.compare(self.t1_contained, self.t1)
        self.assertEqual(result.ccode[0], "c")

        # The valid AS should have a j assigned
        result, _ = mikado_lib.scales.assigner.Assigner.compare(self.t1_as, self.t1)
        self.assertEqual(result.ccode[0], "j")

        # The retained intron AS should have a j assigned
        result, _ = mikado_lib.scales.assigner.Assigner.compare(self.t1_retained, self.t1)
        self.assertEqual(result.ccode[0], "j", result.ccode)

    def testCreate(self):

        """
        Test the creation of the locus
        :return:
        """

        locus = mikado_lib.loci_objects.Locus(self.t1, logger=self.logger)
        locus.json_conf = self.json_conf
        self.assertEqual(len(locus.transcripts), 1)

    def test_exclude_contained(self):

        """Test that we exclude a transcript with a contained class code (c)"""

        locus = mikado_lib.loci_objects.Locus(self.t1, logger=self.logger)
        locus.json_conf = self.json_conf
        self.assertEqual(len(locus.transcripts), 1)
        locus.add_transcript_to_locus(self.t1_contained)
        self.assertEqual(len(locus.transcripts), 1)

    def test_add_contained(self):

        """Test that we add a transcript with a contained class code (c) if
        we explicitly ask for it"""

        locus = mikado_lib.loci_objects.Locus(self.t1, logger=self.logger)
        locus.json_conf = self.json_conf
        locus.json_conf["pick"]["alternative_splicing"]["valid_ccodes"].append("c")
        self.assertEqual(len(locus.transcripts), 1)
        locus.add_transcript_to_locus(self.t1_contained)
        self.assertEqual(len(locus.transcripts), 2)

    def test_addValid(self):

        """Test that we can successfully add a transcript to the locus if
        it passes the muster."""

        locus = mikado_lib.loci_objects.Locus(self.t1, logger=self.logger)
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

        locus = mikado_lib.loci_objects.Locus(self.t1, logger=self.logger)
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

        locus = mikado_lib.loci_objects.Locus(self.t1, logger=self.logger)
        locus.json_conf = self.json_conf
        self.assertEqual(len(locus.transcripts), 1)
        locus.add_transcript_to_locus(self.t1_retained)
        self.assertEqual(len(locus.transcripts), 1)

    def test_add_with_retained_intron(self):
        """Test that a transcript with a retained intron is kept as valid
        if we change the switch."""

        locus = mikado_lib.loci_objects.Locus(self.t1, logger=self.logger)
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
        logger.setLevel(logging.DEBUG)
        locus = mikado_lib.loci_objects.Locus(self.t1, logger=logger)
        locus.json_conf = self.json_conf
        self.assertEqual(len(locus.transcripts), 1)
        locus.add_transcript_to_locus(candidate)
        self.assertEqual(len(locus.transcripts), 1)


unittest.main()
