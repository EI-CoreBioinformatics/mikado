"""
This module tests the alternative splicing and retained intron functionalities.
"""

import unittest
import Mikado.loci
from Mikado.parsers.GTF import GtfLine
import Mikado.scales
import Mikado.utilities
import logging

__author__ = 'Luca Venturini'


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
    unittest.main()
