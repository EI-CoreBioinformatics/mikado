#!/usr/bin/env python3

# import sys
from collections import OrderedDict
import operator
import mikado_lib
from mikado_lib.loci_objects.transcript_methods import splitting
import unittest
import logging

__author__ = 'Luca Venturini'


class TestSplitMonoexonic(unittest.TestCase):

    logger = mikado_lib.utilities.log_utils.create_null_logger("test_mono")
    logger.setLevel(logging.DEBUG)

    def setUp(self):
        
        self.transcript = mikado_lib.loci_objects.Transcript()
        self.transcript.chrom = "Chr1"
        self.transcript.start = 1001
        self.transcript.end = 6000
        self.transcript.strand = "+"
        self.transcript.id = "transcript1"
        self.transcript.id = "gene1"
        self.transcript.source = "Mikado"
        self.transcript.exons = [(1001, 6000)]
        
        self.bed1 = mikado_lib.parsers.bed12.BED12()
        self.header = False
        self.bed1.chrom = "transcript1"
        self.bed1.start = 1
        self.bed1.end = 5000
        self.bed1.name = "Bed1"
        self.bed1.score = 0
        self.bed1.strand = "+"
        self.bed1.thick_start = 101
        self.bed1.thick_end = 3001
        self.bed1.block_counts = 1
        self.bed1.block_sizes = [2901]
        self.bed1.block_starts = [101]
        self.bed1.transcriptomic = True
        self.bed1.has_start_codon = True
        self.bed1.has_stop_codon = True
        self.assertFalse(self.bed1.invalid, self.bed1.invalid_reason)

        self.bed2 = mikado_lib.parsers.bed12.BED12()
        self.header = False
        self.bed2.chrom = "transcript1"
        self.bed2.start = 1
        self.bed2.end = 5000
        self.bed2.name = "Bed2"
        self.bed2.score = 0
        self.bed2.strand = "+"
        self.bed2.thick_start = 4001
        self.bed2.thick_end = 4900
        self.bed2.block_counts = 1
        self.bed2.block_sizes = [900]
        self.bed2.block_starts = [4001]
        self.bed2.transcriptomic = True
        self.bed2.has_start_codon = True
        self.bed2.has_stop_codon = True
        self.assertFalse(self.bed2.invalid)

        self.transcript.logger = self.logger

        self.transcript.json_conf = dict()
        self.transcript.json_conf["pick"] = dict()
        self.transcript.json_conf["pick"]["chimera_split"] = dict()
        self.transcript.json_conf["pick"]["chimera_split"]["blast_check"] = False

        self.transcript.json_conf["pick"]["chimera_split"]["blast_params"] = dict()
        self.transcript.json_conf["pick"][
            "chimera_split"]["blast_params"]["leniency"] = "LENIENT"
        self.transcript.json_conf["pick"][
            "chimera_split"]["blast_params"]["min_overlap_duplication"] = 0.8
        self.transcript.json_conf["pick"][
            "chimera_split"]["blast_params"]["minimal_hsp_overlap"] = 0.8  # 80%
        self.transcript.json_conf["pick"][
            "chimera_split"]["blast_params"]["hsp_evalue"] = 0.0001
        self.transcript.json_conf["pick"][
            "chimera_split"]["blast_params"]["evalue"] = 0.0001

        self.transcript.json_conf["pick"]["orf_loading"] = dict()
        self.transcript.json_conf["pick"]["orf_loading"]["minimal_secondary_orf_length"] = 50

        self.transcript.load_orfs([self.bed1, self.bed2])
        self.assertTrue(self.transcript.is_coding)
        self.assertEqual(self.transcript.number_internal_orfs,
                         2,
                         str(self.transcript))
        self.transcript.finalize()
        self.assertIn("pick", self.transcript.json_conf)

        self.transcript.logger = self.logger

        # Setup faux hits
        hits = []

        hit1 = dict()
        hit1["evalue"] = 10**-6
        hit1["hsps"] = []
        hit1["target"] = "target1"
        hit1["target_length"] = 3000

        hsp1_1 = dict()
        hsp1_1["hsp_evalue"] = 10**(-6)
        hsp1_1["query_hsp_start"] = 50
        hsp1_1["query_hsp_end"] = 3200
        hsp1_1["target_hsp_start"] = 1
        hsp1_1["target_hsp_end"] = 3000
        hit1["hsps"] = [hsp1_1]

        self.transcript.blast_hits = hits

    def test_simple_split(self):

        """
        Check that the original transcript can be split in two
        :return:
        """

        self.assertIsNotNone(self.transcript.json_conf)
        self.assertIn("pick", self.transcript.json_conf)

        self.transcript.json_conf["pick"]["chimera_split"]["blast_check"] = False

        # with self.assertLogs("test_mono", level="DEBUG") as log_split:
        new_transcripts = list(splitting.split_by_cds(self.transcript))

        # print(log_split.output)

        # self.assertIn("DEBUG:test_mono:",
        #               log_split.output)
        self.assertEqual(len(new_transcripts), 2, "\n".join(str(_) for _ in new_transcripts))
        self.assertEqual(new_transcripts[0].start, self.transcript.start)
        self.assertEqual(new_transcripts[1].end, self.transcript.end)

    def test_lenient_split(self):

        # TODO: this test is failing, it is of primary importance to understand why

        self.transcript.json_conf["pick"]["chimera_split"]["blast_check"] = True

        self.transcript.json_conf["pick"][
            "chimera_split"]["blast_params"]["leniency"] = "LENIENT"

        cds_boundaries = OrderedDict()
        for orf in sorted(self.transcript.loaded_bed12,
                          key=operator.attrgetter("thick_start", "thick_end")):
            cds_boundaries[(orf.thick_start, orf.thick_end)] = [orf]

        self.assertEqual(1,
                         len(splitting.check_split_by_blast(self.transcript, cds_boundaries)))

    def test_stringent_split(self):

        self.transcript.json_conf["pick"]["chimera_split"]["blast_check"] = True

        self.transcript.json_conf["pick"][
            "chimera_split"]["blast_params"]["leniency"] = "STRINGENT"

        cds_boundaries = OrderedDict()
        for orf in sorted(self.transcript.loaded_bed12,
                          key=operator.attrgetter("thick_start", "thick_end")):
            cds_boundaries[(orf.thick_start, orf.thick_end)] = [orf]

        self.assertEqual(1,
                         len(splitting.check_split_by_blast(self.transcript, cds_boundaries)))

    def test_permissive_split(self):

        self.transcript.logger = self.logger

        self.transcript.json_conf["pick"]["chimera_split"]["blast_check"] = True

        self.transcript.json_conf["pick"][
            "chimera_split"]["blast_params"]["leniency"] = "PERMISSIVE"

        cds_boundaries = OrderedDict()
        for orf in sorted(self.transcript.loaded_bed12,
                          key=operator.attrgetter("thick_start", "thick_end")):
            cds_boundaries[(orf.thick_start, orf.thick_end)] = [orf]

        self.assertEqual(2,
                         len(splitting.check_split_by_blast(self.transcript, cds_boundaries)))



if __name__ == "__main__":
    unittest.main()
