#!/usr/bin/env python3

import logging
import operator
import unittest
from sys import version_info

import Mikado
from Mikado.transcripts.transcript_methods import splitting

if version_info.minor < 5:
    from sortedcontainers import SortedDict
else:
    from collections import OrderedDict as SortedDict

__author__ = 'Luca Venturini'


class TestSplitMonoexonic(unittest.TestCase):

    logger = Mikado.utilities.log_utils.create_null_logger("test_mono")
    logger.setLevel(logging.DEBUG)

    def setUp(self):
        
        self.transcript = Mikado.loci.Transcript()
        self.transcript.chrom = "Chr1"
        self.transcript.start = 1001
        self.transcript.end = 6000
        self.transcript.strand = "+"
        self.transcript.id = "transcript1"
        self.transcript.parent = "gene1"
        self.transcript.source = "Mikado"
        self.transcript.exons = [(1001, 6000)]
        
        self.bed1 = Mikado.parsers.bed12.BED12()
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

        self.bed2 = Mikado.parsers.bed12.BED12()
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
        hit1["query_start"] = 50
        hit1["query_end"] = 3200
        hit1["evalue"] = 10**(-6)
        hit1["global_identity"] = 100
        hit1["query_aligned_length"] = 3200 - 500
        hit1["target_aligned_length"] = 3000 - 1

        hsp1_1 = dict()
        hsp1_1["hsp_evalue"] = 10**(-6)
        hsp1_1["hsp_bits"] = 1200
        hsp1_1["query_hsp_start"] = 50
        hsp1_1["query_hsp_end"] = 3200
        hsp1_1["target_hsp_start"] = 1
        hsp1_1["target_hsp_end"] = 3000
        hsp1_1["match"] = "|" * (3200 - 500)

        hit1["hsps"] = [hsp1_1]

        hits.append(hit1)

        self.transcript.blast_hits = hits
        self.assertEqual(len(self.transcript.blast_hits), 1)

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

        self.transcript.json_conf["pick"]["chimera_split"]["blast_check"] = True

        self.transcript.json_conf["pick"][
            "chimera_split"]["blast_params"]["leniency"] = "LENIENT"

        cds_boundaries = SortedDict()
        for orf in sorted(self.transcript.loaded_bed12,
                          key=operator.attrgetter("thick_start", "thick_end")):
            cds_boundaries[(orf.thick_start, orf.thick_end)] = [orf]

        self.assertEqual(1,
                         len(splitting.check_split_by_blast(self.transcript, cds_boundaries)))

    def test_stringent_split(self):

        self.transcript.json_conf["pick"]["chimera_split"]["blast_check"] = True

        self.transcript.json_conf["pick"][
            "chimera_split"]["blast_params"]["leniency"] = "STRINGENT"

        cds_boundaries = SortedDict()
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

        cds_boundaries = SortedDict()
        for orf in sorted(self.transcript.loaded_bed12,
                          key=operator.attrgetter("thick_start", "thick_end")):
            cds_boundaries[(orf.thick_start, orf.thick_end)] = [orf]

        self.assertEqual(2,
                         len(splitting.check_split_by_blast(self.transcript, cds_boundaries)))

    @staticmethod
    def get_second_hit():

        hit1 = dict()
        hit1["evalue"] = 10**-6
        hit1["hsps"] = []
        hit1["target"] = "target2"
        hit1["target_length"] = 2000
        hit1["query_start"] = 3700
        hit1["query_end"] = 4850
        hit1["evalue"] = 10**(-6)
        hit1["global_identity"] = 100
        hit1["query_aligned_length"] = 4850 - 3700
        hit1["target_aligned_length"] = 1149

        hsp1_1 = dict()
        hsp1_1["hsp_evalue"] = 10**(-6)
        hsp1_1["hsp_bits"] = 1200
        hsp1_1["query_hsp_start"] = 3700
        hsp1_1["query_hsp_end"] = 4850
        hsp1_1["target_hsp_start"] = 1
        hsp1_1["target_hsp_end"] = 2000
        hsp1_1["match"] = "|" * (4850 - 3700)

        hit1["hsps"] = [hsp1_1]
        return hit1

    def test_permissive_split_twohits(self):

        hit2 = self.get_second_hit()
        self.transcript.blast_hits.append(hit2)
        self.transcript.json_conf["pick"]["chimera_split"]["blast_check"] = True

        self.transcript.json_conf["pick"][
            "chimera_split"]["blast_params"]["leniency"] = "PERMISSIVE"

        self.assertEqual(2,
                         len(list(splitting.split_by_cds(self.transcript))))

    def test_stringent_split_twohits(self):

        hit2 = self.get_second_hit()
        self.transcript.blast_hits.append(hit2)
        self.transcript.json_conf["pick"]["chimera_split"]["blast_check"] = True

        self.transcript.json_conf["pick"][
            "chimera_split"]["blast_params"]["leniency"] = "STRINGENT"

        self.assertEqual(2,
                         len(list(splitting.split_by_cds(self.transcript))))

    def test_lenient_split_twohits(self):

        hit2 = self.get_second_hit()
        self.transcript.blast_hits.append(hit2)
        self.transcript.json_conf["pick"]["chimera_split"]["blast_check"] = True

        self.transcript.json_conf["pick"][
            "chimera_split"]["blast_params"]["leniency"] = "LENIENT"

        self.assertEqual(2,
                         len(list(splitting.split_by_cds(self.transcript))))

    @staticmethod
    def get_spanning_hit():

        hit1 = dict()
        hit1["evalue"] = 10**-6
        hit1["hsps"] = []
        hit1["target"] = "target1"
        hit1["target_length"] = 5000
        hit1["query_start"] = 51
        hit1["query_end"] = 5000
        hit1["evalue"] = 10**(-6)
        hit1["global_identity"] = 100
        hit1["query_aligned_length"] = 4950
        hit1["target_aligned_length"] = 4950

        hsp1_1 = dict()
        hsp1_1["hsp_evalue"] = 10**(-6)
        hsp1_1["hsp_bits"] = 1200
        hsp1_1["query_hsp_start"] = 51
        hsp1_1["query_hsp_end"] = 5000
        hsp1_1["target_hsp_start"] = 1
        hsp1_1["target_hsp_end"] = 4950
        hsp1_1["match"] = "|" * 4950

        hit1["hsps"] = [hsp1_1]
        return hit1

    def test_spanning_hit_lenient(self):

        self.transcript.blast_hits = [self.get_spanning_hit()]
        self.transcript.json_conf["pick"]["chimera_split"]["blast_check"] = True

        self.transcript.json_conf["pick"][
            "chimera_split"]["blast_params"]["leniency"] = "LENIENT"

        self.assertEqual(1,
                         len(list(splitting.split_by_cds(self.transcript))))

    def test_spanning_hit_nocheck(self):
        self.transcript.blast_hits = [self.get_spanning_hit()]
        self.transcript.json_conf["pick"]["chimera_split"]["blast_check"] = False
        cds_boundaries = SortedDict()
        for orf in sorted(self.transcript.loaded_bed12,
                          key=operator.attrgetter("thick_start", "thick_end")):
            cds_boundaries[(orf.thick_start, orf.thick_end)] = [orf]

        self.assertEqual(len(cds_boundaries), 2)
        self.assertEqual(self.transcript.number_internal_orfs, 2)
        self.assertEqual(2,
                         len(list(splitting.split_by_cds(self.transcript))))

    def test_deleted_hits(self):

        delattr(self.transcript, "blast_hits")
        self.transcript.json_conf["pick"]["chimera_split"]["blast_check"] = True
        self.transcript.json_conf["pick"][
            "chimera_split"]["blast_params"]["leniency"] = "LENIENT"
        self.transcript.logger = self.logger
        with self.assertLogs("null", level="WARNING") as log_split:
            self.assertEqual(2,
                             len(list(splitting.split_by_cds(self.transcript))))

        self.assertIn("WARNING:null:BLAST hits store lost for transcript1! Creating a mock one to avoid a crash",
                      log_split.output)

    def test_duplication(self):

        pass


if __name__ == "__main__":
    unittest.main()
