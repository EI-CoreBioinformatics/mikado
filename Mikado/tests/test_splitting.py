#!/usr/bin/env python3

import logging
import operator
import unittest
from sys import version_info
import pysam
from .. import loci, parsers, utilities, configuration
from ..transcripts.transcript_methods import splitting
import tempfile
from ..parsers import bed12
from ..utilities.log_utils import create_default_logger
if version_info.minor < 5:
    from sortedcontainers import SortedDict
else:
    from collections import OrderedDict as SortedDict

__author__ = 'Luca Venturini'


class TestSplitMonoexonic(unittest.TestCase):

    logger = utilities.log_utils.create_null_logger("test_mono")
    logger.setLevel(logging.DEBUG)

    def setUp(self):
        
        self.transcript = loci.Transcript()
        self.transcript.chrom = "Chr1"
        self.transcript.start = 1001
        self.transcript.end = 6000
        self.transcript.strand = "+"
        self.transcript.id = "transcript1"
        self.transcript.parent = "gene1"
        self.transcript.source = "Mikado"
        self.transcript.exons = [(1001, 6000)]
        
        self.bed1 = parsers.bed12.BED12()
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

        self.bed2 = parsers.bed12.BED12()
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

        self.transcript.json_conf = configuration.configurator.to_json(None)
        self.transcript.json_conf["pick"]["chimera_split"]["blast_check"] = False
        self.transcript.json_conf["pick"][
            "chimera_split"]["blast_params"]["leniency"] = "LENIENT"
        self.transcript.json_conf["pick"]["chimera_split"]["blast_params"]["min_overlap_duplication"] = 0.8
        self.transcript.json_conf["pick"]["chimera_split"]["blast_params"]["minimal_hsp_overlap"] = 0.8  # 80%
        self.transcript.json_conf["pick"]["chimera_split"]["blast_params"]["hsp_evalue"] = 0.0001
        self.transcript.json_conf["pick"]["chimera_split"]["blast_params"]["evalue"] = 0.0001
        self.transcript.json_conf["pick"]["orf_loading"]["minimal_secondary_orf_length"] = 50

        self.transcript.load_orfs([self.bed1, self.bed2])
        self.assertTrue(self.transcript.is_coding)
        self.assertEqual(self.transcript.number_internal_orfs,
                         2, str(self.transcript))
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

        self.assertEqual(len(new_transcripts), 2, "\n".join(str(_) for _ in new_transcripts))
        self.assertEqual(new_transcripts[0].start, self.transcript.start)
        self.assertEqual(new_transcripts[1].end, self.transcript.end)

        sl = loci.Superlocus(self.transcript, json_conf=self.transcript.json_conf)
        self.assertFalse(sl.json_conf["pick"]["chimera_split"]["blast_check"])
        self.assertEqual(len(sl.transcripts), 1)
        sl.logger.setLevel("DEBUG")
        sl.load_all_transcript_data(data_dict=dict())
        self.assertEqual(len(sl.transcripts), 2)

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
        sl = loci.Superlocus(self.transcript, json_conf=self.transcript.json_conf)
        self.assertEqual(len(sl.transcripts), 1)
        sl.load_all_transcript_data(data_dict=dict())
        self.assertEqual(len(sl.transcripts), 1)

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
        sl = loci.Superlocus(self.transcript, json_conf=self.transcript.json_conf)
        self.assertEqual(len(sl.transcripts), 1)
        sl.load_all_transcript_data(data_dict=dict())
        self.assertEqual(len(sl.transcripts), 1)

    def test_permissive_split(self):

        self.transcript.logger = self.logger

        self.transcript.json_conf["pick"]["chimera_split"]["blast_check"] = True

        self.transcript.json_conf["pick"][
            "chimera_split"]["blast_params"]["leniency"] = "PERMISSIVE"

        cds_boundaries = SortedDict()
        for orf in sorted(self.transcript.loaded_bed12,
                          key=operator.attrgetter("thick_start", "thick_end")):
            cds_boundaries[(orf.thick_start, orf.thick_end)] = [orf]

        self.assertEqual(2, len(splitting.check_split_by_blast(self.transcript, cds_boundaries)))
        sl = loci.Superlocus(self.transcript, json_conf=self.transcript.json_conf)
        self.assertEqual(len(sl.transcripts), 1)
        sl.load_all_transcript_data(data_dict=dict())
        self.assertEqual(len(sl.transcripts), 2)

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

        self.assertEqual(2, len(list(splitting.split_by_cds(self.transcript))))
        sl = loci.Superlocus(self.transcript, json_conf=self.transcript.json_conf)
        self.assertEqual(len(sl.transcripts), 1)
        sl.load_all_transcript_data(data_dict=dict())
        self.assertEqual(len(sl.transcripts), 2)

    def test_stringent_split_twohits(self):

        hit2 = self.get_second_hit()
        self.transcript.blast_hits.append(hit2)
        self.transcript.json_conf["pick"]["chimera_split"]["blast_check"] = True

        self.transcript.json_conf["pick"][
            "chimera_split"]["blast_params"]["leniency"] = "STRINGENT"

        self.assertEqual(2, len(list(splitting.split_by_cds(self.transcript))))
        sl = loci.Superlocus(self.transcript, json_conf=self.transcript.json_conf)
        self.assertEqual(len(sl.transcripts), 1)
        sl.load_all_transcript_data(data_dict=dict())
        self.assertEqual(len(sl.transcripts), 2)

    def test_no_splitting_by_source(self):
        self.transcript.source = "foo"
        for sources in [[], [self.transcript.source], ["bar"], ["bar", [self.transcript.source]]]:
            with self.subTest(sources=sources):
                self.transcript.json_conf["pick"]["chimera_split"]["skip"] = sources
                if self.transcript.source in sources:
                    final = 1
                else:
                    final = 2
                self.assertEqual(final, len(list(splitting.split_by_cds(self.transcript))))
                sl = loci.Superlocus(self.transcript, json_conf=self.transcript.json_conf)
                sl.load_all_transcript_data(data_dict=dict())
                self.assertEqual(len(sl.transcripts), final)

    def test_one_orf(self):
        self.transcript.strip_cds()
        self.transcript.load_orfs([self.bed1])
        self.assertEqual(1, len(list(splitting.split_by_cds(self.transcript))))

    def test_no_hsps(self):
        self.transcript.blast_hits = []
        self.transcript.json_conf["pick"]["chimera_split"]["blast_check"] = True

        self.transcript.json_conf["pick"]["chimera_split"]["blast_params"]["leniency"] = "LENIENT"

        logger = utilities.log_utils.create_default_logger("test_no_hsps")
        logger.setLevel("DEBUG")
        self.transcript.logger = logger
        self.assertEqual(2, len(list(splitting.split_by_cds(self.transcript))))
        sl = loci.Superlocus(self.transcript, json_conf=self.transcript.json_conf)
        self.assertEqual(len(sl.transcripts), 1)
        sl.load_all_transcript_data(data_dict=dict())
        self.assertEqual(len(sl.transcripts), 2)

    def test_lenient_split_twohits(self):

        hit2 = self.get_second_hit()
        self.transcript.blast_hits.append(hit2)
        self.transcript.json_conf["pick"]["chimera_split"]["blast_check"] = True

        self.transcript.json_conf["pick"][
            "chimera_split"]["blast_params"]["leniency"] = "LENIENT"

        self.assertEqual(2, len(list(splitting.split_by_cds(self.transcript))))
        sl = loci.Superlocus(self.transcript, json_conf=self.transcript.json_conf)
        self.assertEqual(len(sl.transcripts), 1)
        sl.load_all_transcript_data(data_dict=dict())
        self.assertEqual(len(sl.transcripts), 2)

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

        self.assertEqual(1, len(list(splitting.split_by_cds(self.transcript))))
        sl = loci.Superlocus(self.transcript, json_conf=self.transcript.json_conf)
        self.assertEqual(len(sl.transcripts), 1)
        sl.load_all_transcript_data(data_dict=dict())
        self.assertEqual(len(sl.transcripts), 1)

    def test_spanning_hit_nocheck(self):
        self.transcript.blast_hits = [self.get_spanning_hit()]
        self.transcript.json_conf["pick"]["chimera_split"]["blast_check"] = False
        cds_boundaries = SortedDict()
        for orf in sorted(self.transcript.loaded_bed12,
                          key=operator.attrgetter("thick_start", "thick_end")):
            cds_boundaries[(orf.thick_start, orf.thick_end)] = [orf]

        self.assertEqual(len(cds_boundaries), 2)
        self.assertEqual(self.transcript.number_internal_orfs, 2)
        self.assertEqual(2, len(list(splitting.split_by_cds(self.transcript))))
        sl = loci.Superlocus(self.transcript, json_conf=self.transcript.json_conf)
        self.assertEqual(len(sl.transcripts), 1)
        sl.load_all_transcript_data(data_dict=dict())
        self.assertEqual(len(sl.transcripts), 2)

    def test_deleted_hits(self):

        delattr(self.transcript, "blast_hits")
        self.transcript.json_conf["pick"]["chimera_split"]["blast_check"] = True
        self.transcript.json_conf["pick"][
            "chimera_split"]["blast_params"]["leniency"] = "LENIENT"
        self.transcript.logger = self.logger
        with self.assertLogs(logger=self.logger, level="WARNING") as log_split:
            self.assertEqual(2, len(list(splitting.split_by_cds(self.transcript))))

        self.assertIn("WARNING:test_mono:BLAST hits store lost for transcript1! Creating a mock one to avoid a crash",
                      log_split.output)

    def test_duplication(self):

        pass


class TestWithPhase(unittest.TestCase):

    logger = utilities.log_utils.create_default_logger("test_mono")
    logger.setLevel(logging.INFO)

    def setUp(self):

        self.transcript = loci.Transcript()
        self.transcript.chrom = "Chr1"
        self.transcript.start = 1001
        self.transcript.end = 6000
        self.transcript.strand = "."
        self.transcript.id = "transcript1"
        self.transcript.parent = "gene1"
        self.transcript.source = "Mikado"
        self.transcript.exons = [(1001, 6000)]
        self.transcript.logger = self.logger
        self.transcript.logger.setLevel("INFO")

        self.transcript.json_conf = configuration.configurator.to_json("")
        self.transcript.json_conf["pick"]["chimera_split"]["blast_check"] = False
        self.assertIsNotNone(self.transcript.json_conf)
        self.assertIn("pick", self.transcript.json_conf)

    def testPositive(self):

        self.bed1 = parsers.bed12.BED12()
        self.header = False
        self.bed1.chrom = "transcript1"
        self.bed1.start = 1
        self.bed1.end = 5000
        self.bed1.name = "Bed1"
        self.bed1.score = 0
        self.bed1.strand = "+"
        self.bed1.thick_start = 1
        self.bed1.thick_end = 3002
        self.bed1.phase = 2
        self.bed1.block_counts = 1
        self.bed1.block_sizes = [3002]
        self.bed1.block_starts = [1]
        self.bed1.transcriptomic = True
        self.bed1.has_start_codon = False
        self.bed1.has_stop_codon = True
        self.assertFalse(self.bed1.invalid, self.bed1.invalid_reason)

        self.bed2 = parsers.bed12.BED12()
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

        self.transcript.load_orfs([self.bed1, self.bed2])
        self.assertTrue(self.transcript.is_coding)
        self.assertEqual(self.transcript.number_internal_orfs,
                         2, str(self.transcript))
        self.transcript.finalize()

        # with self.assertLogs("test_mono", level="DEBUG") as log_split:
        new_transcripts = list(splitting.split_by_cds(self.transcript))

        # self.assertIn("DEBUG:test_mono:",
        #               log_split.output)
        self.assertEqual(len(new_transcripts), 2, "\n".join(str(_) for _ in new_transcripts))
        self.assertEqual(new_transcripts[0].start, self.transcript.start)
        self.assertEqual(new_transcripts[1].end, self.transcript.end)

    def testNegative(self):
        self.bed1 = parsers.bed12.BED12()
        self.header = False
        self.bed1.chrom = "transcript1"
        self.bed1.start = 1
        self.bed1.end = 5000
        self.bed1.name = "Bed1"
        self.bed1.score = 0
        self.bed1.strand = "-"
        self.bed1.thick_start = 1
        self.bed1.thick_end = 3000
        self.bed1.phase = 0
        self.bed1.block_counts = 1
        self.bed1.block_sizes = [3001]
        self.bed1.block_starts = [1]
        self.bed1.transcriptomic = True
        self.bed1.has_start_codon = True
        self.bed1.has_stop_codon = True
        self.assertFalse(self.bed1.invalid, self.bed1.invalid_reason)

        logger = create_default_logger("testNegative", "DEBUG")
        self.bed2 = parsers.bed12.BED12(logger=logger)
        self.header = False
        self.bed2.chrom = "transcript1"
        self.bed2.start = 1
        self.bed2.end = 5000
        self.bed2.name = "Bed2"
        self.bed2.score = 0
        self.bed2.strand = "-"
        self.bed2.thick_start = 4001
        self.bed2.thick_end = 5000
        self.bed2.phase = 1
        self.bed2.block_counts = 1
        self.bed2.block_sizes = [1000]
        self.bed2.block_starts = [4001]
        self.bed2.transcriptomic = True
        self.bed2.has_start_codon = False
        self.bed2.has_stop_codon = True
        self.assertFalse(self.bed2.invalid, (self.bed2.phase, self.bed2.invalid_reason))

        self.transcript.load_orfs([self.bed1, self.bed2])
        self.assertTrue(self.transcript.is_coding)
        self.assertEqual(self.transcript.number_internal_orfs,
                         2,
                         str(self.transcript))
        self.transcript.finalize()

        self.assertEqual(self.transcript.number_internal_orfs,
                         2,
                         str(self.transcript))

        # with self.assertLogs("test_mono", level="DEBUG") as log_split:
        new_transcripts = list(splitting.split_by_cds(self.transcript))

        self.assertEqual(len(new_transcripts), 2, "\n".join(str(_) for _ in new_transcripts))
        self.assertEqual(new_transcripts[0].start, self.transcript.start)
        self.assertEqual(new_transcripts[1].end, self.transcript.end)

        self.assertEqual(self.transcript.combined_cds_start, 6000)
        self.assertEqual(self.transcript.combined_cds_end, 1001)
        self.assertEqual(self.transcript.selected_cds_end, 1001)
        self.assertEqual(self.transcript.selected_cds_start, 4000)
        self.assertEqual(self.transcript.strand, "-")

        sl = loci.Superlocus(self.transcript, json_conf=self.transcript.json_conf)
        self.assertEqual(len(sl.transcripts), 1)
        sl.load_all_transcript_data(data_dict=dict())
        self.assertEqual(len(sl.transcripts), 2)

    def test_negative_orf_gtg(self):

        fasta = """>sex_morph_FW.stringtie_sex_morph_FW_str.232.2
CACAGTCTCGTGCGGCTATTTTCGTCCGCCGCCTGTCCCTCTAAGAAGAGTTTAAGCTCC
TGGGAGCCGGCGGTAGCCCTAGTAACGTATCGTGATCCGCCGGCGACGGCCGCGAACGCG
GCGCGCTTCACCACGGAGCCCCAGGCACGACAGACGCACGCCGACCAGGGGATGAACCCC
GGCCGACGCACGCCCGCCGCACGCAGGGACGCGGATGGCGCGGCCGCGCCCGACGACCGC
CGTGGACGACGGGCGAACGCGTTCGGGGATACCGGGCCGAGCCGACGGGAACGCGAACAC
GGACGGCCGAAACCGCCCGCGCCGCGCCCACCGCCGACCCGGGTTTACCCGCCTAGTTAG
CAGGACAGAGTCTCGTTCGTTATCGGAATTAACCAGACAGATCGCTCCACCAACTAAGAA
CGGCCATGCACCACCACCCACCGAATCAAGAAAGAGCTCTCAATCTGTCAATCTTTCCGG
TGTCCGGGCCTGGTGAGGTTTCCCGTGTTGAGTCAAATTAAGCCGCAGGCTCCACTCCTG
GTGGTGCCCTTCCGTCAATTCCTTTAAGTTTCAACTTTGCAATCATACTTCCCCCGGAAC
CGAAAAGCTTCGGTTTCCCGGAAGCTGCCCGCCGGGTCGTTAATGAAACGCCGGCGGATC
GCTAGCTGGCATCGTTTACAGTTAGAACTAGGGCGGTATCTGATCGCCTTCGAACCTCTA
ACTTTCGTTCTTGATCATACGAGAACGTACTTGGCAAATGCTTTCGCGTCAGTTCGTCTC
GAGACGATCCAAGAATTTCACCTCTAACGTCTCGGTACGAATGCCCCCGCCCGTCTCTGT
TGATCATTACCCCGGAGGGCGATTTCGCGCGCCCGCGAAGGGCGGAGATGCGCGGGACCA
AGGTCTTGTTCCATTATTCCATGCGACCAGTATTCAGGGCCTTTTGACGAGACGGCCGTG
AAGCCGCCCCGCCAGATTTGAGCCTGCTTTGAGCACTCTAATTTGTTCAAAGTAAACGTG
TCGGCCCGCCGACGGCACTCGGTGAAGAGCACCGCGCAGCAAGATTGGAGTAGGCGGCCG
CCGTCGTCGAACCCCGACGGCCGCGCGACGCGTGGCCGCGCGGCGCGCCGGAAGCACGAG
ACACGTGTCCGCC"""

        prodigal = "\t".join(["sex_morph_FW.stringtie_sex_morph_FW_str.232.2",
                              "Prodigal_v2.6.3",
                              "CDS", "934","1137", "5.7", "-", "0",
                              "ID=12199_3;partial=00;start_type=GTG;rbs_motif=None;rbs_spacer=None;gc_cont=0.657;conf=78.73;score=5.69;cscore=7.02;sscore=-1.33;rscore=0.33;uscore=-0.91;tscore=-0.75;"
                              ])
        temp_fa = tempfile.NamedTemporaryFile(mode="wt", suffix=".fa")
        temp_gff = tempfile.NamedTemporaryFile(mode="wt", suffix=".gff3")
        print("##gff-version\t3", file=temp_gff)
        print(prodigal, file=temp_gff)
        print(fasta, file=temp_fa)
        temp_gff.flush()
        temp_fa.flush()
        handle = open(temp_gff.name, mode="rt")
        fasta_index = pysam.FastaFile(temp_fa.name)
        bed12_parser = bed12.Bed12Parser(handle,
                                         fasta_index=fasta_index,
                                         is_gff=True,
                                         transcriptomic=True,
                                         max_regression=1)
        record = next(bed12_parser)
        self.assertIsInstance(record, bed12.BED12)
        self.assertFalse(record.has_start_codon)
        self.assertFalse(record.invalid)
        self.assertEqual(record.phase, 1, record)


if __name__ == "__main__":
    unittest.main()
