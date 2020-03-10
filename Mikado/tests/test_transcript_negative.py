# coding: utf-8

"""
Unit test for a transcript on the negative strand.
"""

import unittest
import operator
import re
import logging
from .. import parsers, loci, exceptions
from ..utilities.log_utils import create_null_logger  # , create_default_logger


class TranscriptTesterNegative(unittest.TestCase):

    logger = create_null_logger("null")
    logger.setLevel(logging.WARNING)

    tr_gff = """Chr1    TAIR10    mRNA    5928    8737    .    -    .    ID=AT1G01020.1;Parent=AT1G01020
Chr1    TAIR10    five_prime_UTR    8667    8737    .    -    .    Parent=AT1G01020.1
Chr1    TAIR10    CDS    8571    8666    .    -    0    Parent=AT1G01020.1;
Chr1    TAIR10    exon    8571    8737    .    -    .    Parent=AT1G01020.1
Chr1    TAIR10    CDS    8417    8464    .    -    0    Parent=AT1G01020.1;
Chr1    TAIR10    exon    8417    8464    .    -    .    Parent=AT1G01020.1
Chr1    TAIR10    CDS    8236    8325    .    -    0    Parent=AT1G01020.1;
Chr1    TAIR10    exon    8236    8325    .    -    .    Parent=AT1G01020.1
Chr1    TAIR10    CDS    7942    7987    .    -    0    Parent=AT1G01020.1;
Chr1    TAIR10    exon    7942    7987    .    -    .    Parent=AT1G01020.1
Chr1    TAIR10    CDS    7762    7835    .    -    2    Parent=AT1G01020.1;
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

    tr_lines = [line for line in tr_gff.split("\n") if line]
    for pos, line in enumerate(tr_lines):
        tr_lines[pos] = re.sub(r"\s+", "\t", line)
        assert len(tr_lines[pos].split("\t")) == 9, line.split("\t")

    tr_gff_lines = [parsers.GFF.GffLine(line) for line in tr_lines]

    for l in tr_gff_lines:
        assert l.header is False

    def setUp(self):
        """Basic creation test."""

        self.tr = loci.Transcript(self.tr_gff_lines[0], logger=self.logger)
        for line in self.tr_gff_lines[1:]:
            self.tr.add_exon(line)
        self.tr.name = self.tr.id
        self.tr.finalize()
        self.tr.logger = self.logger

        self.orf = parsers.bed12.BED12()
        self.orf.chrom = self.tr.id
        self.orf.start = 1
        self.orf.end = self.tr.cdna_length
        self.orf.name = self.tr.id
        self.orf.strand = "+"
        self.orf.score = 0
        self.orf.thick_start = self.tr.selected_start_distance_from_tss + 1
        self.orf.thick_end = self.tr.cdna_length - self.tr.selected_end_distance_from_tes
        self.orf.block_count = 1
        self.orf.blockSize = self.tr.cdna_length
        self.orf.block_starts = [0]
        self.orf.has_start_codon = True
        self.orf.has_stop_codon = True
        self.orf.transcriptomic = True
        self.assertFalse(self.orf.invalid)
        self.assertEqual(len(self.tr), self.tr.end - self.tr.start + 1)

    def test_print(self):

        self.maxDiff = None
        real_printed = """Chr1\tTAIR10\tmRNA\t5928\t8737\t.\t-\t.\tID=AT1G01020.1;Parent=AT1G01020;Name=AT1G01020.1
Chr1\tTAIR10\texon\t5928\t6263\t.\t-\t.\tID=AT1G01020.1.exon1;Parent=AT1G01020.1
Chr1\tTAIR10\tthree_prime_UTR\t5928\t6263\t.\t-\t.\tID=AT1G01020.1.three_prime_UTR1;Parent=AT1G01020.1
Chr1\tTAIR10\texon\t6437\t7069\t.\t-\t.\tID=AT1G01020.1.exon2;Parent=AT1G01020.1
Chr1\tTAIR10\tthree_prime_UTR\t6437\t6914\t.\t-\t.\tID=AT1G01020.1.three_prime_UTR2;Parent=AT1G01020.1
Chr1\tTAIR10\tCDS\t6915\t7069\t.\t-\t2\tID=AT1G01020.1.CDS1;Parent=AT1G01020.1
Chr1\tTAIR10\tCDS\t7157\t7232\t.\t-\t0\tID=AT1G01020.1.CDS2;Parent=AT1G01020.1
Chr1\tTAIR10\texon\t7157\t7232\t.\t-\t.\tID=AT1G01020.1.exon3;Parent=AT1G01020.1
Chr1\tTAIR10\tCDS\t7384\t7450\t.\t-\t1\tID=AT1G01020.1.CDS3;Parent=AT1G01020.1
Chr1\tTAIR10\texon\t7384\t7450\t.\t-\t.\tID=AT1G01020.1.exon4;Parent=AT1G01020.1
Chr1\tTAIR10\tCDS\t7564\t7649\t.\t-\t0\tID=AT1G01020.1.CDS4;Parent=AT1G01020.1
Chr1\tTAIR10\texon\t7564\t7649\t.\t-\t.\tID=AT1G01020.1.exon5;Parent=AT1G01020.1
Chr1\tTAIR10\tCDS\t7762\t7835\t.\t-\t2\tID=AT1G01020.1.CDS5;Parent=AT1G01020.1
Chr1\tTAIR10\texon\t7762\t7835\t.\t-\t.\tID=AT1G01020.1.exon6;Parent=AT1G01020.1
Chr1\tTAIR10\tCDS\t7942\t7987\t.\t-\t0\tID=AT1G01020.1.CDS6;Parent=AT1G01020.1
Chr1\tTAIR10\texon\t7942\t7987\t.\t-\t.\tID=AT1G01020.1.exon7;Parent=AT1G01020.1
Chr1\tTAIR10\tCDS\t8236\t8325\t.\t-\t0\tID=AT1G01020.1.CDS7;Parent=AT1G01020.1
Chr1\tTAIR10\texon\t8236\t8325\t.\t-\t.\tID=AT1G01020.1.exon8;Parent=AT1G01020.1
Chr1\tTAIR10\tCDS\t8417\t8464\t.\t-\t0\tID=AT1G01020.1.CDS8;Parent=AT1G01020.1
Chr1\tTAIR10\texon\t8417\t8464\t.\t-\t.\tID=AT1G01020.1.exon9;Parent=AT1G01020.1
Chr1\tTAIR10\tCDS\t8571\t8666\t.\t-\t0\tID=AT1G01020.1.CDS9;Parent=AT1G01020.1
Chr1\tTAIR10\texon\t8571\t8737\t.\t-\t.\tID=AT1G01020.1.exon10;Parent=AT1G01020.1
Chr1\tTAIR10\tfive_prime_UTR\t8667\t8737\t.\t-\t.\tID=AT1G01020.1.five_prime_UTR1;Parent=AT1G01020.1"""

        rp = set(real_printed.split("\n"))
        fp = set(str(self.tr).split("\n"))

        diff = "\n====\n".join(["\n".join(sorted(list(rp - set.intersection(rp, fp)))),
                               "\n".join(sorted(list(fp - set.intersection(rp, fp))))])

        self.assertEqual(real_printed,
                         str(self.tr),
                         diff)

        g_bed12 = "Chr1	5927	8737	ID=AT1G01020.1;coding=True;phase=0	0	-	6914	8666	0	10	336,633,76,67,86,74,46,90,48,167	0,509,1229,1456,1636,1834,2014,2308,2489,2643"
        self.assertEqual(g_bed12, self.tr.format("bed12", transcriptomic=False))

        t_bed12 = "AT1G01020.1	0	1623	ID=AT1G01020.1;coding=True;phase=0	0	+	71	809	0	10	167,48,90,46,74,86,67,76,633,336	0,167,215,305,351,425,511,578,654,1287"
        self.assertEqual(t_bed12, self.tr.format("bed12", transcriptomic=True))

    def test_empty(self):

        """
        Test that the inference of exons is valid.
        :return:
        """

        self.tr.finalized = False
        self.tr.exons = []
        self.tr.finalize()
        self.assertEqual(self.tr.strand, "-")
        self.assertEqual(self.tr.number_internal_orfs, 1)
        self.assertEqual(self.tr.exon_num, 10)
        self.assertEqual(self.tr.exon_num, len(self.tr.exons))
        self.assertEqual(self.tr.start, 5928)
        self.assertEqual(self.tr.end, 8737)
        exons = [(5928, 6263), (6437, 7069), (7157, 7232),
                 (7384, 7450), (7564, 7649), (7762, 7835),
                 (7942, 7987), (8236, 8325), (8417, 8464), (8571, 8737)]
        # exons = [intervaltree.Interval(*exon) for exon in exons]
        self.assertEqual(self.tr.exons,
                         exons,
                         self.tr.exons)
        # self.assertRaises(Mikado.exceptions.InvalidTranscript, self.tr.finalize)

    def test_invalid_utr(self):

        """
        Test that a transcript with UTR but no CDS defined will raise an exception.
        :return:
        """

        self.tr.combined_cds = []
        self.tr.finalized = False
        self.assertRaises(exceptions.InvalidTranscript, self.tr.finalize)

    def test_basics(self):

        """
        Test basic assertions about the transcript:

        - chromosome (.chrom) should be Chr1
        - strand should be -
        - number of internal orfs should be 1
        - number of exons should be 10
        - the metric "exon_num" should be 10 as well
        - start should be 5928 (1-based offset)
        - end should be 8737
        - the exons should correspond to those in the original strings (defined here in the list)
          and all of them should be of the "Interval" class

        :return:
        """

        self.assertEqual(self.tr.chrom, "Chr1")
        self.assertEqual(self.tr.strand, "-")
        self.assertEqual(self.tr.number_internal_orfs, 1)
        self.assertEqual(self.tr.exon_num, 10)
        self.assertEqual(self.tr.exon_num, len(self.tr.exons))
        self.assertEqual(self.tr.start, 5928)
        self.assertEqual(self.tr.end, 8737)
        exons = [(5928, 6263), (6437, 7069), (7157, 7232),
                 (7384, 7450), (7564, 7649), (7762, 7835),
                 (7942, 7987), (8236, 8325), (8417, 8464), (8571, 8737)]
        # exons = [intervaltree.Interval(*exon) for exon in exons]
        self.assertEqual(self.tr.exons,
                         exons,
                         self.tr.exons)

    def test_cds(self):
        self.assertEqual(sorted(self.tr.combined_cds),
                         sorted(self.tr.selected_cds))
        cds = [(6915, 7069), (7157, 7232), (7384, 7450), (7564, 7649), (7762, 7835), (7942, 7987),
               (8236, 8325), (8417, 8464), (8571, 8666)]

        self.assertEqual(self.tr.combined_cds,
                         cds,
                         self.tr.combined_cds)
        self.assertEqual(self.tr.selected_cds_start, 8666)
        self.assertEqual(self.tr.selected_cds_end, 6915)

    def test_utr(self):
        self.assertEqual(self.tr.five_utr, [(8667, 8737)])
        self.assertEqual(self.tr.three_utr, [(5928, 6263),
                                             (6437, 6914)])

    def test_utr_metrics(self):

        """Test for UTR exon num, start distance, etc."""

        self.assertEqual(self.tr.five_utr_num, 1)
        self.assertEqual(self.tr.five_utr_num_complete, 0)
        self.assertEqual(self.tr.three_utr_num, 2)
        self.assertEqual(self.tr.three_utr_num_complete, 1)

        self.assertEqual(self.tr.five_utr_length, 8737 + 1 - 8667)
        self.assertEqual(self.tr.three_utr_length, 6263 + 1 - 5928 + 6914 + 1 - 6437)
        self.assertEqual(self.tr.selected_start_distance_from_tss,
                         8738 - 8667,
                         self.tr.selected_end_distance_from_tes)
        self.assertEqual(self.tr.selected_end_distance_from_tes,
                         6263 + 1 - 5928 + 6915 - 6437,
                         self.tr.selected_end_distance_from_tes)
        self.assertEqual(self.tr.selected_end_distance_from_junction,
                         6915 - 6437,
                         self.tr.selected_cds_end)
        self.assertEqual(self.tr.end_distance_from_junction,
                         self.tr.selected_end_distance_from_junction)

    def test_introns(self):

        introns = {(8465, 8570), (8326, 8416), (7988, 8235),
                   (7836, 7941), (7650, 7761), (7451, 7563),
                   (7233, 7383), (7070, 7156), (6264, 6436)}

        self.assertEqual(self.tr.introns,
                         introns,
                         self.tr.introns)

        cds_introns = {(8465, 8570), (8326, 8416), (7988, 8235),
                       (7836, 7941), (7650, 7761), (7451, 7563),
                       (7233, 7383), (7070, 7156)}

        self.assertEqual(self.tr.combined_cds_introns,
                         cds_introns,
                         self.tr.combined_cds_introns)

        selected_cds_introns = {(8465, 8570), (8326, 8416), (7988, 8235),
                                (7836, 7941), (7650, 7761), (7451, 7563),
                                (7233, 7383), (7070, 7156)}

        self.assertEqual(self.tr.selected_cds_introns,
                         selected_cds_introns,
                         self.tr.selected_cds_introns)

    # @unittest.SkipTest
    def test_strip_cds(self):
        """
        Test the "stip_cds" function which (as the name implies) removes completely the CDS
        from a transcript.
        :return:
        """

        with self.assertLogs("null", level="DEBUG") as log_split:
            self.tr.strip_cds()
        self.assertIn("DEBUG:null:Stripping CDS from AT1G01020.1", log_split.output)

        self.assertEqual(self.tr.selected_cds_length, 0)
        self.assertEqual(self.tr.three_utr, [])
        self.assertEqual(self.tr.five_utr, [])
        self.assertEqual(self.tr.selected_cds, [])
        self.assertEqual(self.tr.selected_cds_start, None)
        self.assertEqual(self.tr.selected_cds_end, None)

    def test_remove_utr(self):
        """Test for CDS stripping. We remove the UTRs and verify that start/end
        have moved, no UTR is present, etc."""

        # tr = deepcopy(self.tr)
        self.tr.remove_utrs()

        # tr = deepcopy(self.tr)
        # tr.remove_utrs()
        self.assertEqual(self.tr.selected_cds_start, self.tr.end,
                         ((self.tr.selected_cds_start, self.tr.selected_cds_end),
                          (self.tr.start, self.tr.end)))
        self.assertEqual(self.tr.selected_cds_end, self.tr.start)
        self.assertEqual(self.tr.three_utr, [])
        self.assertEqual(self.tr.five_utr, [])
        combined_cds = [(6915, 7069),
                        (7157, 7232),
                        (7384, 7450),
                        (7564, 7649),
                        (7762, 7835),
                        (7942, 7987),
                        (8236, 8325),
                        (8417, 8464),
                        (8571, 8666)]
        self.assertEqual(self.tr.combined_cds,
                         combined_cds,
                         self.tr.combined_cds)
        self.assertEqual(self.tr.combined_utr, [], self.tr.combined_utr)

    def test_load_orf(self):

        """Test for loading a single ORF. We strip the CDS and reload it."""

        self.tr.strip_cds()
        self.tr.load_orfs([self.orf])

        combined_cds = [(6915, 7069), (7157, 7232), (7384, 7450),
                        (7564, 7649), (7762, 7835), (7942, 7987),
                        (8236, 8325), (8417, 8464), (8571, 8666)]

        self.assertEqual(self.tr.combined_cds,
                         combined_cds,
                         self.tr.combined_cds)
        self.assertEqual(self.tr.selected_cds_start, 8666)
        self.assertEqual(self.tr.selected_cds_end, 6915)

    def test_negative_orf(self):
        """Test loading a negative strand ORF onto a multiexonic transcript.
        This should have no effect.
        """

        self.orf.strand = "-"
        self.tr.strip_cds()
        self.tr.load_orfs([self.orf])
        self.assertEqual(self.tr.selected_cds_start, None)

    def testSegments(self):

        self.assertEqual(self.tr.combined_cds_num, 9)
        self.assertEqual(self.tr.selected_cds_num, 9)
        self.assertEqual(self.tr.highest_cds_exon_number, 9)
        self.assertEqual(self.tr.max_intron_length, 248)
        self.assertEqual(self.tr.number_internal_orfs, 1)

    def test_lengths(self):

        self.assertEqual(self.tr.cdna_length, 1623)
        self.assertEqual(self.tr.selected_cds_length, 738)
        self.assertAlmostEqual(self.tr.combined_cds_fraction, 738 / 1623, delta=0.01)
        self.assertAlmostEqual(self.tr.selected_cds_fraction, 738 / 1623, delta=0.01)

    def test_print_no_cds(self):

        self.maxDiff = None
        # tr = deepcopy(self.tr)
        # tr.finalize()

        real_printed = """Chr1\tTAIR10\tmRNA\t5928\t8737\t.\t-\t.\tID=AT1G01020.1;Parent=AT1G01020;Name=AT1G01020.1
Chr1\tTAIR10\texon\t5928\t6263\t.\t-\t.\tID=AT1G01020.1.exon1;Parent=AT1G01020.1
Chr1\tTAIR10\texon\t6437\t7069\t.\t-\t.\tID=AT1G01020.1.exon2;Parent=AT1G01020.1
Chr1\tTAIR10\texon\t7157\t7232\t.\t-\t.\tID=AT1G01020.1.exon3;Parent=AT1G01020.1
Chr1\tTAIR10\texon\t7384\t7450\t.\t-\t.\tID=AT1G01020.1.exon4;Parent=AT1G01020.1
Chr1\tTAIR10\texon\t7564\t7649\t.\t-\t.\tID=AT1G01020.1.exon5;Parent=AT1G01020.1
Chr1\tTAIR10\texon\t7762\t7835\t.\t-\t.\tID=AT1G01020.1.exon6;Parent=AT1G01020.1
Chr1\tTAIR10\texon\t7942\t7987\t.\t-\t.\tID=AT1G01020.1.exon7;Parent=AT1G01020.1
Chr1\tTAIR10\texon\t8236\t8325\t.\t-\t.\tID=AT1G01020.1.exon8;Parent=AT1G01020.1
Chr1\tTAIR10\texon\t8417\t8464\t.\t-\t.\tID=AT1G01020.1.exon9;Parent=AT1G01020.1
Chr1\tTAIR10\texon\t8571\t8737\t.\t-\t.\tID=AT1G01020.1.exon10;Parent=AT1G01020.1"""

        self.assertEqual(real_printed, self.tr.format("gff", with_cds=False))

        real_printed_gtf = """Chr1\tTAIR10\tmRNA\t5928\t8737\t.\t-\t.\tgene_id "AT1G01020"; transcript_id "AT1G01020.1"; Name "AT1G01020.1";
Chr1\tTAIR10\texon\t5928\t6263\t.\t-\t.\tgene_id "AT1G01020"; transcript_id "AT1G01020.1";
Chr1\tTAIR10\texon\t6437\t7069\t.\t-\t.\tgene_id "AT1G01020"; transcript_id "AT1G01020.1";
Chr1\tTAIR10\texon\t7157\t7232\t.\t-\t.\tgene_id "AT1G01020"; transcript_id "AT1G01020.1";
Chr1\tTAIR10\texon\t7384\t7450\t.\t-\t.\tgene_id "AT1G01020"; transcript_id "AT1G01020.1";
Chr1\tTAIR10\texon\t7564\t7649\t.\t-\t.\tgene_id "AT1G01020"; transcript_id "AT1G01020.1";
Chr1\tTAIR10\texon\t7762\t7835\t.\t-\t.\tgene_id "AT1G01020"; transcript_id "AT1G01020.1";
Chr1\tTAIR10\texon\t7942\t7987\t.\t-\t.\tgene_id "AT1G01020"; transcript_id "AT1G01020.1";
Chr1\tTAIR10\texon\t8236\t8325\t.\t-\t.\tgene_id "AT1G01020"; transcript_id "AT1G01020.1";
Chr1\tTAIR10\texon\t8417\t8464\t.\t-\t.\tgene_id "AT1G01020"; transcript_id "AT1G01020.1";
Chr1\tTAIR10\texon\t8571\t8737\t.\t-\t.\tgene_id "AT1G01020"; transcript_id "AT1G01020.1";"""

        self.assertEqual(real_printed_gtf, self.tr.__str__(print_cds=False, to_gtf=True))

    def testDoubleOrf(self):

        """Test to verify the introduction of multiple ORFs."""

        self.tr.strip_cds()
        self.tr.finalized = False

        first_orf = parsers.bed12.BED12()
        first_orf.chrom = self.tr.id
        first_orf.start = 1
        first_orf.end = self.tr.cdna_length
        first_orf.name = self.tr.id
        first_orf.strand = "+"
        first_orf.score = 0
        first_orf.thick_start = 100
        first_orf.thick_end = 501
        first_orf.block_count = 1
        first_orf.blockSize = self.tr.cdna_length
        first_orf.block_starts = [0]
        first_orf.has_start_codon = True
        first_orf.has_stop_codon = True
        first_orf.transcriptomic = True
        self.assertFalse(first_orf.invalid, (len(first_orf), first_orf.cds_len))

        # This should not be incorporated
        second_orf = parsers.bed12.BED12()
        second_orf.chrom = self.tr.id
        second_orf.start = 0
        second_orf.end = self.tr.cdna_length
        second_orf.name = "second"
        second_orf.strand = "+"
        second_orf.score = 1
        second_orf.thick_start = 300
        second_orf.thick_end = 401
        second_orf.block_count = 1
        second_orf.blockSize = self.tr.cdna_length
        second_orf.block_starts = [0]
        second_orf.has_start_codon = True
        second_orf.has_stop_codon = True
        second_orf.transcriptomic = True
        self.assertFalse(second_orf.invalid, (len(second_orf), second_orf.cds_len))

        self.assertTrue(loci.Transcript.is_overlapping_cds(first_orf, second_orf))

        # This should be added
        third_orf = parsers.bed12.BED12()
        third_orf.chrom = self.tr.id
        third_orf.start = 1
        third_orf.end = self.tr.cdna_length
        third_orf.name = "third"
        third_orf.strand = "+"
        third_orf.score = 0
        third_orf.thick_start = 1000
        third_orf.thick_end = 1602
        third_orf.block_count = 1
        third_orf.blockSize = self.tr.cdna_length
        third_orf.block_starts = [0]
        third_orf.has_start_codon = True
        third_orf.has_stop_codon = True
        third_orf.transcriptomic = True
        self.assertFalse(third_orf.invalid, (len(third_orf), third_orf.cds_len))

        self.assertFalse(
            loci.Transcript.is_overlapping_cds(
                first_orf, third_orf))
        self.assertFalse(
            loci.Transcript.is_overlapping_cds(
                second_orf, third_orf))

        self.assertFalse(third_orf == second_orf)
        self.assertFalse(first_orf == second_orf)
        self.assertFalse(first_orf == third_orf)

        candidates = [first_orf, second_orf, third_orf]
        self.tr.load_orfs(candidates)

        self.assertTrue(self.tr.is_complete)
        self.tr.finalize()
        self.assertEqual(self.tr.start, 5928)

        self.assertEqual(self.tr.number_internal_orfs, 2, "\n".join(
            [str(x) for x in self.tr.internal_orfs]))

        self.assertEqual(self.tr.combined_cds_length, 1005)
        self.assertEqual(self.tr.selected_cds_length, 603)

        new_transcripts = sorted(self.tr.split_by_cds(), key=operator.attrgetter("start"))

        self.assertEqual(len(new_transcripts), 2)
        self.assertEqual(new_transcripts[0].five_utr_length, 0)
        self.assertNotEqual(new_transcripts[0].three_utr_length, 0)
        self.assertEqual(new_transcripts[0].cdna_length, 624, msg="{0}-{1}{2}".format(
            new_transcripts[0].start,
            new_transcripts[0].end,
            new_transcripts[0].strand,
        ))
        self.assertEqual(new_transcripts[0].start, self.tr.start)
        self.assertEqual(new_transcripts[0].end, 6724)

        self.assertEqual(new_transcripts[1].three_utr_length, 0)
        self.assertEqual(new_transcripts[1].end, 8737)

class WheatTest(unittest.TestCase):

    def setUp(self):

        trlines= """Triticum_aestivum_CS42_TGACv1_scaffold_018984_1AS\tta_complete\ttranscript\t217224\t221661\t100\t-\t.\tgene_id "TGAC_Trinity_Triticum_aestivum_CS42_TGACv1_scaffold_018984_1AS:216942_226036:c0_g1_i4.mrna1"; transcript_id "TGAC_Trinity_Triticum_aestivum_CS42_TGACv1_scaffold_018984_1AS:216942_226036:c0_g1_i4.mrna1"; Name "TGAC_Trinity_Triticum_aestivum_CS42_TGACv1_scaffold_018984_1AS:216942_226036:c0_g1_i4.mrna1"; canonical_junctions "1,2,3,4,5"; canonical_proportion "1.0"; canonical_number "5"; target "Triticum_aestivum_CS42_TGACv1_scaffold_018984_1AS:216942_226036:c0_g1_i4 1 1251 +";
Triticum_aestivum_CS42_TGACv1_scaffold_018984_1AS\tta_complete\texon\t217224\t218053\t.\t-\t.\tgene_id "TGAC_Trinity_Triticum_aestivum_CS42_TGACv1_scaffold_018984_1AS:216942_226036:c0_g1_i4.mrna1"; transcript_id "TGAC_Trinity_Triticum_aestivum_CS42_TGACv1_scaffold_018984_1AS:216942_226036:c0_g1_i4.mrna1";
Triticum_aestivum_CS42_TGACv1_scaffold_018984_1AS\tta_complete\texon\t218145\t218384\t.\t-\t.\tgene_id "TGAC_Trinity_Triticum_aestivum_CS42_TGACv1_scaffold_018984_1AS:216942_226036:c0_g1_i4.mrna1"; transcript_id "TGAC_Trinity_Triticum_aestivum_CS42_TGACv1_scaffold_018984_1AS:216942_226036:c0_g1_i4.mrna1";
Triticum_aestivum_CS42_TGACv1_scaffold_018984_1AS\tta_complete\texon\t218474\t219472\t.\t-\t.\tgene_id "TGAC_Trinity_Triticum_aestivum_CS42_TGACv1_scaffold_018984_1AS:216942_226036:c0_g1_i4.mrna1"; transcript_id "TGAC_Trinity_Triticum_aestivum_CS42_TGACv1_scaffold_018984_1AS:216942_226036:c0_g1_i4.mrna1";
Triticum_aestivum_CS42_TGACv1_scaffold_018984_1AS\tta_complete\texon\t219557\t220122\t.\t-\t.\tgene_id "TGAC_Trinity_Triticum_aestivum_CS42_TGACv1_scaffold_018984_1AS:216942_226036:c0_g1_i4.mrna1"; transcript_id "TGAC_Trinity_Triticum_aestivum_CS42_TGACv1_scaffold_018984_1AS:216942_226036:c0_g1_i4.mrna1";
Triticum_aestivum_CS42_TGACv1_scaffold_018984_1AS\tta_complete\texon\t220205\t220327\t.\t-\t.\tgene_id "TGAC_Trinity_Triticum_aestivum_CS42_TGACv1_scaffold_018984_1AS:216942_226036:c0_g1_i4.mrna1"; transcript_id "TGAC_Trinity_Triticum_aestivum_CS42_TGACv1_scaffold_018984_1AS:216942_226036:c0_g1_i4.mrna1";
Triticum_aestivum_CS42_TGACv1_scaffold_018984_1AS\tta_complete\texon\t220411\t221661\t.\t-\t.\tgene_id "TGAC_Trinity_Triticum_aestivum_CS42_TGACv1_scaffold_018984_1AS:216942_226036:c0_g1_i4.mrna1"; transcript_id "TGAC_Trinity_Triticum_aestivum_CS42_TGACv1_scaffold_018984_1AS:216942_226036:c0_g1_i4.mrna1";"""

        trlines = [parsers.GTF.GtfLine(_) for _ in trlines.split("\n")]

        self.tr = loci.Transcript(trlines[0])
        for _ in trlines[1:]:
            self.tr.add_exon(_)

        self.bed1 = parsers.bed12.BED12()
        self.bed1.header = False
        self.bed1.chrom = self.tr.id
        self.bed1.transcriptomic = True
        self.bed1.start = 1
        self.bed1.end = 4009
        self.bed1.name = "g.4208_TGAC_Trinity_Triticum_aestivum_CS42_TGACv1_scaffold_018984_1AS:216942_226036:c0_g1_i4.mrna1|m.4208_type:complete_len:193_(+)"
        self.bed1.strand = "+"
        self.bed1.thick_start = 297
        self.bed1.thick_end = 875
        self.bed1.has_start_codon = True
        self.bed1.has_stop_codon = True
        self.bed1.score = 0
        self.assertEqual(self.bed1.cds_len, 579)

        self.bed2 = parsers.bed12.BED12()
        self.bed2.header = False
        self.bed2.chrom = self.tr.id
        self.bed2.transcriptomic = True
        self.bed2.start = 1
        self.bed2.end = 4009
        self.bed2.name = "TGAC_Trinity_Triticum_aestivum_CS42_TGACv1_scaffold_018984_1AS:216942_226036:c0_g1_i4.mrna1|m.4209_type:complete_len:116_(+)"
        self.bed2.strand = "+"
        self.bed2.thick_start = 2886
        self.bed2.thick_end = 3233
        self.bed2.has_start_codon = True
        self.bed2.has_stop_codon = True
        self.bed2.score = 0
        self.assertEqual(self.bed2.cds_len, 348)

    def test_load(self):
        
        self.tr.load_orfs([self.bed1, self.bed2])
        self.assertEqual(self.tr.number_internal_orfs, 2)
        self.assertEqual(self.tr.selected_cds_start, 221365)
        self.assertEqual(self.tr.selected_cds_end, 220787)

        self.assertEqual(self.tr.combined_cds_start, 221365)
        self.assertEqual(self.tr.combined_cds_end, 218000)

    def test_split(self):

        self.tr.load_orfs([self.bed1, self.bed2])
        new_transcripts = sorted([_ for _ in self.tr.split_by_cds()],
                                 key=operator.attrgetter("start", "end"))
        self.assertEqual(len(new_transcripts), 2)

        self.assertEqual(new_transcripts[0].start, 217224)
        self.assertEqual(new_transcripts[0].end, 218527)

        self.assertEqual(new_transcripts[1].start, 220787)
        self.assertEqual(new_transcripts[1].end, 221661)

        self.assertEqual(new_transcripts[0].selected_cds_start, 218527)
        self.assertEqual(new_transcripts[0].selected_cds_end, 218000)

        self.assertEqual(new_transcripts[1].selected_cds_start, 221365)
        self.assertEqual(new_transcripts[1].selected_cds_end, 220787)


if __name__ == '__main__':
    unittest.main()
