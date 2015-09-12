# coding: utf-8

"""
Unit test for a transcript on the negative strand.
"""

import unittest
import logging
import operator
import re
import mikado_lib.parsers
import mikado_lib.loci_objects


class TranscriptTesterNegative(unittest.TestCase):

    handler = logging.StreamHandler()
    handler.setLevel("DEBUG")
    logger = logging.getLogger("test")
    logger.setLevel("DEBUG")
    logger.propagate = False
    formatter = logging.Formatter(
        "{asctime} - {levelname} - {module}:{lineno} - {funcName} - {name} - {message}",
        style="{"
        )
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    tr_gff = """Chr1    TAIR10    mRNA    5928    8737    .    -    .    ID=AT1G01020.1;Parent=AT1G01020;Name=AT1G01020.1;Index=1
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

    tr_lines = tr_gff.split("\n")
    for pos, line in enumerate(tr_lines):
        tr_lines[pos] = re.sub("\s+", "\t", line)
        assert len(tr_lines[pos].split("\t")) == 9, line.split("\t")

    tr_gff_lines = [mikado_lib.parsers.GFF.GffLine(line) for line in tr_lines]

    for l in tr_gff_lines:
        assert l.header is False
    #         print(l)

    def setUp(self):
        """Basic creation test."""

        self.tr = mikado_lib.loci_objects.transcript.Transcript(self.tr_gff_lines[0])
        for line in self.tr_gff_lines[1:]:
            self.tr.add_exon(line)
        self.tr.finalize()
        self.tr.logger = self.logger

        self.orf = mikado_lib.parsers.bed12.BED12()
        self.orf.chrom = self.tr.id
        self.orf.start = 1
        self.orf.end = self.tr.cdna_length
        self.orf.name = self.tr.id
        self.orf.strand = "+"
        self.orf.score = 0
        self.orf.thickStart = self.tr.selected_start_distance_from_tss + 1
        self.orf.thickEnd = self.tr.cdna_length - self.tr.selected_end_distance_from_tes
        self.orf.blockCount = 1
        self.orf.blockSize = self.tr.cdna_length
        self.orf.blockStarts = 0
        self.orf.has_start_codon = True
        self.orf.has_stop_codon = True
        self.orf.transcriptomic = False
        self.assertFalse(self.orf.invalid)

    def test_basics(self):

        self.assertEqual(self.tr.chrom, "Chr1")
        self.assertEqual(self.tr.strand, "-")
        self.assertEqual(self.tr.exon_num, 10)
        self.assertEqual(self.tr.exon_num, len(self.tr.exons))
        self.assertEqual(self.tr.start, 5928)
        self.assertEqual(self.tr.end, 8737)
        self.assertEqual(self.tr.exons,
                         [(5928, 6263), (6437, 7069), (7157, 7232), (7384, 7450), (7564, 7649), (7762, 7835),
                          (7942, 7987), (8236, 8325), (8417, 8464), (8571, 8737)],
                         self.tr.exons)

    def test_cds(self):
        self.assertEqual(self.tr.combined_cds, self.tr.selected_cds)

        self.assertEqual(self.tr.combined_cds,
                         [(6915, 7069), (7157, 7232), (7384, 7450), (7564, 7649), (7762, 7835), (7942, 7987),
                          (8236, 8325), (8417, 8464), (8571, 8666)],
                         self.tr.combined_cds)
        self.assertEqual(self.tr.selected_cds_start, 8666)
        self.assertEqual(self.tr.selected_cds_end, 6915)

    def test_utr(self):
        self.assertEqual(self.tr.five_utr, [("UTR", 8667, 8737)])
        self.assertEqual(self.tr.three_utr, [("UTR", 5928, 6263), ("UTR", 6437, 6914)])

    def test_utr_metrics(self):

        """Test for UTR exon num, start distance, etc."""

        self.assertEqual(self.tr.five_utr_num, 1)
        self.assertEqual(self.tr.five_utr_num_complete, 0)
        self.assertEqual(self.tr.three_utr_num, 2)
        self.assertEqual(self.tr.three_utr_num_complete, 1)

        self.assertEqual(self.tr.five_utr_length, 8737 + 1 - 8667)
        self.assertEqual(self.tr.three_utr_length, 6263 + 1 - 5928 + 6914 + 1 - 6437)
        self.assertEqual(self.tr.selected_start_distance_from_tss, 8738 - 8667, self.tr.selected_end_distance_from_tes)
        self.assertEqual(self.tr.selected_end_distance_from_tes, 6263 + 1 - 5928 + 6915 - 6437,
                         self.tr.selected_end_distance_from_tes)
        self.assertEqual(self.tr.selected_end_distance_from_junction, 6915 - 6437 + 1)
        self.assertEqual(self.tr.end_distance_from_junction, self.tr.selected_end_distance_from_junction)

    def test_introns(self):

        self.assertEqual(self.tr.introns,
                         {(8465, 8570), (8326, 8416), (7988, 8235), (7836, 7941), (7650, 7761), (7451, 7563),
                          (7233, 7383), (7070, 7156), (6264, 6436)},
                         self.tr.introns
                         )
        self.assertEqual(self.tr.combined_cds_introns,
                         {(8465, 8570), (8326, 8416), (7988, 8235), (7836, 7941), (7650, 7761), (7451, 7563),
                          (7233, 7383), (7070, 7156)},
                         self.tr.combined_cds_introns
                         )
        self.assertEqual(self.tr.selected_cds_introns,
                         {(8465, 8570), (8326, 8416), (7988, 8235), (7836, 7941), (7650, 7761), (7451, 7563),
                          (7233, 7383), (7070, 7156)},
                         self.tr.selected_cds_introns
                         )

    def test_strip_cds(self):

        self.tr.strip_cds()
        self.assertEqual(self.tr.selected_cds_length, 0)
        self.assertEqual(self.tr.three_utr, [])
        self.assertEqual(self.tr.five_utr, [])
        self.assertEqual(self.tr.selected_cds, [])
        self.assertEqual(self.tr.selected_cds_start, None)
        self.assertEqual(self.tr.selected_cds_end, None)

    def test_remove_utr(self):
        """Test for CDS stripping. We remove the UTRs and verify that start/end have moved, no UTR is present, etc."""

        self.tr.remove_utrs()
        self.assertEqual(self.tr.selected_cds_start, self.tr.end)
        self.assertEqual(self.tr.selected_cds_end, self.tr.start)
        self.assertEqual(self.tr.three_utr, [])
        self.assertEqual(self.tr.five_utr, [])
        self.assertEqual(self.tr.combined_cds,
                         [(6915, 7069),
                          (7157, 7232),
                          (7384, 7450),
                          (7564, 7649),
                          (7762, 7835),
                          (7942, 7987),
                          (8236, 8325),
                          (8417, 8464),
                          (8571, 8666)],
                         self.tr.combined_cds)
        self.assertEqual(self.tr.combined_utr, [], self.tr.combined_utr)

    def test_load_orf(self):

        """Test for loading a single ORF. We strip the CDS and reload it."""

        self.tr.strip_cds()
        self.tr.load_orfs([self.orf])
        self.assertEqual(self.tr.combined_cds,
                         [(6915, 7069), (7157, 7232), (7384, 7450), (7564, 7649), (7762, 7835), (7942, 7987),
                          (8236, 8325), (8417, 8464), (8571, 8666)],
                         self.tr.combined_cds)
        self.assertEqual(self.tr.selected_cds_start, 8666)
        self.assertEqual(self.tr.selected_cds_end, 6915)

    def test_negative_orf(self):
        """Test loading a negative strand ORF onto a multiexonic transcript. This should have no effect."""

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

    def testDoubleOrf(self):

        """Test to verify the introduction of multiple ORFs."""

        self.tr.strip_cds()
        self.tr.finalized = False

        first_orf = mikado_lib.parsers.bed12.BED12()
        first_orf.chrom = self.tr.id
        first_orf.start = 1
        first_orf.end = self.tr.cdna_length
        first_orf.name = self.tr.id
        first_orf.strand = "+"
        first_orf.score = 0
        first_orf.thickStart = 100
        first_orf.thickEnd = 501
        first_orf.blockCount = 1
        first_orf.blockSize = self.tr.cdna_length
        first_orf.blockStarts = 0
        first_orf.has_start_codon = True
        first_orf.has_stop_codon = True
        first_orf.transcriptomic = True
        self.assertFalse(first_orf.invalid, (len(first_orf), first_orf.cds_len))

        # This should not be incorporated
        second_orf = mikado_lib.parsers.bed12.BED12()
        second_orf.chrom = self.tr.id
        second_orf.start = 0
        second_orf.end = self.tr.cdna_length
        second_orf.name = "second"
        second_orf.strand = "+"
        second_orf.score = 1
        second_orf.thickStart = 300
        second_orf.thickEnd = 401
        second_orf.blockCount = 1
        second_orf.blockSize = self.tr.cdna_length
        second_orf.blockStarts = 0
        second_orf.has_start_codon = True
        second_orf.has_stop_codon = True
        second_orf.transcriptomic = True
        self.assertFalse(second_orf.invalid, (len(second_orf), second_orf.cds_len))

        self.assertTrue(mikado_lib.loci_objects.transcript.Transcript.is_overlapping_cds(first_orf, second_orf))

        # This should be added
        third_orf = mikado_lib.parsers.bed12.BED12()
        third_orf.chrom = self.tr.id
        third_orf.start = 1
        third_orf.end = self.tr.cdna_length
        third_orf.name = "third"
        third_orf.strand = "+"
        third_orf.score = 0
        third_orf.thickStart = 1000
        third_orf.thickEnd = 1602
        third_orf.blockCount = 1
        third_orf.blockSize = self.tr.cdna_length
        third_orf.blockStarts = 0
        third_orf.has_start_codon = True
        third_orf.has_stop_codon = True
        third_orf.transcriptomic = True
        self.assertFalse(third_orf.invalid, (len(third_orf), third_orf.cds_len))

        self.assertFalse(mikado_lib.loci_objects.transcript.Transcript.is_overlapping_cds(first_orf, third_orf))
        self.assertFalse(mikado_lib.loci_objects.transcript.Transcript.is_overlapping_cds(second_orf, third_orf))

        self.assertFalse(third_orf == second_orf)
        self.assertFalse(first_orf == second_orf)
        self.assertFalse(first_orf == third_orf)

        candidates = [first_orf, second_orf, third_orf]
        # self.assertEqual(len(mikado_lib.loci_objects.transcript.Transcript.find_overlapping_cds(candidates)), 2)
        self.tr.load_orfs(candidates)

        self.assertTrue(self.tr.is_complete)
        self.tr.finalize()
        self.assertEqual(self.tr.start, 5928)

        self.assertEqual(self.tr.number_internal_orfs, 2, "\n".join([str(x) for x in self.tr.internal_orfs]))

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


    # def testDoubleOrfSameExon(self):
    #     """
    #     Test for those cases where the two ORFs share an exon.
    #     :return:
    #     """
    #
    #     self.tr.strip_cds()
    #     self.tr.finalized = False
    #
    #     first_orf = mikado_lib.parsers.bed12.BED12()
    #     first_orf.chrom = self.tr.id
    #     first_orf.start = 1
    #     first_orf.end = self.tr.cdna_length
    #     first_orf.name = self.tr.id
    #     first_orf.strand = "+"
    #     first_orf.score = 0
    #     first_orf.thickStart = 100
    #     first_orf.thickEnd = 501
    #     first_orf.blockCount = 1
    #     first_orf.blockSize = self.tr.cdna_length
    #     first_orf.blockStarts = 0
    #     first_orf.has_start_codon = True
    #     first_orf.has_stop_codon = True
    #     first_orf.transcriptomic = True
    #     self.assertFalse(first_orf.invalid, (len(first_orf), first_orf.cds_len))
    #
    #     # This should not be incorporated
    #     second_orf = mikado_lib.parsers.bed12.BED12()
    #     second_orf.chrom = self.tr.id
    #     second_orf.start = 0
    #     second_orf.end = self.tr.cdna_length
    #     second_orf.name = "second"
    #     second_orf.strand = "+"
    #     second_orf.score = 1
    #     second_orf.thickStart = 300
    #     second_orf.thickEnd = 401
    #     second_orf.blockCount = 1
    #     second_orf.blockSize = self.tr.cdna_length
    #     second_orf.blockStarts = 0
    #     second_orf.has_start_codon = True
    #     second_orf.has_stop_codon = True
    #     second_orf.transcriptomic = True
    #     self.assertFalse(second_orf.invalid, (len(second_orf), second_orf.cds_len))
    #
    #     self.assertTrue(mikado_lib.loci_objects.transcript.Transcript.is_overlapping_cds(first_orf, second_orf))
    #
    #     # This should be added
    #     third_orf = mikado_lib.parsers.bed12.BED12()
    #     third_orf.chrom = self.tr.id
    #     third_orf.start = 1
    #     third_orf.end = self.tr.cdna_length
    #     third_orf.name = "third"
    #     third_orf.strand = "+"
    #     third_orf.score = 0
    #     third_orf.thickStart = 1000
    #     third_orf.thickEnd = 1602
    #     third_orf.blockCount = 1
    #     third_orf.blockSize = self.tr.cdna_length
    #     third_orf.blockStarts = 0
    #     third_orf.has_start_codon = True
    #     third_orf.has_stop_codon = True
    #     third_orf.transcriptomic = True
    #     self.assertFalse(third_orf.invalid, (len(third_orf), third_orf.cds_len))
    #
    #     self.assertFalse(mikado_lib.loci_objects.transcript.Transcript.is_overlapping_cds(first_orf, third_orf))
    #     self.assertFalse(mikado_lib.loci_objects.transcript.Transcript.is_overlapping_cds(second_orf, third_orf))
    #
    #     self.assertFalse(third_orf == second_orf)
    #     self.assertFalse(first_orf == second_orf)
    #     self.assertFalse(first_orf == third_orf)
    #
    #     candidates = [first_orf, second_orf, third_orf]
    #     self.assertEqual(len(mikado_lib.loci_objects.transcript.Transcript.find_overlapping_cds(candidates)), 2)
    #     self.tr.load_orfs(candidates)
    #
    #     self.assertTrue(self.tr.is_complete)
    #     self.tr.finalize()
    #
    #     self.assertEqual(self.tr.number_internal_orfs, 2, "\n".join([str(x) for x in self.tr.internal_orfs]))
    #
    #     self.assertEqual(self.tr.combined_cds_length, 1005)
    #     self.assertEqual(self.tr.selected_cds_length, 603)
    #
    #     new_transcripts = sorted(self.tr.split_by_cds())
    #
    #     self.assertEqual(len(new_transcripts), 2)
    #     self.assertEqual(new_transcripts[1].five_utr_length, 0)
    #     self.assertEqual(new_transcripts[0].three_utr_length, 0)


unittest.main()
