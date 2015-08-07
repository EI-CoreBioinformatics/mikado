# coding: utf-8

"""
Unit test for monoexonic transcripts.
"""

import unittest
import re
import mikado_lib.parsers
import mikado_lib.loci_objects


class TranscriptTester(unittest.TestCase):
    tr_gff = """Chr1    TAIR10    mRNA    5928    8737    .    .    .    ID=AT1G01020.1;Parent=AT1G01020;Name=AT1G01020.1;Index=1
Chr1    TAIR10    three_prime_UTR    8667    8737    .    .    .    Parent=AT1G01020.1
Chr1    TAIR10    CDS    8571    8666    .    .    0    Parent=AT1G01020.1;
Chr1    TAIR10    exon    5928    8737    .    .    .    Parent=AT1G01020.1
Chr1    TAIR10    five_prime_UTR    5928    8570    .    -    .    Parent=AT1G01020.1"""

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
        self.orf.transcriptomic = True
        self.assertFalse(self.orf.invalid)

    def test_basics(self):

        self.assertEqual(self.tr.chrom, "Chr1")
        self.assertEqual(self.tr.strand, None)
        self.assertEqual(self.tr.exon_num, 1)
        self.assertEqual(self.tr.monoexonic, True)
        self.assertEqual(self.tr.exon_num, len(self.tr.exons))
        self.assertEqual(self.tr.start, 5928)
        self.assertEqual(self.tr.end, 8737)
        self.assertEqual(self.tr.exons,
                         [(5928, 8737)],
                         self.tr.exons)

    def test_cds(self):
        """Test the CDS features.
        Note that in a single-exon transcript with no strand, start_codon and stop_codon are defined as False.
        """

        self.assertEqual(self.tr.combined_cds, self.tr.selected_cds)

        self.assertEqual(self.tr.combined_cds,
                         [(8571, 8666)],
                         self.tr.combined_cds)
        self.assertEqual(self.tr.selected_cds_start, 8571)
        self.assertEqual(self.tr.selected_cds_end, 8666)
        self.assertEqual(self.tr.has_start_codon, False)
        self.assertEqual(self.tr.has_stop_codon, False)

    def test_utr(self):

        self.assertEqual(self.tr.selected_internal_orf,
                         [("UTR", 5928, 8570), ("exon", 5928, 8737), ("CDS", 8571, 8666), ("UTR", 8667, 8737)],
                         "Right: {0}\nFound{1}".format([("UTR", 5928, 8570), ("CDS", 8571, 8666), ("UTR", 8667, 8737)],
                                                       self.tr.selected_internal_orf))
        self.assertEqual(self.tr.combined_utr, [(5928, 8570), (8667, 8737)])
        self.assertEqual(self.tr.five_utr, [("UTR", 5928, 8570)], self.tr.five_utr)
        self.assertEqual(self.tr.three_utr, [("UTR", 8667, 8737)])

    def test_utr_metrics(self):

        """Test for UTR exon num, start distance, etc."""

        self.assertEqual(self.tr.five_utr_num, 1)
        self.assertEqual(self.tr.three_utr_num, 1)
        self.assertEqual(self.tr.five_utr_length, 8570 + 1 - 5928)
        self.assertEqual(self.tr.three_utr_length, 8737 + 1 - 8667)
        self.assertEqual(self.tr.selected_start_distance_from_tss, 8571 - 5928, self.tr.selected_end_distance_from_tes)
        self.assertEqual(self.tr.selected_end_distance_from_tes, 8737 - 8666, self.tr.selected_end_distance_from_tes)

    def test_strip_cds(self):

        self.tr.strip_cds()
        self.assertEqual(self.tr.selected_cds_length, 0)
        self.assertEqual(self.tr.three_utr, [])
        self.assertEqual(self.tr.five_utr, [])
        self.assertEqual(self.tr.selected_cds, [])
        self.assertEqual(self.tr.selected_cds_start, None)
        self.assertEqual(self.tr.selected_cds_end, None)

    def test_remove_utr(self):
        """Test for CDS stripping. We remove the UTRs and verify that start/end have moved, no UTR is present, etc.
        """

        self.tr.remove_utrs()
        self.assertEqual(self.tr.selected_cds_start, self.tr.start)
        self.assertEqual(self.tr.selected_cds_end, self.tr.end)
        self.assertEqual(self.tr.three_utr, [])
        self.assertEqual(self.tr.five_utr, [])
        self.assertEqual(self.tr.combined_cds,
                         [(8571, 8666)],
                         self.tr.combined_cds)
        self.assertEqual(self.tr.combined_utr, [], self.tr.combined_utr)

    def test_load_orf(self):

        """Test for loading a single ORF. We strip the CDS and reload it."""

        self.tr.strip_cds()
        self.tr.load_orfs([self.orf])
        self.assertEqual(self.tr.combined_cds,
                         [(8571, 8666)],
                         self.tr.combined_cds)
        self.assertEqual(self.tr.selected_cds_start, 8571)
        self.assertEqual(self.tr.selected_cds_end, 8666)
        self.assertEqual(self.tr.has_start_codon, True)
        self.assertEqual(self.tr.has_stop_codon, True)

    def test_negative_orf(self):
        """Test loading a negative strand ORF onto a monoexonic transcript. This should reverse the ORF."""

        self.orf.strand = "-"
        self.tr.strip_cds()
        self.tr.load_orfs([self.orf])
        self.assertEqual(self.tr.strand, "-")
        self.assertEqual(self.tr.selected_cds_start, 8737 - (8571 - 5928))
        self.assertEqual(self.tr.selected_cds_end, 5928 + (8737 - 8666))

    def test_introns(self):

        self.assertEqual(self.tr.introns,
                         set([
                         ]),
                         self.tr.introns
                         )
        self.assertEqual(self.tr.combined_cds_introns,
                         set([
                         ]),
                         self.tr.combined_cds_introns
                         )
        self.assertEqual(self.tr.selected_cds_introns,
                         set([
                         ]),
                         self.tr.selected_cds_introns
                         )

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
        first_orf.thickStart = 51
        first_orf.thickEnd = 398
        first_orf.blockCount = 1
        first_orf.blockSize = self.tr.cdna_length
        first_orf.blockSizes = [self.tr.cdna_length]
        first_orf.blockStarts = [0]
        first_orf.rgb = 0
        first_orf.has_start_codon = True
        first_orf.has_stop_codon = True
        first_orf.transcriptomic = True
        self.assertFalse(first_orf.invalid)
        # This should not be incorporated
        second_orf = mikado_lib.parsers.bed12.BED12()
        second_orf.chrom = self.tr.id
        second_orf.start = 1
        second_orf.end = self.tr.cdna_length
        second_orf.name = "second"
        second_orf.strand = "+"
        second_orf.score = 0
        second_orf.thickStart = 201
        second_orf.thickEnd = 410
        second_orf.blockCount = 1
        second_orf.blockSize = self.tr.cdna_length
        second_orf.blockSizes = [self.tr.cdna_length]
        second_orf.blockStarts = [0]
        second_orf.rgb = 0
        second_orf.has_start_codon = True
        second_orf.has_stop_codon = True
        second_orf.transcriptomic = True
        self.assertFalse(second_orf.invalid)

        self.assertTrue(mikado_lib.loci_objects.transcript.Transcript.is_overlapping_cds(first_orf, second_orf))

        # This should be added
        third_orf = mikado_lib.parsers.bed12.BED12()
        third_orf.chrom = self.tr.id
        third_orf.start = 1
        third_orf.end = self.tr.cdna_length
        third_orf.name = "third"
        third_orf.strand = "+"
        third_orf.score = 0
        third_orf.thickStart = 501
        third_orf.thickEnd = 800
        third_orf.blockCount = 1
        third_orf.blockSize = self.tr.cdna_length
        third_orf.blockSizes = [self.tr.cdna_length]
        third_orf.blockStarts = [0]
        third_orf.rgb = 0
        third_orf.has_start_codon = True
        third_orf.has_stop_codon = True
        third_orf.transcriptomic = True
        self.assertFalse(third_orf.invalid)

        self.assertFalse(mikado_lib.loci_objects.transcript.Transcript.is_overlapping_cds(first_orf, third_orf))
        self.assertFalse(mikado_lib.loci_objects.transcript.Transcript.is_overlapping_cds(second_orf, third_orf))

        self.assertFalse(third_orf == second_orf)
        self.assertFalse(first_orf == second_orf)
        self.assertFalse(first_orf == third_orf)

        candidates = [first_orf, second_orf, third_orf]

        # self.assertEqual(len(mikado_lib.loci_objects.transcript.Transcript.find_overlapping_cds(candidates)), 2)

        import logging
        handler = logging.StreamHandler()
        handler.setLevel("DEBUG")
        logger = logging.getLogger("test")
        logger.setLevel("DEBUG")
        logger.addHandler(handler)
        logger.debug("TEST")
        logger.propagate = False
        self.tr.logger = logger

        self.tr.load_orfs([first_orf, second_orf, third_orf])

        self.assertTrue(self.tr.is_complete)
        self.tr.finalize()
        self.assertEqual(self.tr.number_internal_orfs, 2, (
            self.tr.cdna_length, self.tr.selected_start_distance_from_tss, self.tr.selected_end_distance_from_tes))

        self.assertEqual(self.tr.combined_cds_length, 648)
        self.assertEqual(self.tr.selected_cds_length, 348)
        self.assertEqual(self.tr.number_internal_orfs, 2, "\n".join([str(x) for x in self.tr.internal_orfs]))

        new_transcripts = sorted(self.tr.split_by_cds())

        self.assertEqual(len(new_transcripts), 2)
        self.assertEqual(new_transcripts[0].three_utr_length, 0)
        self.assertEqual(new_transcripts[1].five_utr_length, 0)

    def testDoubleOrf_negative(self):

        """Test to verify the introduction of multiple ORFs."""

        self.tr.strip_cds()
        self.tr.finalized = False

        first_orf = mikado_lib.parsers.bed12.BED12()
        first_orf.chrom = self.tr.id
        first_orf.start = 1
        first_orf.end = self.tr.cdna_length
        first_orf.name = self.tr.id
        first_orf.strand = "-"
        first_orf.score = 0
        first_orf.thickStart = 51
        first_orf.thickEnd = 398
        first_orf.blockCount = 1
        first_orf.blockSize = self.tr.cdna_length
        first_orf.blockSizes = [self.tr.cdna_length]
        first_orf.blockStarts = [0]
        first_orf.rgb = 0
        first_orf.has_start_codon = True
        first_orf.has_stop_codon = True
        first_orf.transcriptomic = True
        self.assertFalse(first_orf.invalid)
        # This should not be incorporated
        second_orf = mikado_lib.parsers.bed12.BED12()
        second_orf.chrom = self.tr.id
        second_orf.start = 1
        second_orf.end = self.tr.cdna_length
        second_orf.name = "second"
        second_orf.strand = "-"
        second_orf.score = 0
        second_orf.thickStart = 201
        second_orf.thickEnd = 410
        second_orf.blockCount = 1
        second_orf.blockSize = self.tr.cdna_length
        second_orf.blockSizes = [self.tr.cdna_length]
        second_orf.blockStarts = [0]
        second_orf.rgb = 0
        second_orf.has_start_codon = True
        second_orf.has_stop_codon = True
        second_orf.transcriptomic = True
        self.assertFalse(second_orf.invalid)

        self.assertTrue(mikado_lib.loci_objects.transcript.Transcript.is_overlapping_cds(first_orf, second_orf))

        # This should be added
        third_orf = mikado_lib.parsers.bed12.BED12()
        third_orf.chrom = self.tr.id
        third_orf.start = 1
        third_orf.end = self.tr.cdna_length
        third_orf.name = "third"
        third_orf.strand = "-"
        third_orf.score = 0
        third_orf.thickStart = 501
        third_orf.thickEnd = 800
        third_orf.blockCount = 1
        third_orf.blockSize = self.tr.cdna_length
        third_orf.blockSizes = [self.tr.cdna_length]
        third_orf.blockStarts = [0]
        third_orf.rgb = 0
        third_orf.has_start_codon = True
        third_orf.has_stop_codon = True
        third_orf.transcriptomic = True
        self.assertFalse(third_orf.invalid)

        self.assertFalse(mikado_lib.loci_objects.transcript.Transcript.is_overlapping_cds(first_orf, third_orf))
        self.assertFalse(mikado_lib.loci_objects.transcript.Transcript.is_overlapping_cds(second_orf, third_orf))

        self.assertFalse(third_orf == second_orf)
        self.assertFalse(first_orf == second_orf)
        self.assertFalse(first_orf == third_orf)

        candidates = [first_orf, second_orf, third_orf]

        # self.assertEqual(len(mikado_lib.loci_objects.transcript.Transcript.find_overlapping_cds(candidates)), 2)

        import logging
        handler = logging.StreamHandler()
        handler.setLevel("DEBUG")
        logger = logging.getLogger("test")
        logger.setLevel("DEBUG")
        logger.addHandler(handler)
        logger.debug("TEST")
        logger.propagate = False
        self.tr.logger = logger

        self.tr.load_orfs([first_orf, second_orf, third_orf])

        self.assertTrue(self.tr.is_complete)
        self.tr.finalize()
        self.assertEqual(self.tr.number_internal_orfs, 2, (
            self.tr.cdna_length, self.tr.selected_start_distance_from_tss, self.tr.selected_end_distance_from_tes))

        self.assertEqual(self.tr.combined_cds_length, 648)
        self.assertEqual(self.tr.selected_cds_length, 348)
        self.assertEqual(self.tr.number_internal_orfs, 2, "\n".join([str(x) for x in self.tr.internal_orfs]))

        new_transcripts = sorted(self.tr.split_by_cds())

        self.assertEqual(len(new_transcripts), 2)
        self.assertEqual(new_transcripts[0].five_utr_length, 0)
        self.assertEqual(new_transcripts[1].three_utr_length, 0)


unittest.main()
