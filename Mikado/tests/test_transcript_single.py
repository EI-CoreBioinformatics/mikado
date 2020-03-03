# coding: utf-8

"""
Unit test for monoexonic transcripts.
"""

import operator
import re
import unittest
from .. import exceptions, loci, parsers, transcripts
from ..loci import Transcript
from ..utilities.log_utils import create_null_logger, create_default_logger


class TranscriptTester(unittest.TestCase):
    tr_gff = """Chr1    TAIR10    mRNA    5928    8737    .    .    .    ID=AT1G01020.1;Parent=AT1G01020;Name=AT1G01020.1;Index=1
Chr1    TAIR10    exon    5928    8737    .    .    .    Parent=AT1G01020.1"""

    tr_lines = tr_gff.split("\n")
    for pos, line in enumerate(tr_lines):
        tr_lines[pos] = re.sub(r"\s+", "\t", line)
        assert len(tr_lines[pos].split("\t")) == 9, line.split("\t")

    tr_gff_lines = [parsers.GFF.GffLine(line) for line in tr_lines]

    for l in tr_gff_lines:
        assert l.header is False

    logger = create_null_logger("null")

    def setUp(self):
        """Basic creation test."""

        self.tr = Transcript()
        self.tr.logger = self.logger
        self.tr.chrom = "Chr1"
        self.tr.source = "TAIR10"
        self.tr.feature = "mRNA"
        self.tr.start = 5928
        self.tr.end = 8737
        self.tr.strand = "+"
        self.tr.add_exon((5928, 8737))
        self.tr.score = None
        self.tr.id, self.tr.parent, self.tr.name = "AT1G01020.1", "AT1G01020", "AT1G01020.1"
        self.tr.add_exon((8571, 8666), "CDS")
        self.tr.finalize()

        self.orf = parsers.bed12.BED12()
        self.orf.chrom = self.tr.id
        self.orf.start = 1
        self.orf.end = self.tr.cdna_length
        self.orf.name = self.tr.id
        self.orf.strand = "+"
        self.orf.score = 0
        self.orf.thick_start = 8571 - 5928 + 1
        self.orf.thick_end = 8666 - 5928 + 1
        self.orf.block_count = 1
        self.orf.blockSize = self.tr.cdna_length
        self.orf.block_starts = [0]
        self.orf.has_start_codon = True
        self.orf.has_stop_codon = True
        self.orf.transcriptomic = True
        self.assertFalse(self.orf.invalid, self.orf.invalid_reason)
        self.assertEqual((self.orf.thick_end - self.orf.thick_start + 1) % 3, 0)

    def test_basics(self):

        self.assertEqual(self.tr.chrom, "Chr1")
        self.assertEqual(self.tr.exon_num, 1)
        self.assertEqual(self.tr.monoexonic, True)
        self.assertEqual(self.tr.exon_num, len(self.tr.exons))
        self.assertEqual(self.tr.start, 5928)
        self.assertEqual(self.tr.end, 8737)
        self.assertEqual(self.tr.exons,
                         [tuple([5928, 8737])],
                         self.tr.exons)

    def test_cds(self):
        """Test the CDS features.
        Note that in a single-exon transcript with no strand, start_codon and stop_codon are defined as False.
        """

        self.tr.load_orfs([self.orf])
        self.assertEqual(self.tr.combined_cds, self.tr.selected_cds)

        self.assertEqual(self.tr.combined_cds,
                         [tuple([8571, 8666])],
                         self.tr.combined_cds)
        self.assertEqual(self.tr.selected_cds_start, 8571)
        self.assertEqual(self.tr.selected_cds_end, 8666)
        self.assertEqual(self.tr.has_start_codon, True)
        self.assertEqual(self.tr.has_stop_codon, True)

    def test_equality(self):

        new_transcript = self.tr.deepcopy()

        self.assertTrue(new_transcript == self.tr)

        new_transcript.strand = None
        self.assertFalse(new_transcript == self.tr)  # They have now a different strand

        new_transcript.unfinalize()
        new_transcript.strand = "+"  # It becomes a multiexonic transcript, so it must have a strand
        new_transcript.end = 9737

        new_exon = parsers.GFF.GffLine(self.tr_lines[-1])
        new_exon.strand = "+"
        new_exon.start = 9000
        new_exon.end = 9737
        new_transcript.add_exon(new_exon)

        new_transcript.finalize()
        self.assertTrue(new_transcript != self.tr)

    def test_mono_finalising(self):

        transcript_line = [line for line in self.tr_gff_lines if line.feature == "mRNA" ]
        self.assertEqual(len(transcript_line), 1,
                         "\n".join([str(line) for line in self.tr_gff_lines]))

        tr = loci.Transcript(transcript_line[0])
        exon_lines = [line for line in self.tr_gff_lines if
                      line.is_exon is True and "UTR" not in line.feature.upper()]
        tr.add_exons(exon_lines)
        tr.add_exon((8571, 8666), "CDS")

        tr.finalize()
        self.assertGreater(tr.three_utr_length, 0)
        self.assertGreater(tr.five_utr_length, 0)

    def test_invalid_transcript(self):
        lines = """Chr1\tTAIR10\tmRNA\t5928\t8737\t.\t.\t.\tID=AT1G01020.1;Parent=AT1G01020;Name=AT1G01020.1;Index=1
Chr1\tTAIR10\tCDS\t8571\t7500\t.\t.\t0\tParent=AT1G01020.1;
Chr1\tTAIR10\tCDS\t7503\t8666\t.\t.\t0\tParent=AT1G01020.1;
Chr1\tTAIR10\texon\t5928\t8737\t.\t.\t.\tParent=AT1G01020.1"""

        gff_lines = [parsers.GFF.GffLine(line) for line in lines.split("\n")]
        self.assertIsInstance(gff_lines[0], parsers.GFF.GffLine)
        checker = False
        if gff_lines[0].feature.endswith("transcript") or "RNA" in gff_lines[0].feature.upper():
            checker = True
        self.assertTrue(checker)
        self.assertTrue(gff_lines[0].is_transcript)
        transcript = loci.Transcript(gff_lines[0])

        transcript.logger = self.logger
        transcript.add_exons(gff_lines[1:])

        with self.assertRaises(exceptions.InvalidCDS):
            transcripts.transcript_methods.finalizing._check_cdna_vs_utr(transcript)

    def test_utr(self):

        self.assertEqual(self.tr.selected_internal_orf,
                         [("UTR", tuple([5928, 8570])),
                          ("exon", tuple([5928, 8737])),
                          ("CDS", tuple([8571, 8666]), 0),
                          ("UTR", tuple([8667, 8737]))],
                         "Right: {0}\nFound{1}".format([("UTR", 5928, 8570), ("CDS", 8571, 8666), ("UTR", 8667, 8737)],
                                                       self.tr.selected_internal_orf))
        self.assertEqual(self.tr.combined_utr, [tuple([5928, 8570]),
                                                tuple([8667, 8737])])
        self.assertEqual(self.tr.five_utr, [tuple([5928, 8570])],
                         self.tr.five_utr)
        self.assertEqual(self.tr.three_utr, [tuple([8667, 8737])])

    def test_utr_metrics(self):

        """Test for UTR exon num, start distance, etc."""

        self.assertEqual(self.tr.five_utr_num, 1)
        self.assertEqual(self.tr.three_utr_num, 1)
        self.assertEqual(self.tr.five_utr_length, 8570 + 1 - 5928)
        self.assertEqual(self.tr.three_utr_length, 8737 + 1 - 8667)
        self.assertEqual(self.tr.selected_start_distance_from_tss,
                         8571 - 5928,
                         self.tr.selected_end_distance_from_tes)
        self.assertEqual(self.tr.selected_end_distance_from_tes,
                         8737 - 8666,
                         (self.tr.selected_end_distance_from_tes, self.tr.strand))

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
                         [tuple([8571, 8666])],
                         self.tr.combined_cds)
        self.assertEqual(self.tr.combined_utr, [], self.tr.combined_utr)

    def test_negative_orf(self):
        """Test loading a negative strand ORF onto a monoexonic transcript.
        This should reverse the ORF."""

        self.orf.strand = "-"
        self.tr.strip_cds(strand_specific=False)
        self.orf.has_stop_codon = False
        self.tr.load_orfs([self.orf])
        self.assertEqual(self.tr.strand, "-")
        self.assertEqual(self.tr.selected_cds_start, 8666)
        self.assertEqual(self.tr.selected_cds_end, 8571)

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

        first_orf = parsers.bed12.BED12()
        first_orf.chrom = self.tr.id
        first_orf.start = 1
        first_orf.end = self.tr.cdna_length
        first_orf.name = self.tr.id
        first_orf.strand = "+"
        first_orf.score = 0
        first_orf.thick_start = 51
        first_orf.thick_end = 398
        first_orf.block_count = 1
        first_orf.blockSize = self.tr.cdna_length
        first_orf.block_sizes = [self.tr.cdna_length]
        first_orf.block_starts = [0]
        first_orf.rgb = 0
        first_orf.has_start_codon = True
        first_orf.has_stop_codon = True
        first_orf.transcriptomic = True
        self.assertFalse(first_orf.invalid)

        # This should not be incorporated
        second_orf = parsers.bed12.BED12()
        second_orf.chrom = self.tr.id
        second_orf.start = 1
        second_orf.end = self.tr.cdna_length
        second_orf.name = "second"
        second_orf.strand = "+"
        second_orf.score = 0
        second_orf.thick_start = 201
        second_orf.thick_end = 410
        second_orf.block_count = 1
        second_orf.blockSize = self.tr.cdna_length
        second_orf.block_sizes = [self.tr.cdna_length]
        second_orf.block_starts = [0]
        second_orf.rgb = 0
        second_orf.has_start_codon = True
        second_orf.has_stop_codon = True
        second_orf.transcriptomic = True
        self.assertFalse(second_orf.invalid)

        self.assertTrue(loci.Transcript.is_overlapping_cds(
            first_orf, second_orf))

        # This should be added
        third_orf = parsers.bed12.BED12()
        third_orf.chrom = self.tr.id
        third_orf.start = 1
        third_orf.end = self.tr.cdna_length
        third_orf.name = "third"
        third_orf.strand = "+"
        third_orf.score = 0
        third_orf.thick_start = 501
        third_orf.thick_end = 800
        third_orf.block_count = 1
        third_orf.blockSize = self.tr.cdna_length
        third_orf.block_sizes = [self.tr.cdna_length]
        third_orf.block_starts = [0]
        third_orf.rgb = 0
        third_orf.has_start_codon = True
        third_orf.has_stop_codon = True
        third_orf.transcriptomic = True
        self.assertFalse(third_orf.invalid)

        self.assertFalse(loci.Transcript.is_overlapping_cds(first_orf, third_orf))
        self.assertFalse(loci.Transcript.is_overlapping_cds(second_orf, third_orf))

        self.assertFalse(third_orf == second_orf)
        self.assertFalse(first_orf == second_orf)
        self.assertFalse(first_orf == third_orf)

        candidates = [first_orf, second_orf, third_orf]

        self.tr.logger = self.logger

        self.tr.load_orfs([first_orf])
        self.tr.load_orfs([second_orf])
        self.tr.load_orfs([third_orf])

        self.tr.load_orfs([first_orf, second_orf, third_orf])

        self.assertTrue(self.tr.is_complete)
        self.tr.finalize()
        self.assertEqual(self.tr.number_internal_orfs, 2,
                         (self.tr.cdna_length, self.tr.selected_start_distance_from_tss,
                          self.tr.selected_end_distance_from_tes))

        self.assertEqual(self.tr.combined_cds_length, 648)
        self.assertEqual(self.tr.selected_cds_length, 348)
        self.assertEqual(self.tr.number_internal_orfs, 2,
                         "\n".join([str(x) for x in self.tr.internal_orfs]))

        new_transcripts = sorted(self.tr.split_by_cds())

        self.assertEqual(len(new_transcripts), 2)
        self.assertEqual(new_transcripts[0].three_utr_length, 0)
        self.assertEqual(new_transcripts[1].five_utr_length, 0)

    def testDoubleOrf_negative(self):

        """Test to verify the introduction of multiple ORFs."""

        self.tr.strip_cds(strand_specific=False)
        self.tr.finalized = False

        first_orf = parsers.bed12.BED12()
        first_orf.chrom = self.tr.id
        first_orf.start = 1
        first_orf.end = self.tr.cdna_length
        first_orf.name = self.tr.id
        first_orf.strand = "-"
        first_orf.score = 0
        first_orf.thick_start = 51
        first_orf.thick_end = 398
        first_orf.block_count = 1
        first_orf.blockSize = self.tr.cdna_length
        first_orf.block_sizes = [self.tr.cdna_length]
        first_orf.block_starts = [0]
        first_orf.rgb = 0
        first_orf.has_start_codon = True
        first_orf.has_stop_codon = True
        first_orf.transcriptomic = True
        self.assertFalse(first_orf.invalid)

        # This should not be incorporated
        second_orf = parsers.bed12.BED12()
        second_orf.chrom = self.tr.id
        second_orf.start = 1
        second_orf.end = self.tr.cdna_length
        second_orf.name = "second"
        second_orf.strand = "-"
        second_orf.score = 0
        second_orf.thick_start = 201
        second_orf.thick_end = 410
        second_orf.block_count = 1
        second_orf.blockSize = self.tr.cdna_length
        second_orf.block_sizes = [self.tr.cdna_length]
        second_orf.block_starts = [0]
        second_orf.rgb = 0
        second_orf.has_start_codon = True
        second_orf.has_stop_codon = True
        second_orf.transcriptomic = True

        self.assertFalse(second_orf.invalid)

        # self.assertTrue(loci.Transcript.is_overlapping_cds(first_orf,
        #                                                           second_orf))

        # This should be added
        third_orf = parsers.bed12.BED12()
        third_orf.chrom = self.tr.id
        third_orf.start = 1
        third_orf.end = self.tr.cdna_length
        third_orf.name = "third"
        third_orf.strand = "-"
        third_orf.score = 0
        third_orf.thick_start = 501
        third_orf.thick_end = 800
        third_orf.block_count = 1
        third_orf.blockSize = self.tr.cdna_length
        third_orf.block_sizes = [self.tr.cdna_length]
        third_orf.block_starts = [0]
        third_orf.rgb = 0
        third_orf.has_start_codon = True
        third_orf.has_stop_codon = True
        third_orf.transcriptomic = True
        self.assertFalse(third_orf.invalid)

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

        # self.assertEqual(len(self.tr.find_overlapping_cds(candidates)), 2)

        self.tr.logger = self.logger

        self.tr.load_orfs(candidates)

        self.assertTrue(self.tr.is_complete)
        self.tr.finalize()
        self.assertEqual(self.tr.number_internal_orfs, 2, (
            self.tr.cdna_length,
            self.tr.selected_start_distance_from_tss,
            self.tr.selected_end_distance_from_tes))

        # self.assertEqual(self.tr.combined_cds_length, 648)
        self.assertEqual(self.tr.selected_cds_length, 348)
        self.assertEqual(self.tr.number_internal_orfs, 2, "\n".join([str(x) for x in self.tr.internal_orfs]))

        new_transcripts = sorted(self.tr.split_by_cds())

        self.assertEqual(len(new_transcripts), 2)
        self.assertEqual(new_transcripts[0].five_utr_length, 0)
        self.assertEqual(new_transcripts[1].three_utr_length, 0)

    def test_wrong_orf(self):
        # This should be added
        orf = parsers.bed12.BED12()
        orf.chrom = self.tr.id
        orf.start = 1
        orf.end = self.tr.cdna_length + 1
        orf.name = "third"
        orf.strand = "-"
        orf.score = 0
        orf.thick_start = 501
        orf.thick_end = 800
        orf.block_count = 1
        orf.blockSize = self.tr.cdna_length
        orf.block_sizes = [self.tr.cdna_length]
        orf.block_starts = [0]
        orf.rgb = 0
        orf.has_start_codon = True
        orf.has_stop_codon = True
        orf.transcriptomic = True
        self.assertFalse(orf.invalid)

        self.tr.logger = self.logger
        self.tr.strip_cds()
        self.tr.strand = "+"
        self.logger.setLevel("WARNING")
        # self.tr.load_orfs([orf])
        with self.assertLogs("null", level="DEBUG") as cm_out:
            self.tr.load_orfs([orf])

        self.assertFalse(self.tr.is_coding)


class TestWheatRNA(unittest.TestCase):

    def test_negative_orf(self):

        transcript = loci.Transcript()
        transcript.chrom = "Triticum_aestivum_CS42_TGACv1_scaffold_018953_1AS"
        transcript.strand = None
        transcript.start = 215963
        transcript.end = 217518
        transcript.id = "TGAC_Root_Cufflinks_CL_Root.3672.1"
        transcript.add_exons([(215963, 217518)])
        transcript.parent = "foo"
        transcript.score = 23
        transcript.finalize()

        bed = parsers.bed12.BED12()
        bed.chrom = "TGAC_Root_Cufflinks_CL_Root.3672.1"
        bed.transcriptomic = True
        bed.start = 1
        bed.end = 1556
        bed.name = "ID=TGAC_Root_Cufflinks_CL_Root.3672.1|m.294;TGAC_Root_Cufflinks_CL_Root.3672.1|g.294;ORF_TGAC_Root_Cufflinks_CL_Root.3672.1|g.294_TGAC_Root_Cufflinks_CL_Root.3672.1|m.294_type:5prime_partial_len:157_(-)"
        bed.strand = "-"
        bed.thick_start = 1084
        bed.thick_end = 1554
        bed.has_start_codon = False
        bed.has_stop_codon = True
        bed.header = False
        self.assertFalse(bed.invalid)

        bed2 = parsers.bed12.BED12()
        bed2.chrom = "TGAC_Root_Cufflinks_CL_Root.3672.1"
        bed2.transcriptomic = True
        bed2.start = 1
        bed2.end = 1556
        bed2.name = "ID=TGAC_Root_Cufflinks_CL_Root.3672.1|name2"
        bed2.strand = "-"
        bed2.thick_start = 884
        bed2.thick_end = 1054
        bed2.has_start_codon = False
        bed2.has_stop_codon = False
        bed2.header = False
        self.assertFalse(bed2.invalid)

        transcript.load_orfs([bed, bed2])
        self.assertEqual(len(transcript.internal_orfs), 2)
        transcript.finalize()
        self.assertEqual(transcript.start, 215963)
        self.assertEqual(transcript.end, 217518)
        self.assertEqual(transcript.strand, "-")
        self.assertEqual((transcript.start,
                          transcript.selected_cds_start,
                          transcript.selected_cds_end,
                          transcript.end),
                         (215963,
                          217516,
                          217046,
                          217518
                          ))

    def test_orf_sorter(self):
        bed = parsers.bed12.BED12()
        bed.chrom = "TGAC_Root_Cufflinks_CL_Root.3672.1"
        bed.transcriptomic = True
        bed.start = 1
        bed.end = 1556
        bed.name = "ID=TGAC_Root_Cufflinks_CL_Root.3672.1|m.294;TGAC_Root_Cufflinks_CL_Root.3672.1|g.294;ORF_TGAC_Root_Cufflinks_CL_Root.3672.1|g.294_TGAC_Root_Cufflinks_CL_Root.3672.1|m.294_type:5prime_partial_len:157_(-)"
        bed.strand = "-"
        bed.thick_start = 1084
        bed.thick_end = 1554
        bed.has_start_codon = False
        bed.has_stop_codon = True
        self.assertFalse(bed.invalid)

        bed2 = parsers.bed12.BED12()
        bed2.chrom = "TGAC_Root_Cufflinks_CL_Root.3672.1"
        bed2.transcriptomic = True
        bed2.start = 1
        bed2.end = 1556
        bed2.name = "ID=TGAC_Root_Cufflinks_CL_Root.3672.1|2"
        bed2.strand = "-"
        bed2.thick_start = 884
        bed2.thick_end = 1054
        bed2.has_start_codon = False
        bed2.has_stop_codon = False

        after_sorting = sorted([bed, bed2],
                               reverse=True,
                               key=transcripts.transcript_methods.retrieval.orf_sorter
                               )

        self.assertEqual(after_sorting[0], bed)


class TestNegativeSplit(unittest.TestCase):

    def setUp(self):

        tr = """Triticum_aestivum_CS42_TGACv1_scaffold_018974_1AS\tCufflinks\ttranscript\t72914\t76276\t1000\t.\t.\tgene_id "CL_10DPA.2184"; transcript_id "ERP004505_10DPA_Cufflinks_CL_10DPA.2184.1"; exon_number "1"; FPKM "0.1596439944"; conf_hi "0.205257"; conf_lo "0.114031"; frac "1.000000"; cov "3.066911"; Name "ERP004505_10DPA_Cufflinks_CL_10DPA.2184.1";
Triticum_aestivum_CS42_TGACv1_scaffold_018974_1AS\tCufflinks\texon\t72914\t76276\t.\t.\t.\tgene_id "CL_10DPA.2184"; transcript_id "ERP004505_10DPA_Cufflinks_CL_10DPA.2184.1";"""
        tr_lines = [parsers.GTF.GtfLine(_) for _ in tr.split("\n")]
        self.tr = loci.Transcript(tr_lines[0])
        self.tr.add_exon(tr_lines[1])

        self.bed1 = parsers.bed12.BED12()
        self.bed1.header = False
        self.bed1.chrom = "ERP004505_10DPA_Cufflinks_CL_10DPA.2184.1"
        self.bed1.start = 1
        self.bed1.end = 3363
        self.bed1.name = "First|m.2472_type:complete_len:193_(-)"
        self.bed1.strand = "-"
        self.bed1.thick_start = 1423
        self.bed1.thick_end = 2001
        self.bed1.has_stop_codon = True
        self.bed1.has_start_codon = True
        self.bed1.transcriptomic = True
        self.assertEqual(self.bed1.cds_len, 579)
        
        self.bed2 = parsers.bed12.BED12()
        self.bed2.header = False
        self.bed2.chrom = "ERP004505_10DPA_Cufflinks_CL_10DPA.2184.1"
        self.bed2.start = 1
        self.bed2.end = 3363
        self.bed2.name = "Second|m.2472_type:complete_len:137_(-)"
        self.bed2.strand = "-"
        self.bed2.thick_start = 2481
        self.bed2.thick_end = 2891
        self.bed2.has_stop_codon = True
        self.bed2.has_start_codon = True
        self.bed2.transcriptomic = True
        self.assertEqual(self.bed2.cds_len, 411)

    def test_loading(self):
        
        self.tr.load_orfs([self.bed1, self.bed2])
        self.assertEqual(self.tr.number_internal_orfs, 2)
        self.assertEqual(self.tr.selected_cds_start,
                         74914)
        self.assertEqual(self.tr.selected_cds_end,
                         74336)
        self.assertEqual(self.tr.combined_cds_end,
                         74336)
        self.assertEqual(self.tr.combined_cds_start,
                         75804)

    def test_split(self):
        self.tr.load_orfs([self.bed1, self.bed2])
        self.assertEqual(self.tr.number_internal_orfs, 2)
        logger = create_default_logger("splitter")
        logger.setLevel("ERROR")
        self.tr.logger = logger

        new_transcripts = [_ for _ in self.tr.split_by_cds()]

        new_transcripts = sorted(new_transcripts,
                                 key=operator.attrgetter("start", "end"))
        self.assertEqual(len(new_transcripts), 2)

        self.assertEqual(new_transcripts[0].start, 72914)
        self.assertEqual(new_transcripts[0].end, 74914)

        self.assertEqual(new_transcripts[1].end, 76276, self.tr.internal_orfs)
        self.assertEqual(new_transcripts[1].start, 75394)

        self.assertEqual(new_transcripts[0].selected_cds_start, 74914)
        self.assertEqual(new_transcripts[0].selected_cds_end, 74336)

        self.assertEqual(new_transcripts[1].selected_cds_start, 75804)
        self.assertEqual(new_transcripts[1].selected_cds_end, 75394)


class TestTripleNegative(unittest.TestCase):

    def setUp(self):

        trlines = """Triticum_aestivum_CS42_TGACv1_scaffold_019715_1AS\tCufflinks\ttranscript\t44187\t49369\t1000\t.\t.\tgene_id "CL_Spike.6733"; transcript_id "TGAC_Spike_Cufflinks_CL_Spike.6733.1"; exon_number "1"; cov "2963.838405"; frac "0.338459"; FPKM "46.8523047496"; conf_hi "48.870099"; Name "TGAC_Spike_Cufflinks_CL_Spike.6733.1"; conf_lo "44.834511";
Triticum_aestivum_CS42_TGACv1_scaffold_019715_1AS\tCufflinks\texon\t44187\t49369\t.\t.\t.\tgene_id "CL_Spike.6733"; transcript_id "TGAC_Spike_Cufflinks_CL_Spike.6733.1";"""
        trlines = [parsers.GTF.GtfLine(_) for _ in trlines.split("\n")]

        self.tr = loci.Transcript(trlines[0])
        self.tr.add_exon(trlines[1])
        self.tr.finalize()

        self.bed1 = parsers.bed12.BED12(transcriptomic=True)
        self.bed1.header = False
        self.bed1.chrom = self.tr.id
        self.bed1.start = 1
        self.bed1.end = 5183
        self.bed1.strand = "-"
        self.bed1.name = "g.132606_TGAC_Spike_Cufflinks_CL_Spike.6733.1|m.132606_type:3prime_partial_len:144_(-)"
        self.bed1.thick_start = 2
        self.bed1.thick_end = 430
        self.bed1.score = 0
        # self.bed1.has_start_codon = True
        # self.bed1.has_stop_codon = False
        self.assertEqual(self.bed1.cds_len, 429)
        
        self.bed2 = parsers.bed12.BED12(transcriptomic=True)
        self.bed2.header = False
        self.bed2.chrom = self.tr.id
        self.bed2.start = 1
        self.bed2.end = 5183
        self.bed2.strand = "-"
        self.bed2.name = "g.132602_TGAC_Spike_Cufflinks_CL_Spike.6733.1|m.132602_type:complete_len:387_(-)"
        self.bed2.thick_start = 936
        self.bed2.thick_end = 2096
        self.bed2.score = 0
        # self.bed2.has_stop_codon = True
        # self.bed2.has_start_codon = True
        self.assertEqual(self.bed2.cds_len, 1161)
        
        self.bed3 = parsers.bed12.BED12(transcriptomic=True)
        self.bed3.header = False
        self.bed3.chrom = self.tr.id
        self.bed3.start = 1
        self.bed3.end = 5183
        self.bed3.strand = "-"
        self.bed3.name = "g.132601_TGAC_Spike_Cufflinks_CL_Spike.6733.1|m.132601_type:5prime_partial_len:542_(-)"
        # self.bed3.has_start_codon = False
        # self.bed3.has_stop_codon = True
        self.bed3.thick_start = 3558
        self.bed3.thick_end = 5183
        self.bed3.score = 0
        self.assertEqual(self.bed3.cds_len, 1626)

    def test_split(self):

        self.tr.load_orfs([self.bed1, self.bed2, self.bed3])

        new = sorted([_ for _ in self.tr.split_by_cds()], key=operator.attrgetter("start", "end"))
        self.assertEqual(len(new), 3)
        self.assertEqual(new[0].start, 44187)
        self.assertEqual(new[0].end, 44616)
        self.assertEqual(new[0].selected_cds_end, 44188)

        self.assertEqual(new[1].start, new[1].selected_cds_end)
        self.assertEqual(new[1].start, 45122)
        self.assertEqual(new[1].end, new[1].selected_cds_start)
        self.assertEqual(new[1].end, 46282)

        self.assertEqual(new[2].start, 47744)
        self.assertEqual(new[2].end, self.tr.end)
        self.assertEqual(new[2].end, new[2].selected_cds_start)


class TestNegativeMonoexonicLoad(unittest.TestCase):

    def setUp(self):

        trlines="""Chr3\tCufflinks\ttranscript\t2949168\t2952410\t1000\t-\t.\tgene_id "cufflinks_star_at.10687"; transcript_id "cufflinks_cufflinks_star_at.10687.1"; exon_number "1"; conf_lo "0.169545"; FPKM "0.1111324743"; frac "0.781250"; conf_hi "0.308263"; cov "3.129756";
Chr3\tCufflinks\texon\t2949168\t2952410\t.\t-\t.\tgene_id "cufflinks_star_at.10687"; transcript_id "cufflinks_cufflinks_star_at.10687.1";"""

        trlines = [parsers.GTF.GtfLine(_) for _ in trlines.split("\n")]

        self.tr = loci.Transcript(trlines[0])
        self.tr.add_exon(trlines[1])
        self.tr.finalize()

        self.bed1 = parsers.bed12.BED12()
        self.bed1.header = False
        self.bed1.chrom = self.tr.id
        self.bed1.start = 1
        self.bed1.end = 3243
        self.bed1.transcriptomic = True
        self.bed1.thick_start = 203
        self.bed1.thick_end = 2206
        self.bed1.name = "g.135538_cufflinks_cufflinks_star_at.10687.1|m.135538_type:complete_len:668_(+)"
        self.bed1.has_start_codon = True
        self.bed1.has_stop_codon = True
        self.bed1.strand = "+"
        self.assertEqual(self.bed1.cds_len, 2004)
        self.bed1.score = 0
        self.assertFalse(self.bed1.invalid)

        self.bed2 = parsers.bed12.BED12()
        self.bed2.header = False
        self.bed2.chrom = self.tr.id
        self.bed2.start = 1
        self.bed2.end = 3243
        self.bed2.transcriptomic = True
        self.bed2.thick_start = 2543
        self.bed2.thick_end = 3241
        self.bed2.strand = "+"
        self.bed2.name = "g.135539_cufflinks_cufflinks_star_at.10687.1|m.135539_type:3prime_partial_len:234_(+)"
        self.bed2.has_start_codon = True
        self.bed2.has_stop_codon = False
        self.assertEqual(self.bed2.cds_len, 699)
        self.bed2.score = 0
        self.assertFalse(self.bed2.invalid)

    def test_load(self):

        self.tr.load_orfs([self.bed1, self.bed2])
        self.assertEqual(self.tr.strand, "-")

        self.assertEqual(self.tr.number_internal_orfs, 2)
        self.assertEqual(self.tr.combined_cds_end, 2949170)
        self.assertEqual(self.tr.combined_cds_start, 2952208)
        self.assertEqual(self.tr.selected_cds_start, 2952208)
        self.assertEqual(self.tr.selected_cds_end, 2950205)

    # @unittest.skip
    def test_split(self):
# Chr3	Cufflinks	mRNA	2949168	2952410	1000	-	.	ID=cufflinks_cufflinks_star_at.10687.1.orf1;
# Chr3	Cufflinks	three_prime_UTR	2949168	2950204	.	-	.	Parent=cufflinks_cufflinks_star_at.10687.1.orf1
# Chr3	Cufflinks	exon	2949168	2952410	.	-	.	Parent=cufflinks_cufflinks_star_at.10687.1.orf1
# Chr3	Cufflinks	CDS	2,950,205	2,952,208	.	-	0	Parent=cufflinks_cufflinks_star_at.10687.1.orf1
# Chr3	Cufflinks	five_prime_UTR	2952209	2952410	.	-	.	Parent=cufflinks_cufflinks_star_at.10687.1.orf1
#
# Chr3	Cufflinks	mRNA	2949168	2952410	1000	-	.	ID=cufflinks_cufflinks_star_at.10687.1.orf2;
# Chr3	Cufflinks	three_prime_UTR	2949168	2949169	.	-	.	Parent=cufflinks_cufflinks_star_at.10687.1.orf2
# Chr3	Cufflinks	exon	2949168	2952410	.	-	.	Parent=cufflinks_cufflinks_star_at.10687.1.orf2
# Chr3	Cufflinks	CDS	2,949,170	2,949,868	.	-	0	Parent=cufflinks_cufflinks_star_at.10687.1.orf2
# Chr3	Cufflinks	five_prime_UTR	2949869	2952410	.	-	.	Parent=cufflinks_cufflinks_star_at.10687.1.orf2
# =============>
# Chr3	Cufflinks	mRNA	2950205	2952410	1000	-	.	ID=cufflinks_cufflinks_star_at.10687.1.split1;
# Chr3	Cufflinks	exon	2950205 2952410	.	-	.	Parent=cufflinks_cufflinks_star_at.10687.1.split1
# Chr3	Cufflinks	CDS	2950205	2952208	.	-	0	Parent=cufflinks_cufflinks_star_at.10687.1.split1
# Chr3	Cufflinks	five_prime_UTR	2952209	2952410	.	-	.	Parent=cufflinks_cufflinks_star_at.10687.1.split1
#
# Chr3	Cufflinks	mRNA	2949168	2949868	1000	-	.	ID=cufflinks_cufflinks_star_at.10687.1.split2;
# Chr3	Cufflinks	three_prime_UTR	2949168	2949169	.	-	.	Parent=cufflinks_cufflinks_star_at.10687.1.split2
# Chr3	Cufflinks	exon	2949168	2949868	.	-	.	Parent=cufflinks_cufflinks_star_at.10687.1.split2
# Chr3	Cufflinks	CDS	2949170	2949868	.	-	0	Parent=cufflinks_cufflinks_star_at.10687.1.split2

        self.tr.load_orfs([self.bed1, self.bed2])
        self.assertEqual(self.tr.strand, "-")

        self.assertEqual(self.tr.number_internal_orfs, 2)
        self.assertEqual(self.tr.combined_cds_end, 2949170)
        self.assertEqual(self.tr.combined_cds_start, 2952208)
        self.assertEqual(self.tr.selected_cds_start, 2952208)
        self.assertEqual(self.tr.selected_cds_end, 2950205)

        logger = create_default_logger("splitter")
        logger.setLevel("ERROR")
        self.tr.logger = logger

        new_transcripts = sorted([_ for _ in self.tr.split_by_cds()])

        self.assertEqual(new_transcripts[0].start, self.tr.start)
        self.assertEqual(new_transcripts[0].end, 2949868, "\n\n".join([str(_) for _ in new_transcripts]))


if __name__ == '__main__':
    unittest.main()
