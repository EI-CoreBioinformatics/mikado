# coding: utf-8

"""
Unit test for a transcript on the positive strand.
"""

import unittest
import re
import copy
from .. import parsers, exceptions, loci
from ..utilities.log_utils import create_null_logger, create_default_logger


class MonoBaseTester(unittest.TestCase):

    """
    This test verifies the correct ORF loading and splitting
     in the case where the transcript has multiple ORFs and
     in one case it starts exactly at the terminal point of
      a previous exon.
    """

    logger = create_null_logger("null")

    def setUp(self):
        self.tr = loci.Transcript()
        self.tr.chrom = "Chr5"
        self.tr.start = 22597965
        self.tr.end = 22602701
        self.tr.strand = "+"
        self.tr.score = 1000
        self.tr.parent = "StringTie_DN.70115"
        self.tr.id = "StringTie_DN.70115.4"
        self.tr.source = "StringTie"
        self.tr.feature = "transcript"
        self.tr.add_exons([(22597965, 22601782),
                           (22601862, 22601957),
                           (22602039, 22602701)])

        self.tr.logger = self.logger

        # First ORF
        self.bed1 = parsers.bed12.BED12()
        self.bed1.chrom = self.tr.id
        self.bed1.start = 1
        self.bed1.end = 4577
        self.bed1.name = "{0}.1".format(self.tr.id)
        self.bed1.strand = "+"
        self.bed1.score = 0
        self.bed1.thick_start = 434
        self.bed1.thick_end = 3736
        self.bed1.has_start_codon = True
        self.bed1.transcriptomic = True
        self.bed1.has_stop_codon = True
        self.bed1.block_count = 1
        self.bed1.block_sizes = [len(self.bed1)]
        self.bed1.block_starts = [0]

        # Second ORF
        self.bed2 = copy.deepcopy(self.bed1)
        self.bed2.name = "{0}.2".format(self.tr.id)
        self.bed2.thick_start = 2
        self.bed2.thick_end = 388
        self.bed2.has_start_codon = False

        # Third ORF
        self.bed3 = copy.deepcopy(self.bed1)
        self.bed3.name = "{0}.3".format(self.tr.id)
        self.bed3.thick_start = 3914
        self.bed3.thick_end = 4393

    def test_finalise(self):
        self.tr.finalize()
        self.assertTrue(self.tr.finalized)

        self.assertEqual(self.tr.max_exon_length, 3818)
        self.assertEqual(self.tr.min_exon_length, 96)
        self.assertEqual(self.tr.max_intron_length, 81, self.tr.introns)
        self.assertEqual(self.tr.min_intron_length, 79, self.tr.introns)

    def test_load_orfs(self):
        self.assertFalse(self.bed1.invalid)
        self.assertFalse(self.bed2.invalid)
        self.assertFalse(self.bed3.invalid)
        self.assertEqual(self.bed3.cds_len, self.bed3.thick_end-self.bed3.thick_start+1 )

        self.tr.load_orfs([self.bed1, self.bed2, self.bed3])
        self.assertEqual(self.tr.number_internal_orfs, 3)
        self.assertEqual(self.tr.selected_cds_length, self.bed1.cds_len)

    def test_split(self):

        self.tr.load_orfs([self.bed3, self.bed1])
        splitted_transcripts = [l for l in self.tr.split_by_cds()]
        self.assertEqual(len(splitted_transcripts), 2)

    def test_print(self):

        self.tr.logger = self.logger
        self.tr.finalize()
        self.maxDiff = None

        real_printed = """Chr5\tStringTie\ttranscript\t22597965\t22602701\t1000\t+\t.\tID=StringTie_DN.70115.4;Parent=StringTie_DN.70115
Chr5\tStringTie\texon\t22597965\t22601782\t.\t+\t.\tID=StringTie_DN.70115.4.exon1;Parent=StringTie_DN.70115.4
Chr5\tStringTie\texon\t22601862\t22601957\t.\t+\t.\tID=StringTie_DN.70115.4.exon2;Parent=StringTie_DN.70115.4
Chr5\tStringTie\texon\t22602039\t22602701\t.\t+\t.\tID=StringTie_DN.70115.4.exon3;Parent=StringTie_DN.70115.4"""

        self.assertEqual(str(self.tr.format("gff3")),
                         real_printed)

        real_printed_gtf = """Chr5\tStringTie\ttranscript\t22597965\t22602701\t1000\t+\t.\tgene_id "StringTie_DN.70115"; transcript_id "StringTie_DN.70115.4";
Chr5\tStringTie\texon\t22597965\t22601782\t.\t+\t.\tgene_id "StringTie_DN.70115"; transcript_id "StringTie_DN.70115.4";
Chr5\tStringTie\texon\t22601862\t22601957\t.\t+\t.\tgene_id "StringTie_DN.70115"; transcript_id "StringTie_DN.70115.4";
Chr5\tStringTie\texon\t22602039\t22602701\t.\t+\t.\tgene_id "StringTie_DN.70115"; transcript_id "StringTie_DN.70115.4";"""

        self.assertEqual(self.tr.__str__(to_gtf=True),
                         real_printed_gtf)

        pass

    def test_print_cds(self):

        self.tr.load_orfs([self.bed1])
        self.maxDiff = None

        # self.bed1.end = 4577
        # self.bed1.thick_start = 434
        # self.bed1.thick_end = 3736

        real_printed = """Chr5\tStringTie\tmRNA\t22597965\t22602701\t1000\t+\t.\tID=StringTie_DN.70115.4;Parent=StringTie_DN.70115;Name=StringTie_DN.70115.4
Chr5\tStringTie\texon\t22597965\t22601782\t.\t+\t.\tID=StringTie_DN.70115.4.exon1;Parent=StringTie_DN.70115.4
Chr5\tStringTie\tfive_prime_UTR\t22597965\t22598397\t.\t+\t.\tID=StringTie_DN.70115.4.five_prime_UTR1;Parent=StringTie_DN.70115.4
Chr5\tStringTie\tCDS\t22598398\t22601700\t.\t+\t0\tID=StringTie_DN.70115.4.CDS1;Parent=StringTie_DN.70115.4
Chr5\tStringTie\tthree_prime_UTR\t22601701\t22601782\t.\t+\t.\tID=StringTie_DN.70115.4.three_prime_UTR1;Parent=StringTie_DN.70115.4
Chr5\tStringTie\texon\t22601862\t22601957\t.\t+\t.\tID=StringTie_DN.70115.4.exon2;Parent=StringTie_DN.70115.4
Chr5\tStringTie\tthree_prime_UTR\t22601862\t22601957\t.\t+\t.\tID=StringTie_DN.70115.4.three_prime_UTR2;Parent=StringTie_DN.70115.4
Chr5\tStringTie\texon\t22602039\t22602701\t.\t+\t.\tID=StringTie_DN.70115.4.exon3;Parent=StringTie_DN.70115.4
Chr5\tStringTie\tthree_prime_UTR\t22602039\t22602701\t.\t+\t.\tID=StringTie_DN.70115.4.three_prime_UTR3;Parent=StringTie_DN.70115.4"""

        self.assertEqual(str(self.tr),
                         real_printed)

        real_printed_gtf = """Chr5\tStringTie\tmRNA\t22597965\t22602701\t1000\t+\t.\tgene_id "StringTie_DN.70115"; transcript_id "StringTie_DN.70115.4"; Name "StringTie_DN.70115.4";
Chr5\tStringTie\texon\t22597965\t22601782\t.\t+\t.\tgene_id "StringTie_DN.70115"; transcript_id "StringTie_DN.70115.4";
Chr5\tStringTie\t5UTR\t22597965\t22598397\t.\t+\t.\tgene_id "StringTie_DN.70115"; transcript_id "StringTie_DN.70115.4";
Chr5\tStringTie\tCDS\t22598398\t22601700\t.\t+\t0\tgene_id "StringTie_DN.70115"; transcript_id "StringTie_DN.70115.4";
Chr5\tStringTie\t3UTR\t22601701\t22601782\t.\t+\t.\tgene_id "StringTie_DN.70115"; transcript_id "StringTie_DN.70115.4";
Chr5\tStringTie\texon\t22601862\t22601957\t.\t+\t.\tgene_id "StringTie_DN.70115"; transcript_id "StringTie_DN.70115.4";
Chr5\tStringTie\t3UTR\t22601862\t22601957\t.\t+\t.\tgene_id "StringTie_DN.70115"; transcript_id "StringTie_DN.70115.4";
Chr5\tStringTie\texon\t22602039\t22602701\t.\t+\t.\tgene_id "StringTie_DN.70115"; transcript_id "StringTie_DN.70115.4";
Chr5\tStringTie\t3UTR\t22602039\t22602701\t.\t+\t.\tgene_id "StringTie_DN.70115"; transcript_id "StringTie_DN.70115.4";"""

        import itertools

        for lines in itertools.zip_longest(self.tr.__str__(to_gtf=True).split("\n"),
                                           real_printed_gtf.split("\n")):
            self.assertEqual(lines[0], lines[1])

        # self.assertEqual(self.tr.__str__(to_gtf=True),
        #                  real_printed_gtf.rstrip())

    def test_print_multiple_orfs(self):

        self.maxDiff = None
        self.tr.load_orfs([self.bed1, self.bed3])
        
        # self.bed1.end = 4577
        # self.bed1.thick_start = 434
        # self.bed1.thick_end = 3736
        # self.bed3.thick_start = 3914
        # self.bed3.thick_end = 4393
        
        real_printed = """Chr5\tStringTie\tmRNA\t22597965\t22602701\t1000\t+\t.\tID=StringTie_DN.70115.4.orf1;Parent=StringTie_DN.70115;Name=StringTie_DN.70115.4;maximal=True
Chr5\tStringTie\texon\t22597965\t22601782\t.\t+\t.\tID=StringTie_DN.70115.4.orf1.exon1;Parent=StringTie_DN.70115.4.orf1
Chr5\tStringTie\tfive_prime_UTR\t22597965\t22598397\t.\t+\t.\tID=StringTie_DN.70115.4.orf1.five_prime_UTR1;Parent=StringTie_DN.70115.4.orf1
Chr5\tStringTie\tCDS\t22598398\t22601700\t.\t+\t0\tID=StringTie_DN.70115.4.orf1.CDS1;Parent=StringTie_DN.70115.4.orf1
Chr5\tStringTie\tthree_prime_UTR\t22601701\t22601782\t.\t+\t.\tID=StringTie_DN.70115.4.orf1.three_prime_UTR1;Parent=StringTie_DN.70115.4.orf1
Chr5\tStringTie\texon\t22601862\t22601957\t.\t+\t.\tID=StringTie_DN.70115.4.orf1.exon2;Parent=StringTie_DN.70115.4.orf1
Chr5\tStringTie\tthree_prime_UTR\t22601862\t22601957\t.\t+\t.\tID=StringTie_DN.70115.4.orf1.three_prime_UTR2;Parent=StringTie_DN.70115.4.orf1
Chr5\tStringTie\texon\t22602039\t22602701\t.\t+\t.\tID=StringTie_DN.70115.4.orf1.exon3;Parent=StringTie_DN.70115.4.orf1
Chr5\tStringTie\tthree_prime_UTR\t22602039\t22602701\t.\t+\t.\tID=StringTie_DN.70115.4.orf1.three_prime_UTR3;Parent=StringTie_DN.70115.4.orf1
Chr5\tStringTie\tmRNA\t22597965\t22602701\t1000\t+\t.\tID=StringTie_DN.70115.4.orf2;Parent=StringTie_DN.70115;Name=StringTie_DN.70115.4;maximal=False
Chr5\tStringTie\texon\t22597965\t22601782\t.\t+\t.\tID=StringTie_DN.70115.4.orf2.exon1;Parent=StringTie_DN.70115.4.orf2
Chr5\tStringTie\tfive_prime_UTR\t22597965\t22601782\t.\t+\t.\tID=StringTie_DN.70115.4.orf2.five_prime_UTR1;Parent=StringTie_DN.70115.4.orf2
Chr5\tStringTie\texon\t22601862\t22601957\t.\t+\t.\tID=StringTie_DN.70115.4.orf2.exon2;Parent=StringTie_DN.70115.4.orf2
Chr5\tStringTie\tfive_prime_UTR\t22601862\t22601956\t.\t+\t.\tID=StringTie_DN.70115.4.orf2.five_prime_UTR2;Parent=StringTie_DN.70115.4.orf2
Chr5\tStringTie\tCDS\t22601957\t22601957\t.\t+\t0\tID=StringTie_DN.70115.4.orf2.CDS1;Parent=StringTie_DN.70115.4.orf2
Chr5\tStringTie\tCDS\t22602039\t22602517\t.\t+\t2\tID=StringTie_DN.70115.4.orf2.CDS2;Parent=StringTie_DN.70115.4.orf2
Chr5\tStringTie\texon\t22602039\t22602701\t.\t+\t.\tID=StringTie_DN.70115.4.orf2.exon3;Parent=StringTie_DN.70115.4.orf2
Chr5\tStringTie\tthree_prime_UTR\t22602518\t22602701\t.\t+\t.\tID=StringTie_DN.70115.4.orf2.three_prime_UTR1;Parent=StringTie_DN.70115.4.orf2"""

        self.assertEqual(self.tr.format("gff", all_orfs=True),
                         real_printed)

    def test_print_without_cds(self):

        self.maxDiff = None
        self.tr.load_orfs([self.bed1, self.bed3])
        real_printed = """Chr5\tStringTie\tmRNA\t22597965\t22602701\t1000\t+\t.\tID=StringTie_DN.70115.4;Parent=StringTie_DN.70115
Chr5\tStringTie\texon\t22597965\t22601782\t.\t+\t.\tID=StringTie_DN.70115.4.exon1;Parent=StringTie_DN.70115.4
Chr5\tStringTie\texon\t22601862\t22601957\t.\t+\t.\tID=StringTie_DN.70115.4.exon2;Parent=StringTie_DN.70115.4
Chr5\tStringTie\texon\t22602039\t22602701\t.\t+\t.\tID=StringTie_DN.70115.4.exon3;Parent=StringTie_DN.70115.4"""

        self.assertEqual(self.tr.format("gff3", with_cds=False),
                         real_printed)

        real_printed_gtf = """Chr5\tStringTie\tmRNA\t22597965\t22602701\t1000\t+\t.\tgene_id "StringTie_DN.70115"; transcript_id "StringTie_DN.70115.4";
Chr5\tStringTie\texon\t22597965\t22601782\t.\t+\t.\tgene_id "StringTie_DN.70115"; transcript_id "StringTie_DN.70115.4";
Chr5\tStringTie\texon\t22601862\t22601957\t.\t+\t.\tgene_id "StringTie_DN.70115"; transcript_id "StringTie_DN.70115.4";
Chr5\tStringTie\texon\t22602039\t22602701\t.\t+\t.\tgene_id "StringTie_DN.70115"; transcript_id "StringTie_DN.70115.4";"""

        self.assertEqual(self.tr.format("gtf", with_cds=False),
                         real_printed_gtf)


class DrosoTester(unittest.TestCase):

    logger = create_null_logger("droso")

    def setUp(self):

        ref_gtf = """2L\tprotein_coding\ttranscript\t523736\t540560\t.\t+\t.\tgene_id "FBgn0003963"; transcript_id "FBtr0329895"; exon_number "1"; gene_name "ush"; transcript_name "ush-RC"; exon_id "FBgn0003963:5"; gene_biotype "protein_coding";
2L\tprotein_coding\texon\t523736\t524059\t.\t+\t.\tgene_id "FBgn0003963"; transcript_id "FBtr0329895"; exon_number "1"; gene_name "ush"; transcript_name "ush-RC"; exon_id "FBgn0003963:5"; gene_biotype "protein_coding";
2L\tprotein_coding\texon\t525392\t525436\t.\t+\t.\tgene_id "FBgn0003963"; transcript_id "FBtr0329895"; exon_number "2"; gene_name "ush"; transcript_name "ush-RC"; exon_id "FBgn0003963:677"; gene_biotype "protein_coding";
2L\tprotein_coding\texon\t536023\t536966\t.\t+\t.\tgene_id "FBgn0003963"; transcript_id "FBtr0329895"; exon_number "3"; gene_name "ush"; transcript_name "ush-RC"; exon_id "FBgn0003963:7"; gene_biotype "protein_coding";
2L\tprotein_coding\texon\t537037\t537431\t.\t+\t.\tgene_id "FBgn0003963"; transcript_id "FBtr0329895"; exon_number "4"; gene_name "ush"; transcript_name "ush-RC"; exon_id "FBgn0003963:8"; gene_biotype "protein_coding";
2L\tprotein_coding\texon\t537549\t537749\t.\t+\t.\tgene_id "FBgn0003963"; transcript_id "FBtr0329895"; exon_number "5"; gene_name "ush"; transcript_name "ush-RC"; exon_id "FBgn0003963:9"; gene_biotype "protein_coding";
2L\tprotein_coding\texon\t537863\t539249\t.\t+\t.\tgene_id "FBgn0003963"; transcript_id "FBtr0329895"; exon_number "6"; gene_name "ush"; transcript_name "ush-RC"; exon_id "FBgn0003963:10"; gene_biotype "protein_coding";
2L\tprotein_coding\texon\t539310\t539452\t.\t+\t.\tgene_id "FBgn0003963"; transcript_id "FBtr0329895"; exon_number "7"; gene_name "ush"; transcript_name "ush-RC"; exon_id "FBgn0003963:11"; gene_biotype "protein_coding";
2L\tprotein_coding\texon\t539518\t540560\t.\t+\t.\tgene_id "FBgn0003963"; transcript_id "FBtr0329895"; exon_number "8"; gene_name "ush"; transcript_name "ush-RC"; exon_id "FBgn0003963:13"; gene_biotype "protein_coding";
2L\tprotein_coding\tCDS\t524038\t524059\t.\t+\t0\tgene_id "FBgn0003963"; transcript_id "FBtr0329895"; exon_number "1"; gene_name "ush"; transcript_name "ush-RC"; protein_id "FBpp0302929"; gene_biotype "protein_coding";
2L\tprotein_coding\tCDS\t525392\t525436\t.\t+\t2\tgene_id "FBgn0003963"; transcript_id "FBtr0329895"; exon_number "2"; gene_name "ush"; transcript_name "ush-RC"; protein_id "FBpp0302929"; gene_biotype "protein_coding";
2L\tprotein_coding\tCDS\t536023\t536966\t.\t+\t2\tgene_id "FBgn0003963"; transcript_id "FBtr0329895"; exon_number "3"; gene_name "ush"; transcript_name "ush-RC"; protein_id "FBpp0302929"; gene_biotype "protein_coding";
2L\tprotein_coding\tCDS\t537037\t537431\t.\t+\t0\tgene_id "FBgn0003963"; transcript_id "FBtr0329895"; exon_number "4"; gene_name "ush"; transcript_name "ush-RC"; protein_id "FBpp0302929"; gene_biotype "protein_coding";
2L\tprotein_coding\tCDS\t537549\t537749\t.\t+\t1\tgene_id "FBgn0003963"; transcript_id "FBtr0329895"; exon_number "5"; gene_name "ush"; transcript_name "ush-RC"; protein_id "FBpp0302929"; gene_biotype "protein_coding";
2L\tprotein_coding\tCDS\t537863\t539249\t.\t+\t1\tgene_id "FBgn0003963"; transcript_id "FBtr0329895"; exon_number "6"; gene_name "ush"; transcript_name "ush-RC"; protein_id "FBpp0302929"; gene_biotype "protein_coding";
2L\tprotein_coding\tCDS\t539310\t539452\t.\t+\t0\tgene_id "FBgn0003963"; transcript_id "FBtr0329895"; exon_number "7"; gene_name "ush"; transcript_name "ush-RC"; protein_id "FBpp0302929"; gene_biotype "protein_coding";
2L\tprotein_coding\tCDS\t539518\t540016\t.\t+\t1\tgene_id "FBgn0003963"; transcript_id "FBtr0329895"; exon_number "8"; gene_name "ush"; transcript_name "ush-RC"; protein_id "FBpp0302929"; gene_biotype "protein_coding";
"""

        pred_gtf = """2L\tStringTie\ttranscript\t476445\t479670\t1000\t-\t.\tgene_id "Stringtie.63"; transcript_id "Stringtie.63.1"; cov "141.769424"; FPKM "inf";
2L\tStringTie\texon\t476445\t478204\t1000\t-\t.\tgene_id "Stringtie.63"; transcript_id "Stringtie.63.1"; exon_number "1"; cov "149.294586";
2L\tStringTie\texon\t479407\t479670\t1000\t-\t.\tgene_id "Stringtie.63"; transcript_id "Stringtie.63.1"; exon_number "2"; cov "91.601692";"""

        ref_lines = [parsers.GTF.GtfLine(line)
                     for line in filter(lambda x: x!='', ref_gtf.split("\n"))]
        self.ref = loci.Transcript(ref_lines[0])
        self.ref.logger = self.logger
        for l in ref_lines[1:]:
            self.ref.add_exon(l)
        self.ref.finalize()
        
        pred_lines = [parsers.GTF.GtfLine(line)
                      for line in filter(lambda x: x!='', pred_gtf.split("\n"))]
        self.pred = loci.Transcript(pred_lines[0])
        for l in pred_lines[1:]:
            self.pred.add_exon(l)
        self.pred.finalize()
        
    def test_code(self):

        self.ref.finalize()
        self.assertGreater(len(self.ref.combined_cds), 0)
        self.assertEqual(len(self.ref.selected_cds_introns), 7)
        self.assertEqual(len(self.ref.combined_cds_introns), 7)


class TranscriptTesterPositive(unittest.TestCase):

    logger = create_null_logger("test_at")

    tr_gff = """Chr2    TAIR10    mRNA    626642    629176    .    +    .    ID=AT2G02380.1;Parent=AT2G02380
Chr2    TAIR10    exon    626642    626780    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    five_prime_UTR    626642    626780    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    exon    626842    626880    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    five_prime_UTR    626842    626877    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    CDS    626878    626880    .    +    0    Parent=AT2G02380.1
Chr2    TAIR10    exon    626963    627059    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    CDS    626963    627059    .    +    0    Parent=AT2G02380.1
Chr2    TAIR10    exon    627137    627193    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    CDS    627137    627193    .    +    2    Parent=AT2G02380.1
Chr2    TAIR10    exon    627312    627397    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    CDS    627312    627397    .    +    2    Parent=AT2G02380.1
Chr2    TAIR10    exon    627488    627559    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    CDS    627488    627559    .    +    0    Parent=AT2G02380.1
Chr2    TAIR10    exon    627696    627749    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    CDS    627696    627749    .    +    0    Parent=AT2G02380.1
Chr2    TAIR10    exon    627840    627915    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    CDS    627840    627915    .    +    0    Parent=AT2G02380.1
Chr2    TAIR10    exon    628044    628105    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    CDS    628044    628105    .    +    2    Parent=AT2G02380.1
Chr2    TAIR10    exon    628182    628241    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    CDS    628182    628241    .    +    0    Parent=AT2G02380.1
Chr2    TAIR10    exon    628465    628676    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    CDS    628465    628569    .    +    0    Parent=AT2G02380.1
Chr2    TAIR10    three_prime_UTR    628570    628676    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    exon    629070    629176    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    three_prime_UTR    629070    629176    .    +    .    Parent=AT2G02380.1"""

    tr_lines = tr_gff.split("\n")
    for pos, line in enumerate(tr_lines):
        tr_lines[pos] = re.sub(r"\s+", r"\t", line)
        assert len(tr_lines[pos].split("\t")) == 9, line.split("\t")

    tr_gff_lines = [parsers.GFF.GffLine(line) for line in tr_lines]

    for l in tr_gff_lines:
        assert l.header is False

    def setUp(self):
        """Basic creation test."""

        self.tr = loci.Transcript(self.tr_gff_lines[0])
        for line in self.tr_gff_lines[1:]:
            self.tr.add_exon(line)
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
        self.orf.block_starts = 0
        self.orf.has_start_codon = True
        self.orf.has_stop_codon = True
        self.orf.transcriptomic = True

    def test_basics(self):

        self.assertEqual(self.tr.chrom, "Chr2")
        self.assertEqual(self.tr.strand, "+")
        self.assertEqual(self.tr.exon_num, 12)
        self.assertEqual(self.tr.exon_num, len(self.tr.exons))
        self.assertEqual(self.tr.start, 626642)
        self.assertEqual(self.tr.end, 629176)
        exons = [(626642, 626780),
                 (626842, 626880),
                 (626963, 627059),
                 (627137, 627193),
                 (627312, 627397),
                 (627488, 627559),
                 (627696, 627749),
                 (627840, 627915),
                 (628044, 628105),
                 (628182, 628241),
                 (628465, 628676),
                 (629070, 629176)]
        self.assertEqual(self.tr.exons,
                         exons,
                         self.tr.exons)

    def test_no_exons(self):

        self.tr.finalized = False
        self.tr.exons = []
        self.tr.finalize()
        self.assertEqual(self.tr.chrom, "Chr2")
        self.assertEqual(self.tr.strand, "+")
        self.assertEqual(self.tr.exon_num, 12)
        self.assertEqual(self.tr.exon_num, len(self.tr.exons))
        self.assertEqual(self.tr.start, 626642)
        self.assertEqual(self.tr.end, 629176)
        exons = [(626642, 626780),
                 (626842, 626880),
                 (626963, 627059),
                 (627137, 627193),
                 (627312, 627397),
                 (627488, 627559),
                 (627696, 627749),
                 (627840, 627915),
                 (628044, 628105),
                 (628182, 628241),
                 (628465, 628676),
                 (629070, 629176)]
        self.assertEqual(self.tr.exons,
                         exons,
                         self.tr.exons)

    def test_cds(self):
        self.assertEqual(self.tr.combined_cds, self.tr.selected_cds)
        cds = [(626878, 626880),
               (626963, 627059),
               (627137, 627193),
               (627312, 627397),
               (627488, 627559),
               (627696, 627749),
               (627840, 627915),
               (628044, 628105),
               (628182, 628241),
               (628465, 628569)]

        self.assertEqual(self.tr.combined_cds,
                         cds,
                         self.tr.combined_cds)
        self.assertEqual(self.tr.selected_cds_start, 626878)
        self.assertEqual(self.tr.selected_cds_end, 628569)

    def test_secondary_orf(self):

        self.assertEqual(self.tr.cds_not_maximal, 0)
        self.assertEqual(self.tr.cds_not_maximal_fraction, 0)

    def test_utr(self):
        self.assertEqual(self.tr.five_utr, [(626642, 626780),
                                            (626842, 626877)])
        self.assertEqual(self.tr.three_utr, [(628570, 628676),
                                             (629070, 629176)])

    def test_introns(self):

        introns = {(626781, 626841), (626881, 626962), (627060, 627136),
                   (627194, 627311), (627398, 627487), (627560, 627695),
                   (627750, 627839), (627916, 628043), (628106, 628181),
                   (628242, 628464), (628677, 629069)}

        self.assertEqual(self.tr.introns,
                         introns,
                         self.tr.introns)

        introns = {(626881, 626962), (627060, 627136), (627194, 627311),
                   (627398, 627487), (627560, 627695), (627750, 627839),
                   (627916, 628043), (628106, 628181), (628242, 628464)}

        self.assertEqual(self.tr.combined_cds_introns,
                         introns,
                         (sorted(self.tr.combined_cds_introns),
                          sorted(introns)))

        cds_introns = {(626881, 626962), (627060, 627136), (627194, 627311),
                       (627398, 627487), (627560, 627695), (627750, 627839),
                       (627916, 628043), (628106, 628181), (628242, 628464)}

        self.assertEqual(self.tr.selected_cds_introns,
                         cds_introns,
                         self.tr.selected_cds_introns
                         )

    def test_utr_metrics(self):

        """Test for UTR exon num, start distance, etc."""

        self.assertEqual(self.tr.five_utr_num, 2)
        self.assertEqual(self.tr.three_utr_num, 2)
        self.assertEqual(self.tr.five_utr_num_complete, 1)
        self.assertEqual(self.tr.three_utr_num_complete, 1)

        self.assertEqual(self.tr.five_utr_length, 626780 + 1 - 626642 + 626877 + 1 - 626842)
        self.assertEqual(self.tr.three_utr_length, 628676 + 1 - 628570 + 629176 + 1 - 629070)

        self.assertEqual(self.tr.selected_start_distance_from_tss,
                         626780 + 1 - 626642 + 626878 - 626842,
                         self.tr.selected_end_distance_from_tes)
        self.assertEqual(self.tr.selected_start_distance_from_tss, self.tr.start_distance_from_tss)

        self.assertEqual(self.tr.selected_end_distance_from_tes,
                         628676 - 628569 + 629176 + 1 - 629070,
                         self.tr.selected_end_distance_from_tes)
        self.assertEqual(self.tr.selected_end_distance_from_tes, self.tr.end_distance_from_tes)

        self.assertEqual(self.tr.selected_end_distance_from_junction, 628676 - 628569)

    def test_strip_cds(self):

        with self.assertLogs(logger=self.logger, level="DEBUG") as log_split:
            self.tr.strip_cds()

        self.assertIn("DEBUG:{}:Stripping CDS from AT2G02380.1".format(self.logger.name), log_split.output)

        self.assertEqual(self.tr.selected_cds_length, 0)
        self.assertEqual(self.tr.three_utr, [])
        self.assertEqual(self.tr.five_utr, [])
        self.assertEqual(self.tr.selected_cds, [])
        self.assertEqual(self.tr.selected_cds_start, None)
        self.assertEqual(self.tr.selected_cds_end, None)

    def test_with_no_gff_utr(self):

        """
        Test the creation of the transcript without the UTR lines, verify that everything is still alright
        :return:
        """
        tr_gff = """Chr2    TAIR10    mRNA    626642    629176    .    +    .    ID=AT2G02380.1;Parent=AT2G02380
Chr2    TAIR10    exon    626642    626780    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    exon    626842    626880    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    CDS    626878    626880    .    +    0    Parent=AT2G02380.1
Chr2    TAIR10    exon    626963    627059    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    CDS    626963    627059    .    +    0    Parent=AT2G02380.1
Chr2    TAIR10    exon    627137    627193    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    CDS    627137    627193    .    +    2    Parent=AT2G02380.1
Chr2    TAIR10    exon    627312    627397    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    CDS    627312    627397    .    +    2    Parent=AT2G02380.1
Chr2    TAIR10    exon    627488    627559    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    CDS    627488    627559    .    +    0    Parent=AT2G02380.1
Chr2    TAIR10    exon    627696    627749    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    CDS    627696    627749    .    +    0    Parent=AT2G02380.1
Chr2    TAIR10    exon    627840    627915    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    CDS    627840    627915    .    +    0    Parent=AT2G02380.1
Chr2    TAIR10    exon    628044    628105    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    CDS    628044    628105    .    +    2    Parent=AT2G02380.1
Chr2    TAIR10    exon    628182    628241    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    CDS    628182    628241    .    +    0    Parent=AT2G02380.1
Chr2    TAIR10    exon    628465    628676    .    +    .    Parent=AT2G02380.1
Chr2    TAIR10    CDS    628465    628569    .    +    0    Parent=AT2G02380.1
Chr2    TAIR10    exon    629070    629176    .    +    .    Parent=AT2G02380.1"""

        tr_lines = tr_gff.split("\n")
        logger = create_default_logger("test")
        logger.setLevel("INFO")
        for pos, line in enumerate(tr_lines):
            tr_lines[pos] = re.sub(r"\s+", "\t", line)
            assert len(tr_lines[pos].split("\t")) == 9, line.split("\t")

        tr_gff_lines = [parsers.GFF.GffLine(line) for line in tr_lines]

        transcript = loci.Transcript(tr_gff_lines[0],
                                            logger=logger)
        for line in tr_gff_lines[1:]:
            transcript.add_exon(line)

        self.assertEqual(transcript.exons, self.tr.exons)
        self.assertNotEqual([], transcript.combined_cds)
        transcript.finalize()
        self.assertTrue(transcript.is_coding)
        self.assertEqual(transcript.five_utr, self.tr.five_utr)
        self.assertEqual(transcript.three_utr, self.tr.three_utr)

    def test_remove_utr(self):
        """Test for CDS stripping. We remove the UTRs and verify
        that start/end have moved, no UTR is present, etc."""

        self.tr.remove_utrs()
        self.assertEqual(self.tr.selected_cds_start, self.tr.start)
        self.assertEqual(self.tr.selected_cds_end, self.tr.end)
        self.assertEqual(self.tr.three_utr, [])
        self.assertEqual(self.tr.five_utr, [])
        cds = [(626878, 626880),
               (626963, 627059),
               (627137, 627193),
               (627312, 627397),
               (627488, 627559),
               (627696, 627749),
               (627840, 627915),
               (628044, 628105),
               (628182, 628241),
               (628465, 628569)]

        self.assertEqual(self.tr.combined_cds,
                         cds,
                         self.tr.combined_cds)
        self.assertEqual(self.tr.combined_utr, [], self.tr.combined_utr)

    def test_load_orf(self):

        """Test for loading a single ORF. We strip the CDS and reload it."""

        with self.assertLogs(logger=self.logger, level="DEBUG") as cm_out:
            self.tr.strip_cds()
            self.assertIn("Stripping CDS", cm_out.output[0])
        self.tr.load_orfs([self.orf])
        cds = [(626878, 626880), (626963, 627059), (627137, 627193), (627312, 627397), (627488, 627559),
               (627696, 627749), (627840, 627915), (628044, 628105), (628182, 628241), (628465, 628569)]

        self.assertEqual(self.tr.combined_cds,
                         cds,
                         self.tr.combined_cds)
        self.assertEqual(self.tr.selected_cds_start, 626878)
        self.assertEqual(self.tr.selected_cds_end, 628569)

    def test_negative_orf(self):
        """Test loading a negative strand ORF onto a multiexonic transcript.
        This should have no effect."""

        self.orf.strand = "-"
        self.tr.strip_cds()
        self.tr.load_orfs([self.orf])
        self.assertEqual(self.tr.selected_cds_start, None)

    def test_raises_invalid(self):

        self.tr.finalized = False
        self.tr.strand = None

        __current = self.tr.deepcopy()
        self.assertRaises(exceptions.InvalidTranscript,
                          self.tr.finalize)

        self.assertFalse(self.tr.finalized)
        # self.assertTrue(__current is self.tr)
        self.tr.strand = "+"
        self.tr.finalize()
        self.tr.finalized = False
        self.tr.exons += [(625878, 625880)]
        self.assertRaises(exceptions.InvalidTranscript, self.tr.finalize)

    def test_complete(self):

        self.assertTrue(self.tr.has_stop_codon)
        self.assertTrue(self.tr.has_start_codon)
        self.assertTrue(self.tr.is_complete)

    def test_lengths(self):

        self.assertEqual(self.tr.cdna_length, 1061)
        self.assertEqual(self.tr.selected_cds_length, 672)
        self.assertAlmostEqual(self.tr.combined_cds_fraction, 672 / 1061, delta=0.01)
        self.assertAlmostEqual(self.tr.selected_cds_fraction, 672 / 1061, delta=0.01)

    def testSegments(self):

        self.assertEqual(self.tr.combined_cds_num, 10)
        self.assertEqual(self.tr.selected_cds_num, 10)
        self.assertEqual(self.tr.highest_cds_exon_number, 10)
        self.assertEqual(self.tr.max_intron_length, 393)
        self.assertEqual(self.tr.number_internal_orfs, 1)

    def testDoubleOrf(self):

        """Test to verify the introduction of multiple ORFs."""

        self.tr.strip_cds()
        self.tr.finalized = False

        first_orf = parsers.bed12.BED12()
        first_orf.chrom = self.tr.id
        first_orf.start = 1
        first_orf.end = self.tr.cdna_length
        first_orf.name = "first"
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

        self.assertTrue(
            loci.Transcript.is_overlapping_cds(first_orf, second_orf))

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

        self.assertFalse(
            loci.Transcript.is_overlapping_cds(first_orf, third_orf))
        self.assertFalse(
            loci.Transcript.is_overlapping_cds(second_orf, third_orf))

        self.assertFalse(third_orf == second_orf)
        self.assertFalse(first_orf == second_orf)
        self.assertFalse(first_orf == third_orf)

        candidates = [first_orf, second_orf, third_orf]

        # self.assertEqual(len(Mikado.py.loci.transcript.Transcript.find_overlapping_cds(candidates)), 2)

        logger = create_null_logger("null")
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


class AtPositiveWrongSplit(unittest.TestCase):

    def setUp(self):

        trlines="""Chr4\tCufflinks\ttranscript\t15490563\t15495994\t1000\t+\t.\tgene_id "cufflinks_star_at.17370"; transcript_id "cufflinks_cufflinks_star_at.17370.1"; exon_number "1"; conf_lo "14.664134"; FPKM "15.0125792667"; frac "0.939064"; canonical_proportion "1.0"; conf_hi "15.635854"; cov "392.846627";
Chr4\tCufflinks\texon\t15490563\t15491286\t.\t+\t.\tgene_id "cufflinks_star_at.17370"; transcript_id "cufflinks_cufflinks_star_at.17370.1";
Chr4\tCufflinks\texon\t15492353\t15492409\t.\t+\t.\tgene_id "cufflinks_star_at.17370"; transcript_id "cufflinks_cufflinks_star_at.17370.1";
Chr4\tCufflinks\texon\t15492505\t15492567\t.\t+\t.\tgene_id "cufflinks_star_at.17370"; transcript_id "cufflinks_cufflinks_star_at.17370.1";
Chr4\tCufflinks\texon\t15492647\t15492715\t.\t+\t.\tgene_id "cufflinks_star_at.17370"; transcript_id "cufflinks_cufflinks_star_at.17370.1";
Chr4\tCufflinks\texon\t15492819\t15494669\t.\t+\t.\tgene_id "cufflinks_star_at.17370"; transcript_id "cufflinks_cufflinks_star_at.17370.1";
Chr4\tCufflinks\texon\t15495054\t15495174\t.\t+\t.\tgene_id "cufflinks_star_at.17370"; transcript_id "cufflinks_cufflinks_star_at.17370.1";
Chr4\tCufflinks\texon\t15495267\t15495466\t.\t+\t.\tgene_id "cufflinks_star_at.17370"; transcript_id "cufflinks_cufflinks_star_at.17370.1";
Chr4\tCufflinks\texon\t15495556\t15495687\t.\t+\t.\tgene_id "cufflinks_star_at.17370"; transcript_id "cufflinks_cufflinks_star_at.17370.1";
Chr4\tCufflinks\texon\t15495769\t15495908\t.\t+\t.\tgene_id "cufflinks_star_at.17370"; transcript_id "cufflinks_cufflinks_star_at.17370.1";
Chr4\tCufflinks\texon\t15495994\t15495994\t.\t+\t.\tgene_id "cufflinks_star_at.17370"; transcript_id "cufflinks_cufflinks_star_at.17370.1";"""

        trlines = [parsers.GTF.GtfLine(_) for _ in trlines.split("\n")]
        self.tr = loci.Transcript(trlines[0])
        [self.tr.add_exon(_) for _ in trlines[1:]]
        self.tr.finalize()

        self.bed1 = parsers.bed12.BED12()
        self.bed1.header = False
        self.bed1.chrom = self.tr.id
        self.bed1.start = 1
        self.bed1.end = 3358
        self.bed1.name = "g.216019_cufflinks_cufflinks_star_at.17370.1|m.216019_type:complete_len:254_(+)"
        self.assertEqual(self.bed1.end, self.tr.cdna_length, (self.bed1.end, self.tr.cdna_length))
        self.bed1.strand = "+"
        self.bed1.thick_start = 341
        self.bed1.thick_end = 1102
        self.bed1.score = 0
        self.bed1.has_stop_codon = True
        self.bed1.has_start_codon = True
        self.bed1.transcriptomic = True
        self.assertFalse(self.bed1.invalid)
        self.assertEqual(self.bed1.cds_len, 762)

        self.bed2 = parsers.bed12.BED12()
        self.bed2.header = False
        self.bed2.chrom = self.tr.id
        self.bed2.start = 1
        self.bed2.end = 3358
        self.bed2.name = "g.216018_cufflinks_cufflinks_star_at.17370.1|m.216018_type:3prime_partial_len:380_(+)"
        self.assertEqual(self.bed2.end, self.tr.cdna_length, (self.bed2.end, self.tr.cdna_length))
        self.bed2.strand = "+"
        self.bed2.thick_start = 2222
        self.bed2.thick_end = 3358
        self.bed2.score = 0
        self.bed2.has_stop_codon = False
        self.bed2.has_start_codon = True
        self.bed2.transcriptomic = True
        self.assertFalse(self.bed2.invalid)
        self.assertEqual(self.bed2.cds_len, 1137)

    def test_split(self):

        self.tr.load_orfs([self.bed1, self.bed2])

        self.assertEqual(self.tr.selected_cds_start, 15494127)
        self.assertEqual(self.tr.selected_cds_end, 15495994)

        self.assertEqual(self.tr.combined_cds_start, 15490903)
        # The other CDS starts at 15494127

        logger = create_default_logger(self.tr.id)
        logger.setLevel("WARN")
        self.tr.logger = logger

        new_transcripts = [_ for _ in self.tr.split_by_cds()]
        self.assertEqual(len(new_transcripts), 2)


class AugustusTester(unittest.TestCase):

    logger = create_default_logger("augustus")
    logger.setLevel("DEBUG")

    def test_truncated(self):

        lines = """Triticum_aestivum_CS42_TGACv1_scaffold_000043_1AL	Triticum_aestivum_CS42_TGACv1_TRIAE4565_Augustus	mRNA	1	2785	.	+	.	ID=TRIAE4565_1AL_Aug_0021880.1;Parent=TRIAE4565_1AL_Aug_0021880;Name=TRIAE4565_1AL_Aug_0021880.1
Triticum_aestivum_CS42_TGACv1_scaffold_000043_1AL	Triticum_aestivum_CS42_TGACv1_TRIAE4565_Augustus	CDS	1601	2446	.	+	1	ID=TRIAE4565_1AL_Aug_0021880.1.CDS1;Parent=TRIAE4565_1AL_Aug_0021880.1
Triticum_aestivum_CS42_TGACv1_scaffold_000043_1AL	Triticum_aestivum_CS42_TGACv1_TRIAE4565_Augustus	exon	1601	2446	.	+	.	ID=TRIAE4565_1AL_Aug_0021880.1.exon1;Parent=TRIAE4565_1AL_Aug_0021880.1
Triticum_aestivum_CS42_TGACv1_scaffold_000043_1AL	Triticum_aestivum_CS42_TGACv1_TRIAE4565_Augustus	CDS	2540	2654	.	+	1	ID=TRIAE4565_1AL_Aug_0021880.1.CDS2;Parent=TRIAE4565_1AL_Aug_0021880.1
Triticum_aestivum_CS42_TGACv1_scaffold_000043_1AL	Triticum_aestivum_CS42_TGACv1_TRIAE4565_Augustus	exon	2540	2785	.	+	.	ID=TRIAE4565_1AL_Aug_0021880.1.exon2;Parent=TRIAE4565_1AL_Aug_0021880.1
Triticum_aestivum_CS42_TGACv1_scaffold_000043_1AL	Triticum_aestivum_CS42_TGACv1_TRIAE4565_Augustus	three_prime_UTR	2655	2785	.	+	.	ID=TRIAE4565_1AL_Aug_0021880.1.three_prime_UTR1;Parent=TRIAE4565_1AL_Aug_0021880.1"""

        lines = [parsers.GFF.GffLine("\t".join(_.split())) for _ in lines.split("\n")]

        transcript = loci.Transcript(lines[0], logger=self.logger)
        transcript.add_exons(lines[1:])

        with self.assertLogs("augustus", level="WARNING") as cm_out:
            transcript.finalize()
            self.assertTrue(any(
                "The transcript TRIAE4565_1AL_Aug_0021880.1 has coordinates 1:2785" in _ for
            _ in cm_out.output))

        self.assertTrue(transcript.is_coding)

    def test_three_truncated(self):
        lines = """Triticum_aestivum_CS42_TGACv1_scaffold_000112_1AL	Triticum_aestivum_CS42_TGACv1_TRIAE4565_Augustus	mRNA	204336	224434	.	+	.	ID=TRIAE4565_1AL_Aug_0024630.1;Parent=TRIAE4565_1AL_Aug_0024630;Name=TRIAE4565_1AL_Aug_0024630.1
Triticum_aestivum_CS42_TGACv1_scaffold_000112_1AL	Triticum_aestivum_CS42_TGACv1_TRIAE4565_Augustus	exon	204336	205303	.	+	.	ID=TRIAE4565_1AL_Aug_0024630.1.exon1;Parent=TRIAE4565_1AL_Aug_0024630.1
Triticum_aestivum_CS42_TGACv1_scaffold_000112_1AL	Triticum_aestivum_CS42_TGACv1_TRIAE4565_Augustus	five_prime_UTR	204336	204546	.	+	.	ID=TRIAE4565_1AL_Aug_0024630.1.five_prime_UTR1;Parent=TRIAE4565_1AL_Aug_0024630.1
Triticum_aestivum_CS42_TGACv1_scaffold_000112_1AL	Triticum_aestivum_CS42_TGACv1_TRIAE4565_Augustus	CDS	204547	205303	.	+	0	ID=TRIAE4565_1AL_Aug_0024630.1.CDS1;Parent=TRIAE4565_1AL_Aug_0024630.1
Triticum_aestivum_CS42_TGACv1_scaffold_000112_1AL	Triticum_aestivum_CS42_TGACv1_TRIAE4565_Augustus	CDS	206227	207040	.	+	2	ID=TRIAE4565_1AL_Aug_0024630.1.CDS2;Parent=TRIAE4565_1AL_Aug_0024630.1
Triticum_aestivum_CS42_TGACv1_scaffold_000112_1AL	Triticum_aestivum_CS42_TGACv1_TRIAE4565_Augustus	exon	206227	207040	.	+	.	ID=TRIAE4565_1AL_Aug_0024630.1.exon2;Parent=TRIAE4565_1AL_Aug_0024630.1"""

        lines = [parsers.GFF.GffLine("\t".join(_.split())) for _ in lines.split("\n")]

        transcript = loci.Transcript(lines[0], logger=self.logger)
        transcript.add_exons(lines[1:])

        with self.assertLogs("augustus", level="WARNING") as cm_out:
            transcript.finalize()
            self.assertTrue(any(
                "The transcript TRIAE4565_1AL_Aug_0024630.1 has coordinates 204336:224434" in _ for
            _ in cm_out.output))

        self.assertTrue(transcript.is_coding)

    def test_invalid_three_truncated(self):
        lines = """Triticum_aestivum_CS42_TGACv1_scaffold_000112_1AL	Triticum_aestivum_CS42_TGACv1_TRIAE4565_Augustus	mRNA	204336	225434	.	+	.	ID=TRIAE4565_1AL_Aug_0024630.1;Parent=TRIAE4565_1AL_Aug_0024630;Name=TRIAE4565_1AL_Aug_0024630.1
Triticum_aestivum_CS42_TGACv1_scaffold_000112_1AL	Triticum_aestivum_CS42_TGACv1_TRIAE4565_Augustus	exon	204336	205303	.	+	.	ID=TRIAE4565_1AL_Aug_0024630.1.exon1;Parent=TRIAE4565_1AL_Aug_0024630.1
Triticum_aestivum_CS42_TGACv1_scaffold_000112_1AL	Triticum_aestivum_CS42_TGACv1_TRIAE4565_Augustus	five_prime_UTR	204336	204546	.	+	.	ID=TRIAE4565_1AL_Aug_0024630.1.five_prime_UTR1;Parent=TRIAE4565_1AL_Aug_0024630.1
Triticum_aestivum_CS42_TGACv1_scaffold_000112_1AL	Triticum_aestivum_CS42_TGACv1_TRIAE4565_Augustus	CDS	204547	205303	.	+	0	ID=TRIAE4565_1AL_Aug_0024630.1.CDS1;Parent=TRIAE4565_1AL_Aug_0024630.1
Triticum_aestivum_CS42_TGACv1_scaffold_000112_1AL	Triticum_aestivum_CS42_TGACv1_TRIAE4565_Augustus	CDS	206227	207042	.	+	2	ID=TRIAE4565_1AL_Aug_0024630.1.CDS2;Parent=TRIAE4565_1AL_Aug_0024630.1
Triticum_aestivum_CS42_TGACv1_scaffold_000112_1AL	Triticum_aestivum_CS42_TGACv1_TRIAE4565_Augustus	exon	206227	207042	.	+	.	ID=TRIAE4565_1AL_Aug_0024630.1.exon2;Parent=TRIAE4565_1AL_Aug_0024630.1
Triticum_aestivum_CS42_TGACv1_scaffold_000112_1AL	Triticum_aestivum_CS42_TGACv1_TRIAE4565_Augustus	exon	208227	210040	.	+	.	ID=TRIAE4565_1AL_Aug_0024630.1.exon2;Parent=TRIAE4565_1AL_Aug_0024630.1"""

        lines = [parsers.GFF.GffLine("\t".join(_.split())) for _ in lines.split("\n")]

        transcript = loci.Transcript(lines[0], logger=self.logger)
        transcript.add_exons(lines[1:])

        with self.assertLogs("augustus", level="WARNING") as cm_out:
            transcript.finalize()
            self.assertTrue(any(
                "The transcript TRIAE4565_1AL_Aug_0024630.1 has coordinates 204336:225434" in _ for
            _ in cm_out.output))
            # self.assertTrue(any(
            #     "strip_cds" in _ for
            # _ in cm_out.output))

        self.assertFalse(transcript.is_coding)

    def test_valid_three_truncated(self):
        
        """
        Picked from the EnsEMBL Human GTF (v. 70)
        :return: 
        """

        lines = """11\tnonsense_mediated_decay\texon\t134177086\t134177102\t.\t+\t.\tgene_id "ENSG00000166105"; transcript_id "ENST00000455971"; exon_number "1"; gene_name "GLB1L3"; gene_biotype "protein_coding"; transcript_name "GLB1L3-006"; exon_id "ENSE00002461794";
11\tnonsense_mediated_decay\tCDS\t134177086\t134177102\t.\t+\t2\tgene_id "ENSG00000166105"; transcript_id "ENST00000455971"; exon_number "1"; gene_name "GLB1L3"; gene_biotype "protein_coding"; transcript_name "GLB1L3-006"; protein_id "ENSP00000397929";
11\tnonsense_mediated_decay\texon\t134179522\t134179657\t.\t+\t.\tgene_id "ENSG00000166105"; transcript_id "ENST00000455971"; exon_number "2"; gene_name "GLB1L3"; gene_biotype "protein_coding"; transcript_name "GLB1L3-006"; exon_id "ENSE00002147723";
11\tnonsense_mediated_decay\tCDS\t134179522\t134179657\t.\t+\t0\tgene_id "ENSG00000166105"; transcript_id "ENST00000455971"; exon_number "2"; gene_name "GLB1L3"; gene_biotype "protein_coding"; transcript_name "GLB1L3-006"; protein_id "ENSP00000397929";
11\tnonsense_mediated_decay\texon\t134180465\t134180545\t.\t+\t.\tgene_id "ENSG00000166105"; transcript_id "ENST00000455971"; exon_number "3"; gene_name "GLB1L3"; gene_biotype "protein_coding"; transcript_name "GLB1L3-006"; exon_id "ENSE00001278318";
11\tnonsense_mediated_decay\tCDS\t134180465\t134180545\t.\t+\t2\tgene_id "ENSG00000166105"; transcript_id "ENST00000455971"; exon_number "3"; gene_name "GLB1L3"; gene_biotype "protein_coding"; transcript_name "GLB1L3-006"; protein_id "ENSP00000397929";
11\tnonsense_mediated_decay\texon\t134180958\t134181064\t.\t+\t.\tgene_id "ENSG00000166105"; transcript_id "ENST00000455971"; exon_number "4"; gene_name "GLB1L3"; gene_biotype "protein_coding"; transcript_name "GLB1L3-006"; exon_id "ENSE00001140726";
11\tnonsense_mediated_decay\tCDS\t134180958\t134181064\t.\t+\t2\tgene_id "ENSG00000166105"; transcript_id "ENST00000455971"; exon_number "4"; gene_name "GLB1L3"; gene_biotype "protein_coding"; transcript_name "GLB1L3-006"; protein_id "ENSP00000397929";
11\tnonsense_mediated_decay\texon\t134182243\t134182383\t.\t+\t.\tgene_id "ENSG00000166105"; transcript_id "ENST00000455971"; exon_number "5"; gene_name "GLB1L3"; gene_biotype "protein_coding"; transcript_name "GLB1L3-006"; exon_id "ENSE00003177017";
11\tnonsense_mediated_decay\tCDS\t134182243\t134182383\t.\t+\t0\tgene_id "ENSG00000166105"; transcript_id "ENST00000455971"; exon_number "5"; gene_name "GLB1L3"; gene_biotype "protein_coding"; transcript_name "GLB1L3-006"; protein_id "ENSP00000397929";
11\tnonsense_mediated_decay\texon\t134182710\t134182781\t.\t+\t.\tgene_id "ENSG00000166105"; transcript_id "ENST00000455971"; exon_number "6"; gene_name "GLB1L3"; gene_biotype "protein_coding"; transcript_name "GLB1L3-006"; exon_id "ENSE00003096760";
11\tnonsense_mediated_decay\tCDS\t134182710\t134182781\t.\t+\t0\tgene_id "ENSG00000166105"; transcript_id "ENST00000455971"; exon_number "6"; gene_name "GLB1L3"; gene_biotype "protein_coding"; transcript_name "GLB1L3-006"; protein_id "ENSP00000397929";
11\tnonsense_mediated_decay\texon\t134183835\t134183922\t.\t+\t.\tgene_id "ENSG00000166105"; transcript_id "ENST00000455971"; exon_number "7"; gene_name "GLB1L3"; gene_biotype "protein_coding"; transcript_name "GLB1L3-006"; exon_id "ENSE00003040614";
11\tnonsense_mediated_decay\tCDS\t134183835\t134183837\t.\t+\t0\tgene_id "ENSG00000166105"; transcript_id "ENST00000455971"; exon_number "7"; gene_name "GLB1L3"; gene_biotype "protein_coding"; transcript_name "GLB1L3-006"; protein_id "ENSP00000397929";
11\tnonsense_mediated_decay\tstop_codon\t134183838\t134183840\t.\t+\t0\tgene_id "ENSG00000166105"; transcript_id "ENST00000455971"; exon_number "7"; gene_name "GLB1L3"; gene_biotype "protein_coding"; transcript_name "GLB1L3-006";
11\tnonsense_mediated_decay\texon\t134184224\t134184335\t.\t+\t.\tgene_id "ENSG00000166105"; transcript_id "ENST00000455971"; exon_number "8"; gene_name "GLB1L3"; gene_biotype "protein_coding"; transcript_name "GLB1L3-006"; exon_id "ENSE00003002659";
11\tnonsense_mediated_decay\texon\t134188525\t134188641\t.\t+\t.\tgene_id "ENSG00000166105"; transcript_id "ENST00000455971"; exon_number "9"; gene_name "GLB1L3"; gene_biotype "protein_coding"; transcript_name "GLB1L3-006"; exon_id "ENSE00003191545";
11\tnonsense_mediated_decay\texon\t134188771\t134189178\t.\t+\t.\tgene_id "ENSG00000166105"; transcript_id "ENST00000455971"; exon_number "10"; gene_name "GLB1L3"; gene_biotype "protein_coding"; transcript_name "GLB1L3-006"; exon_id "ENSE00001441085";"""

        lines = [parsers.GTF.GtfLine("\t".join(_.split("\t"))) for _ in lines.split("\n")]
        assert all([line.header is False for line in lines])

        transcript = loci.Transcript(lines[0], logger=self.logger)

        transcript.add_exons(lines[1:])
        with self.assertLogs("augustus", level="DEBUG") as cm_out:
            transcript.finalize()

        self.assertTrue(transcript.is_coding)
        self.assertEqual(560, transcript.selected_cds_length,
                         sum(
                             _[1][1] - _[1][0] + 1 for _ in transcript.selected_internal_orf if _[0] == "CDS")
                         )

if __name__ == '__main__':
    unittest.main()
