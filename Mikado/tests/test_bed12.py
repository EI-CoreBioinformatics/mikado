from Mikado.parsers.bed12 import BED12
# from Mikado.transcripts import Transcript
import unittest
from Bio.Seq import Seq

class Bed12GenToTrans(unittest.TestCase):

    def setUp(self):
        pass

    def test_positive_mono_transfer(self):

        string_bed = "1\t10\t500\ttest\t0\t+\t300\t390\t0\t1\t490\t0"


        bed = BED12(string_bed)
        self.assertFalse(bed.invalid)
        self.assertFalse(bed.header)
        self.assertEqual(bed.thick_start, 301)
        self.assertEqual(bed.thick_end, 390)
        self.assertEqual(bed.blocks, [(11, 500)])

        seq = "AAA" * 96 + "AA"
        seq += "ATG" * 29
        seq += "TGA"
        seq += "A" * 110
        self.assertEqual(len(seq), len(bed))

        tbed = bed.to_transcriptomic(sequence=seq)
        self.assertEqual(tbed.start, 1)
        self.assertEqual(tbed.end, 500 - 11 + 1)
        self.assertEqual(tbed.thick_end - tbed.thick_start + 1, 90)
        self.assertEqual(tbed.thick_start, 291)
        self.assertEqual(tbed.thick_end, 380)
        self.assertTrue(tbed.has_start_codon)
        self.assertTrue(tbed.has_stop_codon)
        self.assertEqual(seq[tbed.thick_start - 1:tbed.thick_end], "ATG" * 29 + "TGA")

    def test_negative_mono_transfer(self):

        string_bed = "1\t10\t500\ttest\t0\t-\t300\t390\t0\t1\t490\t0"

        bed = BED12(string_bed)
        self.assertFalse(bed.invalid)
        self.assertFalse(bed.header)
        self.assertEqual(bed.thick_start, 301)
        self.assertEqual(bed.thick_end, 390)
        self.assertEqual(bed.blocks, [(11, 500)])

        seq = "A" * 110
        seq += "ATG" * 29
        seq += "TGA"
        seq += "A" * (490 - 110 - 90)
        self.assertEqual(len(seq), len(bed))

        tbed = bed.to_transcriptomic(sequence=seq)
        self.assertEqual(tbed.start, 1)
        self.assertEqual(tbed.end, 500 - 11 + 1)

        self.assertEqual(tbed.thick_end - tbed.thick_start + 1, 90)
        self.assertEqual(tbed.thick_start, 111)
        self.assertEqual(tbed.thick_end, 200)
        self.assertTrue(tbed.has_start_codon)
        self.assertTrue(tbed.has_stop_codon)

        self.assertEqual(seq[tbed.thick_start - 1:tbed.thick_end], "ATG" * 29 + "TGA")

    def test_diexonic_pos_transfer(self):

        string_bed = "1\t10\t1000\ttest\t0\t+\t80\t920\t0\t2\t190,200\t0,790"
        bed = BED12(string_bed)
        self.assertFalse(bed.invalid or bed.header)
        self.assertEqual(bed.start, 11)
        self.assertEqual(bed.end, 1000)
        self.assertEqual(bed.blocks, [(11, 200), (801, 1000)])
        self.assertEqual(bed.thick_start, 81)
        self.assertEqual(bed.thick_end, 920)

        seq = "A" * 70
        seq += "ATG" * 79
        seq += "TAA" * 1
        seq += "A" * 80
        self.assertEqual(len(seq), sum(bed.block_sizes))
        tbed = bed.to_transcriptomic(sequence=seq)
        self.assertEqual(tbed.start, 1)
        self.assertEqual(tbed.end, 390)
        self.assertEqual(tbed.thick_end - tbed.thick_start + 1, 240)
        self.assertEqual(tbed.thick_start, 71)
        self.assertEqual(tbed.thick_end, 310)
        self.assertTrue(tbed.has_start_codon)
        self.assertTrue(tbed.has_stop_codon)

        self.assertEqual(seq[tbed.thick_start - 1:tbed.thick_end], "ATG" * 79 + "TAA")

    def test_diexonic_neg_transfer(self):

        string_bed = "1\t10\t1000\ttest\t0\t-\t80\t920\t0\t2\t190,200\t0,790"
        bed = BED12(string_bed)
        self.assertFalse(bed.invalid or bed.header)
        self.assertEqual(bed.start, 11)
        self.assertEqual(bed.end, 1000)
        self.assertEqual(bed.blocks, [(11, 200), (801, 1000)])
        self.assertEqual(bed.thick_start, 81)
        self.assertEqual(bed.thick_end, 920)

        seq = "A" * 80
        seq += "ATG" * 79
        seq += "TAA" * 1
        seq += "A" * 70
        self.assertEqual(len(seq), sum(bed.block_sizes))
        tbed = bed.to_transcriptomic(sequence=seq)
        self.assertEqual(tbed.start, 1)
        self.assertEqual(tbed.end, 390)
        self.assertEqual(tbed.thick_end - tbed.thick_start + 1, 240)
        self.assertEqual(tbed.thick_start, 81)
        self.assertEqual(tbed.thick_end, 320)
        self.assertTrue(tbed.has_start_codon)
        self.assertTrue(tbed.has_stop_codon)

        self.assertEqual(seq[tbed.thick_start - 1:tbed.thick_end], "ATG" * 79 + "TAA")

    def test_wheat_1(self):

        string_bed = "chr7A\t207087445\t207089574\tTraesCS7A01G235400.1\t0\t-\t207087615\t207088433\t0\t3\t457,393,30\t0,603,2099"
        string_seq = """CGCGTCGGTGCATCCGGATACGTCGCCTGGGCTACACAATGGCGCTGATCGATTGGATAG
AACTGAGTGATGATGCAGAGATTATTGAATTGAGCAGTAGCGAGGAGAATGTCGAAGAAT
ATCAGGGAGCTTCCACACAGAGCAACTCAGCTCAGCACCAAGCTACATTGCATGATGATC
AGACCATGTTTGTCACAGAAGGTGAACGAAGGGAAGAGGCTACCGAATCTGGAAATGCAG
AGGAAGCTACTGCATCCTCCTCCGTTACTGAAAAAGAGTCTGGAAATCCAGAGGAAGCTA
CTGTATCAGAGACTTGTCCTCACACTCCAACAGCCGTTGTTACTCTAAAAGAGTCCAGAA
ATCCAGATGAAGCTACTGTATCACAAACTTGTCCTCACACTACAACAGCTGTGACGTTCC
CCGGACCTGATCCCTTCTCTCCGAAGGCTCTGACTTTTCAAGCCTGTGTCGCCAAGCCGC
GGAGAGCGAAGGCGAAACGTCGCAGGAAGCTCTTCCACGCCGATACACCTTGGACTTGGA
GGAGCCCTAGGCTTGAAGACAAGCATAAAGGTCGTAGAAGGCCCGTGGAAGAGCTGGCAG
TTGAGCGCAAGGCTTCTGCCATGCCGGAAGAAGAAGAACAGAAGAACACCTCGGCTAGGA
GCCGGAGGAAATGCAAGGTGTTCGGCAAACGGTGTAACCTCGGCAGATAGCCTCCACTGA
TCCCGCAAAACTAAGGTAGTCGACTGAAGATGAAAGATGATCGTCTCGAAACTGTCTGGA
TTTTGGCGGTCCAGTTTGCCTCTCCTGTGACGCGGGAAGAATCACTAGTGTTGGTGGGGT
TTGGTGGGGTCGTGTTGCCCTTGCTGCGCTGTCCATATGA"""

        string_seq = "".join(string_seq.split("\n"))

        string_cds = """ATGGCGCTGATCGATTGGATAGAACTGAGTGATGATGCAGAGATTATTGAATTGAGCAGT
AGCGAGGAGAATGTCGAAGAATATCAGGGAGCTTCCACACAGAGCAACTCAGCTCAGCAC
CAAGCTACATTGCATGATGATCAGACCATGTTTGTCACAGAAGGTGAACGAAGGGAAGAG
GCTACCGAATCTGGAAATGCAGAGGAAGCTACTGCATCCTCCTCCGTTACTGAAAAAGAG
TCTGGAAATCCAGAGGAAGCTACTGTATCAGAGACTTGTCCTCACACTCCAACAGCCGTT
GTTACTCTAAAAGAGTCCAGAAATCCAGATGAAGCTACTGTATCACAAACTTGTCCTCAC
ACTACAACAGCTGTGACGTTCCCCGGACCTGATCCCTTCTCTCCGAAGGCTCTGACTTTT
CAAGCCTGTGTCGCCAAGCCGCGGAGAGCGAAGGCGAAACGTCGCAGGAAGCTCTTCCAC
GCCGATACACCTTGGACTTGGAGGAGCCCTAGGCTTGAAGACAAGCATAAAGGTCGTAGA
AGGCCCGTGGAAGAGCTGGCAGTTGAGCGCAAGGCTTCTGCCATGCCGGAAGAAGAAGAA
CAGAAGAACACCTCGGCTAGGAGCCGGAGGAAATGCAAGGTGTTCGGCAAACGGTGTAAC
CTCGGCAGATAG"""
        string_cds = "".join(string_cds.split("\n"))

        bed = BED12(string_bed)
        self.assertFalse(bed.invalid or bed.header)
        self.assertEqual(bed.start, 207087445+1)
        self.assertEqual(bed.end, 207089574)
        self.assertEqual(bed.strand, "-")
        self.assertEqual(bed.blocks, [(207087446,207087902), (207088049,207088441), (207089545,207089574) ])
        self.assertEqual(bed.thick_start, 207087616)
        self.assertEqual(bed.thick_end, 207088433)

        self.assertEqual(len(string_seq), sum(bed.block_sizes))
        tbed = bed.to_transcriptomic(sequence=string_seq)
        self.assertEqual(tbed.thick_end - tbed.thick_start + 1, 672)
        self.assertEqual(tbed.thick_start, string_seq.index("ATGGCGCTGATCGATTGGA") + 1)
        self.assertEqual(tbed.thick_start, 39)
        self.assertEqual(tbed.thick_end, string_seq.index("CTCGGCAGATAG") + len("CTCGGCAGATAG"))
        self.assertEqual(string_cds, string_seq[tbed.thick_start - 1:tbed.thick_end])

        print(tbed)


if __name__ == "__main__":
    unittest.main()