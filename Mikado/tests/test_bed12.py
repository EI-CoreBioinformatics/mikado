from ..parsers.bed12 import BED12, _translate_str, standard, get_tables
from ..transcripts import Transcript
import unittest
from Bio.Data import CodonTable


class TestTranslate(unittest.TestCase):

    """
    >>>
    >>> table = CodonTable.ambiguous_dna_by_id[1]
    >>> _translate_str("AAA", table)
    'K'
    >>> _translate_str("TAR", table)
    '*'
    >>> _translate_str("TAN", table)
    'X'
    >>> _translate_str("TAN", table, pos_stop="@")
    '@'
    >>> _translate_str("TA?", table)

    >>> _translate_str("ATGCCCTAG", table, cds=True)
    'MP'
    >>> _translate_str("AAACCCTAG", table, cds=True)
    Traceback (most recent call last):
       ...
    Bio.Data.CodonTable.TranslationError: First codon 'AAA' is not a start codon
    >>> _translate_str("ATGCCCTAGCCCTAG", table, cds=True)
    """

    def test_from_Bio(self):

        self.assertEqual(_translate_str("TAN", standard, pos_stop="X"), "X")
        for codon, amino in standard.forward_table.forward_table.items():
            self.assertEqual(_translate_str(codon, standard), amino)

        self.assertEqual(_translate_str("TAR", standard), "*")
        self.assertEqual(_translate_str("TAR", standard, stop_symbol="?"), "?")

        self.assertEqual(_translate_str("TAN", standard, pos_stop="U"), "U")
        self.assertEqual(_translate_str("TAN", standard, pos_stop=b"U"), "U")
        with self.assertRaises(CodonTable.TranslationError):
            _translate_str("TA?", standard)

        ncbi_standard = CodonTable.ambiguous_dna_by_id[1]
        self.assertNotEqual(ncbi_standard.start_codons, standard.start_codons)
        self.assertEqual(ncbi_standard.stop_codons, standard.stop_codons)
        self.assertEqual(_translate_str("ATGCCCTAG", standard, cds=True), "MP*")
        self.assertEqual(_translate_str("ATGCCCTAG", standard, to_stop=True), "MP")
        with self.assertRaises(CodonTable.TranslationError):
            _translate_str("AAACCCTAG", standard, cds=True)
        with self.assertRaises(CodonTable.TranslationError) as exc:
            _translate_str("ATGCCCTAGCCCTAG", standard, cds=True)
        self.assertTrue(str(exc.exception).startswith("Extra in-frame stop codon found."),
                        str(exc.exception))
        with self.assertRaises(CodonTable.TranslationError) as exc:
            _translate_str("ATGCCCTAGCCCTAT", standard, cds=True)
        self.assertTrue(str(exc.exception).startswith("Extra in-frame stop codon. Sequence:"),
                        str(exc.exception))
        for invalid in (10, "AB", b"NT"):
            with self.assertRaises(ValueError):
                _translate_str("ATGCCCTAG", standard, cds=True, gap=invalid)

        with self.assertRaises(CodonTable.TranslationError):
            _translate_str("ATGC?G", standard)

    def test_invalid_pos_stop(self):
        for invalid in (10, "AB", b"NT"):
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                _translate_str("ATGCCCTAG", standard, cds=True, pos_stop=invalid)

    def test_ncbi_standard(self):
        standard = CodonTable.ambiguous_dna_by_id[1]
        for codon, amino in standard.forward_table.forward_table.items():
            self.assertEqual(_translate_str(codon, standard), amino)

        self.assertEqual(_translate_str("TAR", standard), "*")
        self.assertEqual(_translate_str("TAR", standard, stop_symbol="?"), "?")

        self.assertEqual(_translate_str("TAN", standard, pos_stop="U"), "U")
        with self.assertRaises(CodonTable.TranslationError):
            _translate_str("TA?", standard)

        self.assertIn("CTG", standard.start_codons)
        self.assertEqual(_translate_str("ATGCCCTAG", standard, cds=True), "MP*")
        self.assertEqual(_translate_str("ATGCCCTAG", standard, to_stop=True), "MP")
        self.assertEqual(_translate_str("CTGCCCTAG", standard, cds=True), "MP*")
        self.assertEqual(_translate_str("CTGCCCTAG", standard, cds=True, to_stop=True), "MP")

    def test_ambiguous(self):
        ambigouous = None
        for key, table in CodonTable.ambiguous_dna_by_id.items():
            amb = 0 < len([c for c in table._codon_table.forward_table.keys() if c in table.stop_codons])
            if amb:
                amb = table
                break
        if ambigouous is None:
            return  # Nothing to test

    def test_set_table(self):
        b = BED12()
        for invalid in (True, list(), "Inexistent", b"Standard2"):
            with self.assertRaises(ValueError):
                b.table = invalid
        self.assertNotIn(0, CodonTable.ambiguous_dna_by_id.keys())
        for num in range(0, max(CodonTable.ambiguous_dna_by_id.keys()) + 10):
            if num in CodonTable.ambiguous_dna_by_id.keys():
                b.table = num
                self.assertEqual(b.table, CodonTable.ambiguous_dna_by_id[num])
            elif num == 0:
                b.table = num
                self.assertEqual(b.table, standard)
            else:
                with self.assertRaises(ValueError):
                    b.table = num

        for valid in list(CodonTable.ambiguous_dna_by_name.keys()):
            b.table = valid
            self.assertEqual(b.table, CodonTable.ambiguous_dna_by_name[valid])


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
        self.assertEqual(len(seq), bed.block_sizes.sum())
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
        self.assertEqual(len(seq), bed.block_sizes.sum())
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

        self.assertEqual(len(string_seq), bed.block_sizes.sum())
        tbed = bed.to_transcriptomic(sequence=string_seq)
        self.assertEqual(tbed.thick_end - tbed.thick_start + 1, 672)
        self.assertEqual(tbed.thick_start, string_seq.index("ATGGCGCTGATCGATTGGA") + 1)
        self.assertEqual(tbed.thick_start, 39)
        self.assertEqual(tbed.thick_end, string_seq.index("CTCGGCAGATAG") + len("CTCGGCAGATAG"))
        self.assertEqual(string_cds, string_seq[tbed.thick_start - 1:tbed.thick_end])

    def test_mono_pos_bed_with_phase(self):

        string = "1\t10\t101\tID=test;phase=2;coding=True\t0\t+\t10\t101\t0\t1\t91\t0"

        seq = "A" + "CGG" * 29 + "TAA"

        bed = BED12(string, transcriptomic=False)
        self.assertFalse(bed.header)
        self.assertEqual(bed.thick_start, 11)
        self.assertFalse(bed.invalid)
        self.assertEqual(bed.name, "test")
        self.assertEqual(bed.phase, 2)

        tbed = bed.to_transcriptomic(sequence=seq)

        self.assertTrue(tbed.has_stop_codon)
        self.assertFalse(tbed.has_start_codon)
        self.assertEqual(tbed.thick_start, 1)
        self.assertEqual(tbed.thick_end, 91)
        self.assertFalse(tbed.invalid)
        self.assertEqual(tbed.phase, 2)
        self.assertTrue(tbed.transcriptomic)

    def test_mono_neg_bed_with_phase(self):

        string = "1\t10\t101\tID=test;phase=2;coding=True\t0\t-\t10\t101\t0\t1\t91\t0"

        seq = "A" + "CGG" * 29 + "TAA"

        bed = BED12(string, transcriptomic=False)
        self.assertFalse(bed.header)
        self.assertEqual(bed.thick_start, 11)
        self.assertFalse(bed.invalid)
        self.assertEqual(bed.name, "test")
        self.assertEqual(bed.phase, 2)

        tbed = bed.to_transcriptomic(sequence=seq)

        self.assertTrue(tbed.has_stop_codon)
        self.assertFalse(tbed.has_start_codon)
        self.assertEqual(tbed.thick_start, 1)
        self.assertEqual(tbed.thick_end, 91)
        self.assertFalse(tbed.invalid)
        self.assertEqual(tbed.phase, 2)
        self.assertTrue(tbed.transcriptomic)

    def test_diex_pos_bed_with_phase_one(self):

        string = "1\t10\t111\tID=test;phase=1;coding=True\t0\t+\t10\t101\t0\t1\t101\t0"

        seq = "A" + "CGG" * 29 + "TAA" + "A" * 10

        bed = BED12(string, transcriptomic=False)
        self.assertFalse(bed.header)
        self.assertEqual(bed.thick_start, 11)
        self.assertFalse(bed.invalid)
        self.assertEqual(bed.name, "test")
        self.assertEqual(bed.phase, 1)

        tbed = bed.to_transcriptomic(sequence=seq)

        self.assertFalse(tbed.has_start_codon)
        self.assertEqual(tbed.thick_start, 1)
        self.assertEqual(tbed.thick_end, 91)
        self.assertFalse(tbed.invalid, tbed.invalid_reason)
        self.assertEqual(tbed.phase, 1)
        self.assertTrue(tbed.transcriptomic)
        self.assertTrue(tbed.has_stop_codon)

    def test_diex_neg_bed_with_phase_one(self):

        string = "1\t10\t300\tID=test;phase=1;coding=True\t0\t-\t70\t300\t0\t2\t90,100\t0,190"

        seq = "A" + "CGG" * 42 + "TAA" + "A" * 60

        bed = BED12(string, transcriptomic=False)
        self.assertFalse(bed.header)
        self.assertEqual(bed.thick_start, 71)
        self.assertFalse(bed.invalid)
        self.assertEqual(bed.name, "test")
        self.assertEqual(bed.phase, 1)

        tbed = bed.to_transcriptomic(sequence=seq)

        self.assertFalse(tbed.has_start_codon)
        self.assertEqual(tbed.thick_start, 1)
        self.assertEqual(tbed.thick_end, 130)
        self.assertEqual(tbed.phase, 1)
        self.assertFalse(tbed.invalid, tbed.invalid_reason)
        self.assertTrue(tbed.transcriptomic)
        self.assertTrue(tbed.has_stop_codon)

    def test_diex_pos_bed_with_phase_two(self):

        string = "1\t9\t111\tID=test;phase=2;coding=True\t0\t+\t9\t101\t0\t1\t102\t0"

        seq = "GA" + "CGG" * 29 + "TAA" + "A" * 10

        bed = BED12(string, transcriptomic=False)
        self.assertFalse(bed.header)
        self.assertEqual(bed.thick_start, 10)
        self.assertFalse(bed.invalid)
        self.assertEqual(bed.name, "test")
        self.assertEqual(bed.phase, 2)

        tbed = bed.to_transcriptomic(sequence=seq)

        self.assertFalse(tbed.has_start_codon)
        self.assertEqual(tbed.thick_start, 1)
        self.assertEqual(tbed.thick_end, 92)
        self.assertFalse(tbed.invalid, tbed.invalid_reason)
        self.assertEqual(tbed.phase, 2)
        self.assertTrue(tbed.transcriptomic)
        self.assertTrue(tbed.has_stop_codon)

    def test_diex_neg_bed_with_phase_two(self):

        string = "1\t10\t301\tID=test;phase=2;coding=True\t0\t-\t70\t301\t0\t2\t90,101\t0,190"

        seq = "GA" + "CGG" * 42 + "TAA" + "A" * 60

        bed = BED12(string, transcriptomic=False)
        self.assertFalse(bed.header)
        self.assertEqual(bed.thick_start, 71)
        self.assertFalse(bed.invalid)
        self.assertEqual(bed.name, "test")
        self.assertEqual(bed.phase, 2)

        tbed = bed.to_transcriptomic(sequence=seq)

        self.assertFalse(tbed.has_start_codon)
        self.assertEqual(tbed.thick_start, 1)
        self.assertEqual(tbed.thick_end, 131)
        self.assertEqual(tbed.phase, 2)
        self.assertFalse(tbed.invalid, tbed.invalid_reason)
        self.assertTrue(tbed.transcriptomic)
        self.assertTrue(tbed.has_stop_codon)

    def test_tran_to_bed12_neg(self):

        for end, phase in [(299, 0), (300, 1), (301, 2)]:
            with self.subTest():
                t = Transcript()
                t.chrom = "1"
                t.start, t.end, t.strand, t.id = 11, end, "-", "test"

                t.add_exons([(11, 100), (201, end)])
                t.add_exons([(71, 100), (201, end)], features="CDS", phases=[0, phase])
                t.finalize()
                r = t.as_bed12()
                self.assertEqual(r.name, "ID={};coding={};phase={}".format(t.id, True, phase))
                self.assertEqual(r.phase, phase, (end, phase))
                self.assertEqual(r.thick_end, end)
                self.assertFalse(r.invalid)

    def test_touching_exons(self):

        bed12line = ["chr1", 172601, 175626, "ID=foo.1", 100, "-", 172601, 175626, "0,0,0", 3, "199,1281,861,",
                     "0,199,2164"]
        bed = BED12(bed12line, transcriptomic=False)
        self.assertFalse(bed.invalid, bed.invalid_reason)
        t = Transcript(bed)
        t.finalize()
        self.assertEqual(t.exon_num, 2)  # The touching exons should have been merged
        self.assertEqual(t.exons, [(172602, 174081), (174766, 175626)], t.exons)


if __name__ == "__main__":
    unittest.main()