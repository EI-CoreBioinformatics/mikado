import operator
import os
import pickle
import re
import unittest
import random
import pandas as pd
from pytest import mark

from Mikado.exceptions import InvalidTranscript, InvalidCDS
from Mikado._transcripts.transcript_methods.finalizing import _check_completeness, _check_internal_orf, \
    _check_phase_correctness, _calculate_phases, _calculate_introns, _basic_final_checks, _verify_boundaries, \
    _check_cdna_vs_utr, _fix_stop_codon
from Mikado.configuration import MikadoConfiguration
from sqlalchemy.engine import reflection
import itertools
from Mikado.utilities import Interval, IntervalTree
from ..configuration.configurator import load_and_validate_config
from ..loci import Transcript
from ..parsers.bed12 import BED12
from ..parsers.GTF import GtfLine
from ..parsers.GFF import GffLine
from ..transcripts.transcript import Metric
from Mikado.transcripts.transcript_methods import retrieval
from ..utilities.log_utils import create_default_logger


class TestMetricClass(unittest.TestCase):

    def test_wrong_set(self):
        class o(object):
            def __init__(self):
                self.__test = None

            def getter(self):
                return self.__test

            def setter(self, value):
                if not isinstance(value, int):
                    raise ValueError
                self.__test = value

            def deleter(self):
                self.__test = None

            test = Metric(fget=getter, fset=setter, fdel=deleter)

        p = o()

        with self.assertRaises(ValueError):
            p.test = '20'

        self.assertIs(p.test, None)
        p.test = 50
        self.assertEqual(p.test, 50)
        del p.test
        self.assertIs(p.test, None)

    def test_wrong_raw(self):

        with self.assertRaises(ValueError):
            class o:
                def __init__(self):
                    self.__test = 0

                @Metric
                def test(self):
                    return self.__test

                test.usable_raw = "foo"

        class o2:
            def __init__(self):
                self.__test = 0

            @Metric
            def test(self):
                return self.__test

            test.usable_raw = True

        p = o2()
        self.assertTrue(hasattr(p, "test"))
        self.assertTrue(getattr(o2, "test").usable_raw)

    def test_wrong_rtype(self):
        valid = ["int", "float", "str", "bytes",
                     "complex", "dict", "bytearray", "set", "frozenset",
                     "list", "bool", "object", "slice", "tuple", None]

        for typ in valid:
            for setter in (str, bytes):
                with self.subTest(typ=typ, setter=setter):
                    if setter == bytes and typ is not None:
                        typ = typ.encode()

                    class o:
                        def __init__(self):
                            self.__test = 0

                        @Metric
                        def test(self):
                            return self.__test
                        test.rtype = typ

        for invalid in [0, "foo", dict(), unittest.TestCase, 0j]:
            with self.subTest(typ=invalid):
                with self.assertRaises(ValueError):
                    class o:
                        def __init__(self):
                            self.__test = 0

                        @Metric
                        def test(self):
                            return self.__test

                        test.rtype = invalid

    def test_external_score(self):
        conf = MikadoConfiguration()
        conf.prepare.files.source_score = {"at": 5, "tr": -1, "pb": 1, "st": 0}
        t = Transcript(configuration=conf)
        for source in conf.prepare.files.source_score.keys():
            t.original_source = source
            self.assertEqual(t.source_score, conf.prepare.files.source_score[source])

        t.original_source = "foo"
        self.assertEqual(t.source_score, 0)

    def test_pickling_unpickling_metrics(self):
        bed_line = "Chr5\t26585506\t26586850\tID=c58_g1_i2.mrna1.35;coding=False\t99.0\t+\t26585506\t26585507\t0\t5\t383,121,78,105,213\t0,475,710,913,1131"
        conf = MikadoConfiguration()
        conf.prepare.files.source_score = {"at": 5, "tr": -1, "pb": 1, "st": 0}
        t = Transcript(bed_line, source="tr", configuration=conf)
        t.finalize()
        for metrics in t.get_available_metrics():
            try:
                rtype = operator.attrgetter("{metric}.rtype")(Transcript)
            except AttributeError:
                continue
            if rtype == "float":
                value = random.random()
            elif rtype == "bool":
                value = random.choice([True, False])
            elif rtype == "int":
                value = random.randint(0, 1000)
            try:
                setattr(t, metrics, value)
            except AttributeError:
                continue
        u = pickle.loads(pickle.dumps(t))
        for metric in t.get_available_metrics():
            original, new = getattr(t, metric), getattr(u, metric)
            self.assertEqual(original, new, (metric, original, new))

    def test_undefined_strand_multi(self):
        bed_line = "Chr5\t26585506\t26586850\tID=c58_g1_i2.mrna1.35;coding=False\t99.0\t.\t26585506\t26585507\t0\t5\t383,121,78,105,213\t0,475,710,913,1131"
        for tentative in [False, True]:
            with self.subTest(tentative=tentative):
                if tentative is False:
                    with self.assertRaises(InvalidTranscript):
                        t = Transcript(bed_line, source="tr", accept_undefined_multi=tentative)
                else:
                    t = Transcript(bed_line, source="tr", accept_undefined_multi=tentative)
                    t.finalize()
                    self.assertEqual(t.strand, None)

    def test_cds_bridging_two_exons(self):
        t = Transcript()
        t.chrom, t.start, t.end, t.strand, t.id = "Chr1", 1001, 1500, "+", "foo"
        t.exons = [(1001, 1200), (1301, 1500)]
        t.combined_cds = [(1101, 1400)]
        t.logger = create_default_logger("test_cds_bridging_two_exons", level="DEBUG")
        with self.assertLogs("test_cds_bridging_two_exons") as cmo:
            t.finalize()
        self.assertFalse(t.is_coding, cmo.output)

    def test_overlapping_exons(self):
        t = Transcript()
        t.chrom, t.start, t.end, t.strand, t.id = "Chr1", 1001, 1500, "+", "foo"
        t.exons = [(1001, 1400), (1301, 1500)]
        t.logger = create_default_logger("test_overlapping_exons", level="DEBUG")
        with self.assertRaises(InvalidTranscript):
            t.finalize()

    def test_category(self):

        correct = ["External", None, "Internal", "SuperCali", "My_test"]
        for corr, setter in itertools.product(correct, (str, bytes)):
            with self.subTest(corr=corr, setter=setter):
                if corr is not None and setter == bytes:
                    corr = corr.encode()

                class o:
                    def __init__(self):
                        self.__test = 0

                    @Metric
                    def test(self):
                        return self.__test

                    test.category = corr

        incorr = [{}, dict, 0, 50j, 0.3, []]
        for inc in incorr:
            with self.subTest(inc=inc):
                with self.assertRaises(TypeError):
                    class o:
                        def __init__(self):
                            self.__test = 0

                        @Metric
                        def test(self):
                            return self.__test
                        test.category = inc


class TestTranscriptInit(unittest.TestCase):

    def test_wrong_data(self):

        incorrect = [0, dict(), 0j, [], tuple([]), set(), ""]

        for inc in incorrect:
            with self.subTest():
                with self.assertRaises(TypeError):
                    t = Transcript(inc)

    def test_initialize_bed12_string(self):

        t = "Chr1\t100\t1000\tID=t1;coding=False\t0\t+\t100\t1000\t0\t1\t900\t0"
        self.assertFalse(BED12(t).header)
        tr = Transcript(t)
        self.assertEqual(tr.start, BED12(t).start)
        self.assertEqual(tr.chrom, BED12(t).chrom)
        tr = Transcript(t.encode())
        self.assertEqual(tr.start, BED12(t).start)
        self.assertEqual(tr.chrom, BED12(t).chrom)

    def test_equality(self):
        t = "Chr1\t100\t1000\tID=t1;coding=False\t0\t+\t100\t1000\t0\t1\t900\t0"
        self.assertFalse(BED12(t).header)
        tr = Transcript(t)
        self.assertNotEqual(tr, "")
        self.assertNotEqual(tr, t)
        self.assertNotEqual(tr, BED12(t))
        self.assertEqual(tr, Transcript(t))

    def test_comparison(self):
        t = "Chr1\t100\t1000\tID=t1;coding=False\t0\t+\t100\t1000\t0\t1\t900\t0"
        tr = Transcript(t)
        with self.assertRaises(TypeError):
            _ = (tr < "")

        t2 = "Chr2\t100\t1000\tID=t1;coding=False\t0\t+\t100\t1000\t0\t1\t900\t0"
        self.assertLess(tr, Transcript(t2))
        t2 = "Chr1\t101\t1000\tID=t1;coding=False\t0\t+\t101\t1000\t0\t1\t899\t0"
        self.assertLess(tr, Transcript(t2))
        t2 = "Chr1\t100\t1001\tID=t1;coding=False\t0\t+\t100\t1001\t0\t1\t901\t0"
        self.assertLess(tr, Transcript(t2))

        # Now when they are coding.
        # To prevent the whole sublocus situation to crumble, we have to avoid considering
        # the coding part ..
        t2 = "Chr1\t100\t1000\tID=t1;coding=True\t0\t+\t600\t900\t0\t1\t900\t0"
        self.assertFalse(tr < Transcript(t2))

        t = "Chr1\t100\t1000\tID=t1;coding=True\t0\t+\t200\t5000\t0\t1\t900\t0"
        tr = Transcript(t)
        self.assertFalse(tr < Transcript(t2))

    def test_get_orfs(self):
        t = "Chr1\t100\t1000\tID=t1;coding=False\t0\t+\t100\t1000\t0\t1\t900\t0"
        tr = Transcript(t)
        orfs = list(tr.get_internal_orf_beds())
        self.assertEqual(len(orfs), 1)
        orf = orfs.pop()
        self.assertIsInstance(orf, BED12)
        self.assertFalse(orf.coding)
        self.assertEqual(tr._internal_orfs_transcripts, [])


class WrongLoadedOrf(unittest.TestCase):

    def setUp(self):

        self.tr = Transcript()
        self.tr.start, self.tr.end, self.tr.chrom, self.tr.strand = (101, 1000, "Chr1", "+")
        self.tr.id = "test1"
        self.tr.add_exons([(101, 400), (701, 1000)])
        self.tr.finalize()

    def test_load_invalid_length(self):

        b_invalid = BED12(transcriptomic=True)
        b_invalid.chrom = self.tr.id
        self.assertTrue(b_invalid.transcriptomic)
        # b_invalid.name = self.tr.id
        b_invalid.start = 0
        b_invalid.strand = "+"
        b_invalid.end = self.tr.cdna_length + 10
        b_invalid.thick_start = 101
        b_invalid.thick_end = 190
        self.assertEqual(b_invalid.chrom,
                         b_invalid.id,
                         b_invalid.id)

        with self.assertLogs("null", "WARNING") as cm:
            retrieval.load_orfs(self.tr, [b_invalid])

        found_message = False
        for _ in cm.output:
            if "Wrong ORF for {}:".format(self.tr.id) in _:
                found_message = True
                break

        self.assertTrue(found_message, cm.output)

    def test_load_invalid_multiple(self):

        b_valid = BED12(transcriptomic=True)
        b_valid.chrom = self.tr.id
        b_valid.name = "valid"
        b_valid.start, b_valid.end, b_valid.strand = 0, self.tr.cdna_length - 1, "+"
        b_valid.thick_start, b_valid.thick_end = 101, 190

        b_invalid = b_valid.copy()
        b_invalid.name = "invalid"
        b_invalid.thick_start = 1
        b_invalid.thick_end = 89
        b_invalid.phase = 0

        self.assertTrue(b_invalid.invalid)
        self.assertFalse(b_valid.invalid, b_valid.invalid_reason)

        with self.assertLogs("null", "DEBUG") as _:
            retrieval.load_orfs(self.tr, [b_valid, b_invalid])

        self.assertEqual(self.tr.number_internal_orfs, 1)

    def test_filter_non_transcriptomic(self):

        b_valid = BED12(transcriptomic=True)
        b_valid.chrom = self.tr.id
        b_valid.name = "valid"
        b_valid.start, b_valid.end, b_valid.strand = 0, self.tr.cdna_length - 1, "+"
        b_valid.thick_start, b_valid.thick_end = 101, 190

        b_invalid = b_valid.copy()
        b_invalid.name = "non-transcriptomic"
        b_invalid.transcriptomic = False

        retained = retrieval.find_overlapping_cds(self.tr, [b_invalid, b_valid])
        self.assertEqual(retained, [b_valid])


class TestAsBed12(unittest.TestCase):

    def test_casePositive(self):

        tr = Transcript()
        tr.chrom, tr.start, tr.end, tr.strand = "Chr1", 101, 3000, "+"
        tr.id = "test1"
        tr.add_exons([(101, 300),
                      (401, 600),
                      (801, 1200),
                      (2501, 3000)
                      ])

        tr.add_exons([(421, 600),  # 180
                      (801, 1200),  # 400
                      (2501, 2700)  # 200  = 780 % 3 == 0
                      ], features="CDS")
        with self.assertLogs("null", "DEBUG") as _:
            tr.finalize()
        self.assertTrue(tr.is_coding)

        b12 = tr.as_bed12()
        self.assertEqual(b12.thick_start, tr.combined_cds_start)
        self.assertEqual(b12.thick_end, tr.combined_cds_end)
        self.assertEqual(len(b12.block_sizes), tr.exon_num)
        self.assertTrue((b12.block_sizes == [200, 200, 400, 500]).all(),
                         b12.block_sizes)
        self.assertEqual(b12.strand, "+")
        self.assertTrue((b12.block_starts == [0, 300, 700, 2400]).all(),
                         b12.block_starts)
        self.assertEqual(str(b12),
                         "\t".join([str(_) for _ in
                                    ["Chr1", 100, 3000,
                                     "ID={};coding=True;phase=0".format(tr.id),
                                     0, tr.strand,
                                     b12.thick_start - 1, b12.thick_end,
                                     0, 4,
                                     ",".join([str(__) for __ in [200, 200, 400, 500]]),
                                     ",".join([str(___) for ___ in [0, 300, 700, 2400]])]]
                                   ))

    def test_caseNegative(self):
        tr = Transcript()
        tr.chrom, tr.start, tr.end, tr.strand = "Chr1", 101, 3000, "-"
        tr.id = "test1"
        tr.add_exons([(101, 300),
                      (401, 600),
                      (801, 1200),
                      (2501, 3000)
                      ])

        tr.add_exons([(421, 600),  # 180
                      (801, 1200),  # 400
                      (2501, 2700)  # 200  = 780 % 3 == 0
                      ], features="CDS")
        with self.assertLogs("null", "DEBUG") as _:
            tr.finalize()
        self.assertTrue(tr.is_coding)

        b12 = tr.as_bed12()
        self.assertEqual(b12.thick_start, tr.combined_cds_end)
        self.assertEqual(b12.thick_end, tr.combined_cds_start)
        self.assertEqual(len(b12.block_sizes), tr.exon_num)
        self.assertTrue((b12.block_sizes == [200, 200, 400, 500]).all(), b12.block_sizes)
        self.assertEqual(b12.strand, "-")
        self.assertTrue((b12.block_starts == [0, 300, 700, 2400]).all(), b12.block_starts)

        self.assertEqual(tr.format("bed12"), str(b12))
        self.assertEqual(str(b12),
                         "\t".join([str(_) for _ in
                                    ["Chr1", 100, 3000,
                                     "ID={};coding=True;phase=0".format(tr.id),
                                     0, tr.strand,
                                     b12.thick_start - 1, b12.thick_end,
                                     0, 4,
                                     ",".join([str(__) for __ in [200, 200, 400, 500]]),
                                     ",".join([str(___) for ___ in [0, 300, 700, 2400]])]]
                                   ))


class TestPrintTranscriptomic(unittest.TestCase):
    
    def setUp(self):
        transcript_lines = """
Chr5\t26581217\t26581531\tID=st_Stringtie_STAR.21709.1;coding=False\t0\t+\t26581217\t26581218\t0\t1\t314\t0
Chr5\t26584773\t26587782\tID=at_AT5G66610.2;coding=True;phase=0\t0\t+\t26585222\t26587755\t0\t10\t106,54,545,121,78,105,213,63,119,496\t0,446,571,1208,1443,1646,1864,2160,2310,2513
Chr5\t26574999\t26578012\tID=at_AT5G66600.2;coding=True;phase=0\t0\t-\t26575104\t26577954\t0\t9\t411,126,87,60,100,809,126,72,157\t0,495,711,885,1035,1261,2163,2378,2856
Chr5\t26574999\t26578625\tID=at_AT5G66600.3;coding=True;phase=0\t0\t-\t26575104\t26578315\t0\t11\t411,126,87,60,100,809,126,72,82,188,107\t0,495,711,885,1035,1261,2163,2378,2856,3239,3519
        """

        self.transcripts = dict() 
        for line in transcript_lines.split("\n"):
            if line.strip() == "":
                continue
            transcript = Transcript(BED12(line))
            transcript.finalize()
            transcript.parent = transcript.id + "_gene"
            transcript.name = transcript.id
            self.transcripts[transcript.id] = transcript

        self.assertEqual(len(self.transcripts), 4)

    def test_bed12(self):

        results = {
            "st_Stringtie_STAR.21709.1": [
                "st_Stringtie_STAR.21709.1\t0\t314\tID=st_Stringtie_STAR.21709.1;coding=False\t0\t+\t0\t1\t0\t1\t314\t0",
                "st_Stringtie_STAR.21709.1\t0\t314\tID=st_Stringtie_STAR.21709.1;coding=False\t0\t+\t0\t1\t0\t1\t314\t0"
            ],
            "at_AT5G66610.2": [
                "at_AT5G66610.2\t0\t1900\tID=at_AT5G66610.2;coding=True;phase=0\t0\t+\t109\t1873\t0\t10\t106,54,545,121,78,105,213,63,119,496\t0,106,160,705,826,904,1009,1222,1285,1404",
                "at_AT5G66610.2\t0\t1900\tID=at_AT5G66610.2;coding=False\t0\t+\t0\t1\t0\t10\t106,54,545,121,78,105,213,63,119,496\t0,106,160,705,826,904,1009,1222,1285,1404",
            ],
            "at_AT5G66600.2": [
                "at_AT5G66600.2\t0\t1948\tID=at_AT5G66600.2;coding=True;phase=0\t0\t+\t58\t1843\t0\t9\t157,72,126,809,100,60,87,126,411\t0,157,229,355,1164,1264,1324,1411,1537",
                "at_AT5G66600.2\t0\t1948\tID=at_AT5G66600.2;coding=False\t0\t+\t0\t1\t0\t9\t157,72,126,809,100,60,87,126,411\t0,157,229,355,1164,1264,1324,1411,1537"
            ],
            "at_AT5G66600.3": [
                "at_AT5G66600.3\t0\t2168\tID=at_AT5G66600.3;coding=True;phase=0\t0\t+\t218\t2063\t0\t11\t107,188,82,72,126,809,100,60,87,126,411\t0,107,295,377,449,575,1384,1484,1544,1631,1757",
                "at_AT5G66600.3\t0\t2168\tID=at_AT5G66600.3;coding=False\t0\t+\t0\t1\t0\t11\t107,188,82,72,126,809,100,60,87,126,411\t0,107,295,377,449,575,1384,1484,1544,1631,1757"
            ]
        }

        for tid, transcript in self.transcripts.items():
            with self.subTest(tid=tid):
                with_cds = transcript.format(format_name="bed12", with_cds=True, transcriptomic=True)
                without_cds = transcript.format(format_name="bed12", with_cds=False, transcriptomic=True)
                self.assertEqual(with_cds, results[tid][0])
                self.assertEqual(without_cds, results[tid][1])

    def test_gff3(self):

        results = {
            "st_Stringtie_STAR.21709.1": [
                """st_Stringtie_STAR.21709.1\tbed12\ttranscript\t1\t314\t0.0\t+\t.\tID=st_Stringtie_STAR.21709.1;Parent=st_Stringtie_STAR.21709.1_gene;Name=st_Stringtie_STAR.21709.1
st_Stringtie_STAR.21709.1\tbed12\texon\t1\t314\t.\t+\t.\tID=st_Stringtie_STAR.21709.1.exon1;Parent=st_Stringtie_STAR.21709.1""",
                """st_Stringtie_STAR.21709.1\tbed12\ttranscript\t1\t314\t0.0\t+\t.\tID=st_Stringtie_STAR.21709.1;Parent=st_Stringtie_STAR.21709.1_gene;Name=st_Stringtie_STAR.21709.1
st_Stringtie_STAR.21709.1\tbed12\texon\t1\t314\t.\t+\t.\tID=st_Stringtie_STAR.21709.1.exon1;Parent=st_Stringtie_STAR.21709.1"""
            ],
            "at_AT5G66610.2":
                [
                    """at_AT5G66610.2\tbed12\tmRNA\t1\t1900\t0.0\t+\t.\tID=at_AT5G66610.2;Parent=at_AT5G66610.2_gene;Name=at_AT5G66610.2
at_AT5G66610.2\tbed12\texon\t1\t1900\t.\t+\t.\tID=at_AT5G66610.2.exon1;Parent=at_AT5G66610.2
at_AT5G66610.2\tbed12\tfive_prime_UTR\t1\t109\t.\t+\t.\tID=at_AT5G66610.2.five_prime_UTR1;Parent=at_AT5G66610.2
at_AT5G66610.2\tbed12\tCDS\t110\t1873\t.\t+\t0\tID=at_AT5G66610.2.CDS1;Parent=at_AT5G66610.2
at_AT5G66610.2\tbed12\tthree_prime_UTR\t1874\t1900\t.\t+\t.\tID=at_AT5G66610.2.three_prime_UTR1;Parent=at_AT5G66610.2""",
                    """at_AT5G66610.2\tbed12\ttranscript\t1\t1900\t0.0\t+\t.\tID=at_AT5G66610.2;Parent=at_AT5G66610.2_gene;Name=at_AT5G66610.2
at_AT5G66610.2\tbed12\texon\t1\t1900\t.\t+\t.\tID=at_AT5G66610.2.exon1;Parent=at_AT5G66610.2"""
                ],
            "at_AT5G66600.2": ["""at_AT5G66600.2\tbed12\tmRNA\t1\t1948\t0.0\t+\t.\tID=at_AT5G66600.2;Parent=at_AT5G66600.2_gene;Name=at_AT5G66600.2
at_AT5G66600.2\tbed12\texon\t1\t1948\t.\t+\t.\tID=at_AT5G66600.2.exon1;Parent=at_AT5G66600.2
at_AT5G66600.2\tbed12\tfive_prime_UTR\t1\t58\t.\t+\t.\tID=at_AT5G66600.2.five_prime_UTR1;Parent=at_AT5G66600.2
at_AT5G66600.2\tbed12\tCDS\t59\t1843\t.\t+\t0\tID=at_AT5G66600.2.CDS1;Parent=at_AT5G66600.2
at_AT5G66600.2\tbed12\tthree_prime_UTR\t1844\t1948\t.\t+\t.\tID=at_AT5G66600.2.three_prime_UTR1;Parent=at_AT5G66600.2""",
                 """at_AT5G66600.2\tbed12\ttranscript\t1\t1948\t0.0\t+\t.\tID=at_AT5G66600.2;Parent=at_AT5G66600.2_gene;Name=at_AT5G66600.2
at_AT5G66600.2\tbed12\texon\t1\t1948\t.\t+\t.\tID=at_AT5G66600.2.exon1;Parent=at_AT5G66600.2"""],
            "at_AT5G66600.3": ["""at_AT5G66600.3\tbed12\tmRNA\t1\t2168\t0.0\t+\t.\tID=at_AT5G66600.3;Parent=at_AT5G66600.3_gene;Name=at_AT5G66600.3
at_AT5G66600.3\tbed12\texon\t1\t2168\t.\t+\t.\tID=at_AT5G66600.3.exon1;Parent=at_AT5G66600.3
at_AT5G66600.3\tbed12\tfive_prime_UTR\t1\t218\t.\t+\t.\tID=at_AT5G66600.3.five_prime_UTR1;Parent=at_AT5G66600.3
at_AT5G66600.3\tbed12\tCDS\t219\t2063\t.\t+\t0\tID=at_AT5G66600.3.CDS1;Parent=at_AT5G66600.3
at_AT5G66600.3\tbed12\tthree_prime_UTR\t2064\t2168\t.\t+\t.\tID=at_AT5G66600.3.three_prime_UTR1;Parent=at_AT5G66600.3""",
                 """at_AT5G66600.3\tbed12\ttranscript\t1\t2168\t0.0\t+\t.\tID=at_AT5G66600.3;Parent=at_AT5G66600.3_gene;Name=at_AT5G66600.3
at_AT5G66600.3\tbed12\texon\t1\t2168\t.\t+\t.\tID=at_AT5G66600.3.exon1;Parent=at_AT5G66600.3"""]
        }

        self.maxDiff = None
        for tid, transcript in self.transcripts.items():
            with self.subTest(tid=tid):
                with_cds = transcript.format(format_name="gff3", with_cds=True, transcriptomic=True).strip()
                without_cds = transcript.format(format_name="gff3", with_cds=False, transcriptomic=True).strip()
                self.assertEqual(with_cds, results[tid][0])
                self.assertEqual(without_cds, results[tid][1])

#     def test_gtf(self):
#         output = """
# Chr5    bed12   transcript      26581218        26581531        0.0     +       .       gene_id "gene_1"; transcript_id "st_Stringtie_STAR.21709.1"; Name "st_Stringtie_STAR.21709.1";
# Chr5    bed12   exon    26581218        26581531        .       +       .       gene_id "gene_1"; transcript_id "st_Stringtie_STAR.21709.1";
# at_AT5G66610.2  bed12   mRNA    1       1900    0.0     +       .       gene_id "gene_2"; transcript_id "at_AT5G66610.2"; Name "at_AT5G66610.2";
# at_AT5G66610.2  bed12   exon    1       1900    0.0     +       .       gene_id "gene_2"; transcript_id "at_AT5G66610.2"; Name "at_AT5G66610.2";
# Chr5    bed12   5UTR    1       109     .       +       .       gene_id "gene_2"; transcript_id "at_AT5G66610.2";
# Chr5    bed12   CDS     110     1873    .       +       0       gene_id "gene_2"; transcript_id "at_AT5G66610.2";
# Chr5    bed12   3UTR    1874    1900    .       +       .       gene_id "gene_2"; transcript_id "at_AT5G66610.2";
# at_AT5G66600.2  bed12   mRNA    1       1948    0.0     +       .       gene_id "gene_3"; transcript_id "at_AT5G66600.2"; Name "at_AT5G66600.2";
# at_AT5G66600.2  bed12   exon    1       1948    0.0     +       .       gene_id "gene_3"; transcript_id "at_AT5G66600.2"; Name "at_AT5G66600.2";
# Chr5    bed12   3UTR    1       58      .       -       .       gene_id "gene_3"; transcript_id "at_AT5G66600.2";
# Chr5    bed12   CDS     59      1843    .       -       0       gene_id "gene_3"; transcript_id "at_AT5G66600.2";
# Chr5    bed12   5UTR    1844    1948    .       -       .       gene_id "gene_3"; transcript_id "at_AT5G66600.2";
# at_AT5G66600.3  bed12   mRNA    1       2168    0.0     +       .       gene_id "gene_4"; transcript_id "at_AT5G66600.3"; Name "at_AT5G66600.3";
# at_AT5G66600.3  bed12   exon    1       2168    0.0     +       .       gene_id "gene_4"; transcript_id "at_AT5G66600.3"; Name "at_AT5G66600.3";
# Chr5    bed12   3UTR    1       218     .       -       .       gene_id "gene_4"; transcript_id "at_AT5G66600.3";
# Chr5    bed12   CDS     219     2063    .       -       0       gene_id "gene_4"; transcript_id "at_AT5G66600.3";
# Chr5    bed12   5UTR    2064    2168    .       -       .       gene_id "gene_4"; transcript_id "at_AT5G66600.3";
#     """


class TestPrintIntrons(unittest.TestCase):

    def test_not_coding(self):

        tr = Transcript()
        tr.chrom = "Chr1"
        tr.start = 101
        tr.end = 2000
        tr.strand = "+"
        tr.add_exons([(101, 300),
                      (1701, 2000)])
        tr.id = "test1"
        tr.parent = "gene1"
        tr.finalize()

        gff = tr.format("gff", with_introns=True)
        self.maxDiff = None
        res = """Chr1\tMikado\ttranscript\t101\t2000\t.\t+\t.\tID=test1;Parent=gene1
Chr1\tMikado\texon\t101\t300\t.\t+\t.\tID=test1.exon1;Parent=test1
Chr1\tMikado\tintron\t301\t1700\t.\t+\t.\tID=test1.intron1;Parent=test1
Chr1\tMikado\texon\t1701\t2000\t.\t+\t.\tID=test1.exon2;Parent=test1"""
        self.assertEqual(gff,
                         res,
                         "++++\n\n"+"\n+++\n".join([gff, res]))

        gtf = tr.format("gtf", with_introns=True)
        res = """Chr1\tMikado\ttranscript\t101\t2000\t.\t+\t.\tgene_id "gene1"; transcript_id "test1";
Chr1\tMikado\texon\t101\t300\t.\t+\t.\tgene_id "gene1"; transcript_id "test1";
Chr1\tMikado\tintron\t301\t1700\t.\t+\t.\tgene_id "gene1"; transcript_id "test1";
Chr1\tMikado\texon\t1701\t2000\t.\t+\t.\tgene_id "gene1"; transcript_id "test1";"""
        self.assertEqual(gtf, res,
                         "++++\n\n" + "\n+++\n".join([gtf, res]))

    def test_non_coding_negative(self):
        tr = Transcript()
        tr.chrom = "Chr1"
        tr.start = 101
        tr.end = 2000
        tr.strand = "-"
        tr.add_exons([(101, 300),
                      (1701, 2000)])
        tr.id = "test1"
        tr.parent = "gene1"
        tr.finalize()

        gff = tr.format("gff", with_introns=True)
        self.maxDiff = None
        res = """Chr1\tMikado\ttranscript\t101\t2000\t.\t-\t.\tID=test1;Parent=gene1
Chr1\tMikado\texon\t101\t300\t.\t-\t.\tID=test1.exon1;Parent=test1
Chr1\tMikado\tintron\t301\t1700\t.\t-\t.\tID=test1.intron1;Parent=test1
Chr1\tMikado\texon\t1701\t2000\t.\t-\t.\tID=test1.exon2;Parent=test1"""
        self.assertEqual(gff,
                         res,
                         "++++\n\n"+"\n+++\n".join([gff, res]))

        gtf = tr.format("gtf", with_introns=True)
        res = """Chr1\tMikado\ttranscript\t101\t2000\t.\t-\t.\tgene_id "gene1"; transcript_id "test1";
Chr1\tMikado\texon\t101\t300\t.\t-\t.\tgene_id "gene1"; transcript_id "test1";
Chr1\tMikado\tintron\t301\t1700\t.\t-\t.\tgene_id "gene1"; transcript_id "test1";
Chr1\tMikado\texon\t1701\t2000\t.\t-\t.\tgene_id "gene1"; transcript_id "test1";"""
        self.assertEqual(gtf, res,
                         "++++\n\n" + "\n+++\n".join([gtf, res]))

    def test_coding_positive(self):
        tr = Transcript()
        tr.chrom = "Chr1"
        tr.start = 101
        tr.end = 2000
        tr.strand = "+"
        tr.add_exons([(101, 300),
                      (1701, 2000)])
        tr.add_exons([(101, 300),
                      (1701, 2000)], features="CDS")
        tr.id = "test1"
        tr.parent = "gene1"

        gff = tr.format("gff", with_introns=True)
        self.maxDiff = None
        res = """Chr1\tMikado\tmRNA\t101\t2000\t.\t+\t.\tID=test1;Parent=gene1;Name=test1
Chr1\tMikado\tCDS\t101\t300\t.\t+\t0\tID=test1.CDS1;Parent=test1
Chr1\tMikado\texon\t101\t300\t.\t+\t.\tID=test1.exon1;Parent=test1
Chr1\tMikado\tintron\t301\t1700\t.\t+\t.\tID=test1.intron1;Parent=test1
Chr1\tMikado\tCDS\t1701\t2000\t.\t+\t1\tID=test1.CDS2;Parent=test1
Chr1\tMikado\texon\t1701\t2000\t.\t+\t.\tID=test1.exon2;Parent=test1"""
        self.assertEqual(gff,
                         res,
                         "++++\n\n" + "\n+++\n".join([gff, res]))

        gtf = tr.format("gtf", with_introns=True)
        res = """Chr1\tMikado\tmRNA\t101\t2000\t.\t+\t.\tgene_id "gene1"; transcript_id "test1"; Name "test1";
Chr1\tMikado\tCDS\t101\t300\t.\t+\t0\tgene_id "gene1"; transcript_id "test1";
Chr1\tMikado\texon\t101\t300\t.\t+\t.\tgene_id "gene1"; transcript_id "test1";
Chr1\tMikado\tintron\t301\t1700\t.\t+\t.\tgene_id "gene1"; transcript_id "test1";
Chr1\tMikado\tCDS\t1701\t2000\t.\t+\t1\tgene_id "gene1"; transcript_id "test1";
Chr1\tMikado\texon\t1701\t2000\t.\t+\t.\tgene_id "gene1"; transcript_id "test1";"""
        self.assertEqual(gtf, res,
                         "++++\n\n" + "\n+++\n".join([gtf, res]))
        
    def test_coding_negative(self):
        tr = Transcript()
        tr.chrom = "Chr1"
        tr.start = 101
        tr.end = 2000
        tr.strand = "-"
        tr.add_exons([(101, 300),
                      (1701, 2000)])
        tr.add_exons([(101, 300),
                      (1701, 2000)], features="CDS")
        tr.id = "test1"
        tr.parent = "gene1"

        # Phase 0, 0 because the first CDS exon is 300bp
        gff = tr.format("gff", with_introns=True)
        self.maxDiff = None
        res = """Chr1\tMikado\tmRNA\t101\t2000\t.\t-\t.\tID=test1;Parent=gene1;Name=test1
Chr1\tMikado\tCDS\t101\t300\t.\t-\t0\tID=test1.CDS1;Parent=test1
Chr1\tMikado\texon\t101\t300\t.\t-\t.\tID=test1.exon1;Parent=test1
Chr1\tMikado\tintron\t301\t1700\t.\t-\t.\tID=test1.intron1;Parent=test1
Chr1\tMikado\tCDS\t1701\t2000\t.\t-\t0\tID=test1.CDS2;Parent=test1
Chr1\tMikado\texon\t1701\t2000\t.\t-\t.\tID=test1.exon2;Parent=test1"""
        self.assertEqual(gff,
                         res,
                         "++++\n\n"+"\n+++\n".join([gff, res,
                                                    ",\t".join([str(_) for _ in tr.internal_orfs])
                                                    ]
                                                   )
                         )

        gtf = tr.format("gtf", with_introns=True)
        res = """Chr1\tMikado\tmRNA\t101\t2000\t.\t-\t.\tgene_id "gene1"; transcript_id "test1"; Name "test1";
Chr1\tMikado\tCDS\t101\t300\t.\t-\t0\tgene_id "gene1"; transcript_id "test1";
Chr1\tMikado\texon\t101\t300\t.\t-\t.\tgene_id "gene1"; transcript_id "test1";
Chr1\tMikado\tintron\t301\t1700\t.\t-\t.\tgene_id "gene1"; transcript_id "test1";
Chr1\tMikado\tCDS\t1701\t2000\t.\t-\t0\tgene_id "gene1"; transcript_id "test1";
Chr1\tMikado\texon\t1701\t2000\t.\t-\t.\tgene_id "gene1"; transcript_id "test1";"""
        self.assertEqual(gtf, res,
                         "++++\n\n" + "\n+++\n".join([gtf, res]))

    def test_coding_negative_2(self):
        tr = Transcript()
        tr.chrom = "Chr1"
        tr.start = 101
        tr.end = 2000
        tr.strand = "-"
        tr.add_exons([(102, 300),
                      (1701, 1999)])
        tr.add_exons([(102, 300),
                      (1701, 1999)], features="CDS")
        tr.id = "test1"
        tr.parent = "gene1"
        tr.finalize()
        self.assertTrue(tr.is_coding)
        self.maxDiff = 10000

        gtf_res = """Chr1\tMikado\tmRNA\t102\t1999\t.\t-\t.\tgene_id "gene1"; transcript_id "test1"; Name "test1";
Chr1\tMikado\tCDS\t102\t300\t.\t-\t1\tgene_id "gene1"; transcript_id "test1";
Chr1\tMikado\texon\t102\t300\t.\t-\t.\tgene_id "gene1"; transcript_id "test1";
Chr1\tMikado\tCDS\t1701\t1999\t.\t-\t0\tgene_id "gene1"; transcript_id "test1";
Chr1\tMikado\texon\t1701\t1999\t.\t-\t.\tgene_id "gene1"; transcript_id "test1";"""

        gtf = tr.format("gtf")
        self.assertEqual(gtf, gtf_res)

        gff3_res = """Chr1\tMikado\tmRNA\t102\t1999\t.\t-\t.\tID=test1;Parent=gene1;Name=test1
Chr1\tMikado\tCDS\t102\t300\t.\t-\t1\tID=test1.CDS1;Parent=test1
Chr1\tMikado\texon\t102\t300\t.\t-\t.\tID=test1.exon1;Parent=test1
Chr1\tMikado\tCDS\t1701\t1999\t.\t-\t0\tID=test1.CDS2;Parent=test1
Chr1\tMikado\texon\t1701\t1999\t.\t-\t.\tID=test1.exon2;Parent=test1"""

        gff3 = tr.format("gff3")
        self.assertEqual(gff3, gff3_res)

        gff3_cds = [GffLine(_) for _ in gff3.split("\n") if GffLine(_).feature == "CDS"]
        gtf_cds = [GtfLine(_) for _ in gtf.split("\n") if GtfLine(_).feature == "CDS"]
        gff3_cds = dict(((_.start, _.end), _) for _ in gff3_cds)
        gtf_cds = dict(((_.start, _.end), _) for _ in gtf_cds)

        self.assertEqual(gff3_cds.keys(), gtf_cds.keys())
        for key in gff3_cds:
            gtf_8th = str(gtf_cds[key]).split("\t")[7]
            gff_8th = str(gff3_cds[key]).split("\t")[7]
            self.assertEqual(gff_8th, gtf_8th)


class TestRetrieval(unittest.TestCase):
    
    def setUp(self):
        self.tr = Transcript()
        self.tr.chrom = "Chr1"
        self.tr.start = 101
        self.tr.end = 2000
        self.tr.strand = None
        self.tr.add_exons([(101, 2000)])
        self.tr.id = "test1"
        self.tr.parent = "gene1"
        self.tr.finalize()
        conf = load_and_validate_config(os.path.join(
            os.path.dirname(__file__),
            "configuration.yaml"
        ))
        self.assertTrue(conf.pick.chimera_split.blast_check)
        self.assertTrue(conf.pick.chimera_split.execute)
        self.assertEqual(conf.pick.chimera_split.blast_params.leniency, "LENIENT")
        conf.pick.orf_loading.minimal_secondary_orf_length = 50

        self.tr.configuration = conf

    def test_load_pos_and_neg(self):
        
        b1 = BED12(transcriptomic=True)
        b1.chrom = self.tr.id
        b1.start = 0
        b1.end = self.tr.cdna_length - 1
        b1.strand = "+"
        b1.name = "first"
        b1.thick_start = 101
        b1.thick_end = 190
        self.assertFalse(b1.invalid)

        b2 = b1.copy()
        b2.strand = "-"
        b2.thick_start = 1
        b2.thick_end = 87
        b2.name = "second"
        self.assertFalse(b2.invalid)
        with self.assertLogs("null", "DEBUG") as _:
            after_overlap_check = retrieval.find_overlapping_cds(self.tr, [b1, b2])

        self.assertEqual(len(after_overlap_check), 2, self.tr.configuration.pick.orf_loading)
        self.assertEqual(after_overlap_check,
                         [b1, b2],
                         [_.name for _ in after_overlap_check])
        retrieval.load_orfs(self.tr, [b1, b2])
        self.assertEqual(self.tr.number_internal_orfs, 1)
        self.assertEqual(self.tr.combined_cds_start, 201, self.tr.combined_cds_start)
        self.assertEqual(self.tr.combined_cds_length, 90)

    def test_connect(self):

        retrieval._connect_to_db(self.tr)
        reflector = reflection.Inspector.from_engine(self.tr.engine)


class TestMiscellanea(unittest.TestCase):

    def setUp(self):
        gtf = """Chr5	TAIR10	mRNA	5256	5891	.	-	.	ID=AT5G01015.1;Parent=AT5G01015;Name=AT5G01015.1;Index=1
Chr5	TAIR10	five_prime_UTR	5770	5891	.	-	.	Parent=AT5G01015.1
Chr5	TAIR10	CDS	5697	5769	.	-	0	Parent=AT5G01015.1;
Chr5	TAIR10	exon	5697	5891	.	-	.	Parent=AT5G01015.1
Chr5	TAIR10	CDS	5335	5576	.	-	2	Parent=AT5G01015.1;
Chr5	TAIR10	three_prime_UTR	5256	5334	.	-	.	Parent=AT5G01015.1
Chr5	TAIR10	exon	5256	5576	.	-	.	Parent=AT5G01015.1"""

        gtf_lines = [GffLine(_) for _ in gtf.split("\n")]

        self.t1 = Transcript(gtf_lines[0])
        self.t1.add_exons(gtf_lines[1:])
        self.t1.finalize()
        self.assertIs(self.t1.is_coding, True)

    def test_exon_intron_lengths(self):

        t = Transcript()
        t.chrom, t.start, t.end, t.strand, t.id = "Chr1", 101, 1000, "+", "foo.1"
        self.assertEqual(t.max_exon_length, 0)
        self.assertEqual(t.min_exon_length, 0)
        self.assertEqual(t.max_intron_length, 0)
        self.assertEqual(t.min_intron_length, 0)
        t.add_exons([(101, 1000)])
        self.assertEqual(t.max_exon_length, 900)
        self.assertEqual(t.min_exon_length, 900)
        self.assertEqual(t.max_intron_length, 0)
        self.assertEqual(t.min_intron_length, 0)

    def test_suspicious_splicing(self):
        t = Transcript()
        t.chrom, t.start, t.end, t.strand, t.id = "Chr1", 101, 1000, "+", "foo.1"
        t.attributes["canonical_on_reverse_strand"] = "True"
        t.attributes["mixed_splices"] = False
        self.assertTrue(t.suspicious_splicing)

    def test_frames(self):

        self.assertIsInstance(self.t1.frames, dict)
        self.assertEqual(len(self.t1.frames), 3)
        self.assertEqual(set(self.t1.frames.keys()), {0, 1, 2})
        for key in self.t1.frames:
            self.assertEqual(len(self.t1.frames[key]), self.t1.combined_cds_length / 3)

        correct_frames = []
        current = []
        for num in itertools.chain(range(5769, 5697 -1, -1), range(5576, 5335 -1, -1)):
            current.append(num)
            if len(current) == 3:
                current = tuple(current)
                correct_frames.append(current)
                current = []
        current = tuple(current)
        if current:
            correct_frames.append(current)

        self.assertEqual(len(correct_frames), int(self.t1.combined_cds_length / 3), correct_frames)
        self.assertEqual(self.t1.framed_codons, correct_frames)

    def test_wrong_proportions(self):
        t = Transcript()
        t.chrom, t.start, t.end, t.strand, t.id = "Chr1", 101, 1000, "+", "foo.1"
        t.add_exons([(101, 1000)])
        t.finalize()
        for invalid in (10, -0.1, "a"):
            with self.subTest(invalid=invalid):
                with self.assertRaises(TypeError):
                    t.proportion_verified_introns_inlocus = invalid
                t.proportion_verified_introns_inlocus = 1
                with self.assertRaises(TypeError):
                    t.cds_disrupted_by_ri = invalid
                t.cds_disrupted_by_ri = True
                with self.assertRaises(TypeError):
                    t.selected_cds_intron_fraction = invalid
                with self.assertRaises(TypeError):
                    t.combined_cds_intron_fraction = invalid
                with self.assertRaises(TypeError):
                    t.selected_cds_locus_fraction = invalid
                with self.assertRaises(ValueError):
                    # t is non-coding
                    t.selected_cds_locus_fraction = .4

    def test_set_codons(self):
        t = Transcript()
        t.chrom, t.start, t.end, t.strand, t.id = "Chr1", 101, 1000, "+", "foo.1"
        t.add_exons([(101, 1000)])
        t.add_exons([(101, 1000)], features=["CDS"])
        t.finalize()
        for valid in (0, "0", "false", False):
            t.has_start_codon = valid
            t.has_stop_codon = valid
            self.assertFalse(t.has_start_codon)
            self.assertFalse(t.has_stop_codon)
        for valid in (1, "1", "true", "True"):
            t.has_start_codon = valid
            t.has_stop_codon = valid
            self.assertTrue(t.has_start_codon)
            self.assertTrue(t.has_stop_codon)

    def test_conversion(self):

        asbed = self.t1.format("bed12")
        t1 = Transcript(BED12(asbed))
        t1.finalize()
        self.assertEqual(t1, self.t1)

        asgtf = self.t1.format("gtf")
        asgtf = [GtfLine(_) for _ in asgtf.split("\n")]
        t1 = Transcript(asgtf[0])
        t1.add_exons(asgtf[1:])
        t1.finalize()
        self.assertEqual(t1, self.t1)

        asgtf = self.t1.format("gff3")
        asgtf = [GffLine(_) for _ in asgtf.split("\n")]
        t1 = Transcript(asgtf[0])
        t1.add_exons(asgtf[1:])
        t1.finalize()
        self.assertEqual(t1, self.t1)

    def test_get_coding_exons(self):
        correct = sorted([(5697,5891), (5256, 5576)])
        self.t1.finalize()
        self.assertIsNotNone(self.t1.coding_exons)
        self.assertEqual(sorted(self.t1.coding_exons), sorted(correct))
        t1 = self.t1.copy()
        t1.strip_cds()
        self.assertEqual([], t1.coding_exons)
        t1.unfinalize()
        t1.start = 4256
        t1.add_exon((4256, 5000))
        t1.finalize()
        self.assertEqual(sorted(correct), sorted(self.t1.coding_exons))

    def test_wrong_exons(self):
        with self.assertRaises(NotImplementedError):
            self.t1.exons = self.t1.exons
        self.t1.unfinalize()
        with self.assertRaises(TypeError):
            self.t1.exons = [[5697,5891], [5256,5576]]
        with self.assertRaises(TypeError):
            self.t1.exons = [[5697,5891], ["5256",5576]]

    def test_comparisons(self):

        t1 = self.t1.deepcopy()
        t1.unfinalize()
        t1.start = 4256
        t1.add_exon((4256, 5000))
        t1.finalize()

        self.assertEqual(t1.end, self.t1.end)
        self.assertLess(t1.start, self.t1.start)
        self.assertLess(t1, self.t1)
        self.assertEqual(t1.combined_cds_start, self.t1.combined_cds_start)
        self.assertGreater(self.t1, t1)
        self.assertEqual(t1, t1)
        self.assertEqual(self.t1, self.t1)

    def test_to_and_from_dict(self):

        d = self.t1.as_dict()
        self.assertIsInstance(d, dict)
        self.assertIn("chrom", d)
        t1 = Transcript()
        t1.load_dict(d)
        self.assertIs(t1.is_coding, True)
        self.assertEqual(t1, self.t1)

    def test_to_and_from_dict_truncated(self):

        gtf = GtfLine("""Chr5\tfoo\tCDS\t100\t200\t.\t+\t2\ttranscript_id "test.1"; gene_id "test.2";""")
        logger = create_default_logger("test_to_and_from_dict_truncated", level="WARNING")
        t = Transcript(gtf)
        t.finalize()
        self.assertIs(t.is_coding, True)
        self.assertEqual(t.selected_internal_orf, [("CDS", (100, 200), 2), ("exon", (100, 200))])

        d = t.as_dict()
        self.assertEqual(d["orfs"]["0"], [["CDS", [100, 200], 2], ["exon", [100, 200]]])

        new_t = Transcript(logger=logger)
        new_t.load_dict(d)
        self.assertEqual(t.selected_internal_orf, new_t.selected_internal_orf)

    def test_to_and_from_dict_multiple_truncated(self):

        gtf = GtfLine("""Chr5\texon\texon\t1\t1200\t.\t+\t.\ttranscript_id "test.1"; gene_id "test.2";""")
        logger = create_default_logger("tdmt", level="WARNING")
        t = Transcript(gtf, logger=logger)
        t.finalize()
        orf = "test.1\t0\t1200\tID=test.1.orf1;coding=True;phase=2\t0\t+\t0\t320\t.\t1\t1200\t0"
        self.assertEqual(len(orf.split("\t")), 12)
        b = BED12(orf, transcriptomic=True, coding=True)
        self.assertFalse(b.header)
        self.assertFalse(b.invalid, b.invalid_reason)
        orf2 = "test.1\t0\t1200\tID=test.1.orf1;coding=True;phase=0\t0\t+\t400\t1000\t.\t1\t1200\t0"
        b2 = BED12(orf2, transcriptomic=True, coding=True)
        self.assertFalse(b2.header)
        self.assertFalse(b2.invalid, b2.invalid_reason)
        t.load_orfs([b, b2])
        self.assertTrue(t.is_coding)
        self.assertEqual(t.number_internal_orfs, 2)
        with self.subTest():
            new_t = Transcript(logger=logger)
            d = t.as_dict()
            self.assertEqual(len(d["orfs"]), 2)
            new_t.load_dict(d)
            self.assertTrue(new_t.is_coding)
            self.assertEqual(new_t.number_internal_orfs, 2)
            self.assertEqual(new_t.combined_utr, [(321, 400), (1001, 1200)])

        # Now let us test out what happens with three ORFs

        with self.subTest():
            new_t = Transcript(logger=logger)
            d = t.as_dict()
            d["orfs"]["2"] = [("exon", (1, 1200)), ("UTR", (1, 1099)), ("CDS", [1100, 1200], 0)]
            d["finalized"] = False
            new_t.load_dict(d)
            self.assertTrue(new_t.is_coding)
            self.assertEqual(new_t.number_internal_orfs, 3)
            self.assertEqual(new_t.combined_utr, [(321, 400), (1001, 1099)])

        # Now let us check what happens with the addition of an incompatible ORF
        new_t = Transcript()
        d = t.as_dict()
        d["finalized"] = False
        d["orfs"]["2"] = [("CDS", (900, 1200), 0)]
        with self.assertLogs(logger=new_t.logger, level="DEBUG") as cmo:
            new_t.load_dict(d)
        self.assertTrue(new_t.is_coding)
        self.assertEqual(new_t.number_internal_orfs, 2)
        self.assertNotIn(("CDS", (900, 1200), 0), new_t.as_dict()["orfs"].get(2, []))
        s = pd.Series(cmo.output)
        self.assertTrue(any(s.str.contains(
            "ORF 2 of test.1 is invalid. Exception:")), cmo.output)

    def test_as_and_load_dict_verified_introns(self):
        self.t1.verified_introns = {(5577, 5696)}
        self.t1.proportion_verified_introns_inlocus = 1.0
        t1_state = self.t1.as_dict()
        new_t = Transcript()
        new_t.load_dict(t1_state)
        assert new_t.verified_introns == self.t1.verified_introns

    def test_as_and_load_dict_blast_hits(self):
        hit = dict()
        hit["evalue"] = 10**-6
        hit["hsps"] = []
        hit["target"] = "target1"
        hit["target_length"] = 636
        hit["query_start"] = 1
        hit["query_end"] = 636
        hit["target_cov"] = 100
        hit["bits"] = 1200
        hit["query_cov"] = 100
        hit["evalue"] = 10**(-6)
        hit["global_identity"] = 100
        hit["query_aligned_length"] = 636 - 1
        hit["target_aligned_length"] = 636 - 1
        self.t1.blast_hits.append(hit)
        t1_state = self.t1.as_dict()
        new_t = Transcript()
        new_t.load_dict(t1_state)
        assert new_t.blast_hits == self.t1.blast_hits


class TestPicklingAndToFromDict(unittest.TestCase):

    def test_as_dict_and_load(self):
        bed_line = "Chr5\t26584773\t26587782\tID=at_AT5G66610.2;coding=True;phase=0\t0\t+\t26585222\t26587755\t0\t10\t106,54,545,121,78,105,213,63,119,496\t0,446,571,1208,1443,1646,1864,2160,2310,2513"
        original = Transcript(bed_line)
        original.external.foo = 10

        original.finalize()
        recovered = Transcript()
        dumped = original.as_dict()
        recovered.load_dict(dumped, trust_orf=False)

        missed_keys = []
        new_keys = []
        different_keys = []

        for key in original.__dict__:
            if not key in recovered.__dict__:
                missed_keys.append(key)
                continue
            new, old = recovered.__dict__[key], original.__dict__[key]
            if not isinstance(new, type(old)):
                different_keys.append((key, str(new), str(old)))
                continue
            elif old != new:
                if "external" in key:
                    previous = dict((key, value) for key, value in old.items())
                    now = dict((key, value) for key, value in new.items())
                    different = False
                    for key in set.union(set(previous.keys()), set(now.keys())):
                        if key not in previous or key not in now or now[key] != previous[key]:
                            different = True
                            break
                    if different:  # and max(len(now), len(previous)) > 0:
                        different_keys.append((key, str(previous), str(now)))
                elif isinstance(recovered.__dict__[key], IntervalTree):
                    old_ivs = []
                    new_ivs = []
                    old.traverse(lambda iv: old_ivs.append(iv))
                    new.traverse(lambda iv: new_ivs.append(iv))
                    old_ivs = sorted(old_ivs)
                    new_ivs = sorted(new_ivs)
                    if old_ivs != new_ivs:
                        different_keys.append((key, str(old_ivs), str(new_ivs)))
                else:
                    different_keys.append((key, str(old), str(new)))

        for key in recovered.__dict__:
            if not key in original.__dict__:
                new_keys.append(key)

        self.assertEqual(len(missed_keys), 0, "Missed keys: {}".format(", ".join([str(_) for _ in missed_keys])))
        self.assertEqual(len(new_keys), 0, "New keys: {}".format(", ".join([str(_) for _ in new_keys])))
        self.assertEqual(len(different_keys), 0, "New keys: {}".format(", ".join([str(_) for _ in different_keys])))

    def test_pickling_comprehensive(self):
        bed_line = "Chr5\t26584773\t26587782\tID=at_AT5G66610.2;coding=True;phase=0\t0\t+\t26585222\t26587755\t0\t10\t106,54,545,121,78,105,213,63,119,496\t0,446,571,1208,1443,1646,1864,2160,2310,2513"
        original = Transcript(bed_line)
        original.external.foo = 10

        original.finalize()
        recovered = pickle.loads(pickle.dumps(original))
        missed_keys = []
        new_keys = []
        different_keys = []

        for key in original.__dict__:
            if not key in recovered.__dict__:
                if key in ["engine", "session", "sessionmaker"]:
                    # We expect these keys to be missing, they cannot be pickled and we are not trasmitting
                    # the configuration.
                    continue
                missed_keys.append(key)
                continue
            new, old = recovered.__dict__[key], original.__dict__[key]
            if not isinstance(new, type(old)):
                different_keys.append((key, str(new), str(old)))
                continue
            elif old != new:
                if "external" in key:
                    previous = dict((key, value) for key, value in old.items())
                    now = dict((key, value) for key, value in new.items())
                    different = False
                    for key in set.union(set(previous.keys()), set(now.keys())):
                        if key not in previous or key not in now or now[key] != previous[key]:
                            different = True
                            break
                    if different:  # and max(len(now), len(previous)) > 0:
                        different_keys.append((key, str(previous), str(now)))
                elif isinstance(recovered.__dict__[key], IntervalTree):
                    old_ivs = []
                    new_ivs = []
                    old.traverse(lambda iv: old_ivs.append(iv))
                    new.traverse(lambda iv: new_ivs.append(iv))
                    old_ivs = sorted(old_ivs)
                    new_ivs = sorted(new_ivs)
                    if old_ivs != new_ivs:
                        different_keys.append((key, str(old_ivs), str(new_ivs)))
                else:
                    different_keys.append((key, str(old), str(new)))

        for key in recovered.__dict__:
            if not key in original.__dict__:
                new_keys.append(key)

        self.assertEqual(len(missed_keys), 0, "Missed keys: {}".format(", ".join([str(_) for _ in missed_keys])))
        self.assertEqual(len(new_keys), 0, "New keys: {}".format(", ".join([str(_) for _ in new_keys])))
        self.assertEqual(len(different_keys), 0, "New keys: {}".format(", ".join([str(_) for _ in different_keys])))


class FinalizeTests(unittest.TestCase):

    # from Mikado._transcripts.transcript_methods.finalizing import _check_completeness, _check_internal_orf, \
    #     _check_phase_correctness, _calculate_phases, _calculate_introns, _basic_final_checks, _verify_boundaries, \
    #     _check_cdna_vs_utr, _fix_stop_codon

    def test_calculate_introns(self):
        transcript = Transcript()
        transcript.chrom, transcript.start, transcript.end, transcript.strand = "Chr1", 101, 1000, "+"
        transcript.add_exon((101, 1000))
        transcript = _calculate_introns(transcript)
        self.assertEqual(transcript.introns, set())
        self.assertEqual(transcript.splices, set())
        self.assertEqual(transcript.combined_cds_introns, set())
        self.assertEqual(transcript.selected_cds_introns, set())

        transcript = Transcript()
        transcript.chrom, transcript.start, transcript.end, transcript.strand = "Chr1", 101, 1000, "+"
        transcript.add_exons([(101, 500), (601, 1000)])
        transcript = _calculate_introns(transcript)
        self.assertEqual(transcript.introns, {(501, 600)})
        self.assertEqual(transcript.splices, {501, 600})
        self.assertEqual(transcript.combined_cds_introns, set())
        self.assertEqual(transcript.selected_cds_introns, set())

        transcript = Transcript()
        transcript.chrom, transcript.start, transcript.end, transcript.strand = "Chr1", 101, 1000, "+"
        transcript.add_exons([(101, 500), (601, 1000)])
        transcript.add_exons([(201, 500)], features="CDS")
        transcript = _calculate_introns(transcript)
        self.assertEqual(transcript.introns, {(501, 600)})
        self.assertEqual(transcript.splices, {501, 600})
        self.assertEqual(transcript.combined_cds_introns, set())
        self.assertEqual(transcript.selected_cds_introns, set())

        # transcript = Transcript()
        # transcript.chrom, transcript.start, transcript.end, transcript.strand = "Chr1", 101, 1000, "+"
        # transcript.add_exons([(101, 500), (601, 1000)])
        # transcript.add_exons([(201, 500), (601, 900)], features="CDS")
        # self.assertGreater(len(transcript.combined_cds), 0)
        # self.assertGreater(transcript.number_internal_orfs, 0)
        # transcript = _calculate_introns(transcript)
        # self.assertEqual(transcript.introns, {(501, 600)})
        # self.assertEqual(transcript.splices, {501, 600})
        # self.assertEqual(transcript.combined_cds_introns, {(501, 600)}, transcript.combined_cds_introns)
        # self.assertEqual(transcript.selected_cds_introns, {(501, 600)}, transcript.selected_cds_introns)

    def test_basic_final_checks(self):
        transcript = Transcript()
        transcript.start, transcript.end = 101, 1000
        transcript.id = "foo"
        with self.assertRaises(InvalidTranscript) as exc:
            _basic_final_checks(transcript)
        self.assertIsNotNone(re.search(r"No exon defined for the transcript foo. Aborting",
                             str(exc.exception)))

        transcript = Transcript()
        transcript.start, transcript.end = 101, 1000
        transcript.id = "foo"
        transcript._possibly_without_exons = True
        _basic_final_checks(transcript)
        self.assertEqual(transcript.exons, [(101, 1000)])

        transcript = Transcript()
        transcript.start, transcript.end, transcript.id = 101, 1000, "foo"
        transcript.exons = [(91, 1000)]
        _basic_final_checks(transcript)
        self.assertEqual(transcript.start, 91)
        self.assertEqual(transcript.end, 1000)

        transcript = Transcript()
        transcript.start, transcript.end, transcript.id = 101, 1000, "foo"
        transcript.exons = [(101, 5000)]
        _basic_final_checks(transcript)
        self.assertEqual(transcript.start, 101)
        self.assertEqual(transcript.end, 5000)

        # Infer exons from CDS + UTR
        transcript = Transcript()
        transcript.start, transcript.end, transcript.id, transcript.strand = 101, 1000, "foo", "+"
        utr = [(101, 200), (201, 400), (801, 900), (950, 1000)]
        cds = [(401, 600), (701, 800)]
        exons = [(101, 200), (201, 600), (701, 900), (950, 1000)]
        transcript.add_exons(cds, features="CDS")
        transcript.add_exons(utr, features="UTR")
        self.assertEqual(transcript.combined_utr, utr)
        self.assertEqual(transcript.combined_cds, cds)
        _basic_final_checks(transcript)
        self.assertEqual(transcript.exons, exons)

        # Mangled UTR/CDS
        transcript = Transcript()
        transcript.start, transcript.end, transcript.id, transcript.strand = 101, 1000, "foo", "+"
        transcript.add_exons([(201, 500), (601, 900)], features="CDS")
        transcript.add_exons([(101, 210), (850, 1000)], features="UTR")
        self.assertEqual(transcript.combined_utr, [(101, 210), (850, 1000)])
        self.assertEqual(transcript.combined_cds, [(201, 500), (601, 900)])
        with self.assertRaises(InvalidTranscript) as exc:
            _basic_final_checks(transcript)
        self.assertIsNotNone(re.search(r"Transcript.* has a mangled CDS/UTR", str(exc.exception)),
                             str(exc.exception))

        # CDS absent, UTR present
        transcript = Transcript()
        transcript.start, transcript.end, transcript.id, transcript.strand = 101, 1000, "foo", "+"
        transcript.add_exons([(101, 1000)])
        transcript.add_exons([(101, 210), (850, 1000)], features="UTR")
        self.assertEqual(transcript.combined_utr, [(101, 210), (850, 1000)])
        with self.assertRaises(InvalidTranscript) as exc:
            _basic_final_checks(transcript)
        self.assertIsNotNone(re.search(r"Transcript.* has defined UTRs but no CDS.*", str(exc.exception)),
                             str(exc.exception))

        # Overlapping exons
        transcript = Transcript()
        transcript.start, transcript.end, transcript.id, transcript.strand = 101, 1000, "foo", "+"
        transcript.add_exons([(101, 500), (499, 1000)])
        with self.assertRaises(InvalidTranscript) as exc:
            _basic_final_checks(transcript)
        self.assertIsNotNone(re.search(r"Overlapping exons found", str(exc.exception)))

    def test_check_with_double_internal_orf(self):
        transcript = Transcript()
        transcript.start, transcript.end, transcript.id, transcript.strand = 101, 1000, "foo", "+"
        transcript.add_exons([(101, 470), (499, 1000)])
        # Now add the internal ORFs
        transcript.internal_orfs.append(
            [('UTR', (101, 150)), ('exon', (101, 470)), ('CDS', (151, 450), 0), ("UTR", (451, 470)),
             ('exon', (499, 1000)), ('UTR', (499, 1000))])
        transcript.internal_orfs.append(
            [('UTR', (101, 470)), ('exon', (101, 470)), ('UTR', (499, 500)), ('exon', (499, 1000)),
             ('CDS', (501, 560), 0), ('UTR', (561, 1000))])
        with self.assertLogs(transcript.logger, level="DEBUG") as cmo:
            transcript.finalize()
        self.assertTrue(transcript.is_coding, cmo.output)
        self.assertEqual(transcript.combined_cds, [(151, 450), (501, 560)])


if __name__ == '__main__':
    unittest.main()
