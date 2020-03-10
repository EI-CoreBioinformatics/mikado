import os
import unittest
import builtins
from sqlalchemy.engine import reflection
import itertools
from ..configuration.configurator import to_json
from ..loci import Transcript
from ..parsers.bed12 import BED12
from ..parsers.GTF import GtfLine
from ..parsers.GFF import GffLine
from ..transcripts.transcript import Metric
from ..transcripts.transcript_methods import retrieval
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
        gtf_cds = [GtfLine(_) for _ in gtf.split("\n") if GffLine(_).feature == "CDS"]
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
        conf = to_json(os.path.join(
            os.path.dirname(__file__),
            "configuration.yaml"
        ))
        self.assertTrue(conf["pick"]["chimera_split"]["blast_check"])
        self.assertTrue(conf["pick"]["chimera_split"]["execute"])
        self.assertEqual(conf["pick"]["chimera_split"]["blast_params"]["leniency"],
                         "LENIENT")

        conf["pick"]["orf_loading"]["minimal_secondary_orf_length"] = 50

        self.tr.json_conf = conf

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

        self.assertEqual(len(after_overlap_check), 2, self.tr.json_conf["pick"]["orf_loading"])
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
        import pandas as pd
        s = pd.Series(cmo.output)
        self.assertTrue(any(s.str.contains(
            "ORF 2 of test.1 is invalid, removing.*")), cmo.output)


if __name__ == '__main__':
    unittest.main()
