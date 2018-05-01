import os
import unittest

from sqlalchemy.engine import reflection

from Mikado.configuration.configurator import to_json
from Mikado.loci import Transcript
from Mikado.parsers.bed12 import BED12
from Mikado.transcripts.transcript_methods import retrieval


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

        # print(*cm.output, sep="\n")

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
        self.assertEqual(b12.block_sizes,
                         [200, 200, 400, 500],
                         b12.block_sizes)
        self.assertEqual(b12.strand, "+")
        self.assertEqual(b12.block_starts,
                         [0, 300, 700, 2400],
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
        self.assertEqual(b12.block_sizes,
                         [200, 200, 400, 500],
                         b12.block_sizes)
        self.assertEqual(b12.strand, "-")
        self.assertEqual(b12.block_starts,
                         [0, 300, 700, 2400],
                         b12.block_starts)

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
Chr1\tMikado\tCDS\t1701\t2000\t.\t+\t2\tgene_id "gene1"; transcript_id "test1";
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
        # print(*_.output, sep="\n")

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

if __name__ == '__main__':
    unittest.main()
