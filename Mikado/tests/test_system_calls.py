import csv
import glob
import gzip
import itertools
import logging
import os
import random
import sys
import tempfile
import unittest
import pkg_resources
import pyfaidx
import yaml
import shutil
import Mikado.daijin
import Mikado.subprograms.configure
from Mikado.configuration import configurator, daijin_configurator
from Mikado.picking import picker
from Mikado.preparation import prepare
from Mikado.scales.compare import compare, load_index
from Mikado.subprograms.util.stats import Calculator
from Mikado.subprograms.prepare import prepare_launcher
from Mikado.subprograms.prepare import setup as prepare_setup
from Mikado.transcripts.transcript import Namespace
from Mikado.utilities.log_utils import create_null_logger, create_default_logger
from Mikado.parsers.GFF import GffLine
from Mikado.parsers.GTF import GtfLine
from Mikado.parsers import to_gff
from Mikado.transcripts import Transcript


class PrepareCheck(unittest.TestCase):

    __genomefile__ = None

    @classmethod
    def setUpClass(cls):
        cls.__genomefile__ = tempfile.NamedTemporaryFile(mode="wb", delete=False, suffix=".fa", prefix="prepare")

        with pkg_resources.resource_stream("Mikado.tests", "chr5.fas.gz") as _:
            cls.__genomefile__.write(gzip.decompress(_.read()))
        cls.__genomefile__.flush()

        cls.trinity_res = dict((_[0], _[1]) for _ in [("tr_c73_g1_i1.mrna1.160", 286),
                                                      ("tr_c11_g1_i2.mrna1.111", 844),
                                                      ("tr_c11_g1_i1.mrna1.350", 659),
                                                      ("tr_c3_g1_i1.mrna1.716", 292),
                                                      ("tr_c3_g1_i2.mrna1.181", 213),
                                                      ("tr_c58_g1_i3.mrna1.19", 1718),
                                                      ("tr_c58_g1_i2.mrna1.35", 900),
                                                      ("tr_c58_g1_i7.mrna1.1", 504),
                                                      ("tr_c58_g1_i1.mrna1.190", 417),
                                                      ("tr_c108_g1_i1.mrna1.104", 458),
                                                      ("tr_c58_g1_i11.mrna1", 1744),
                                                      ("tr_c58_g1_i5.mrna1.3", 4707),
                                                      ("tr_c58_g1_i10.mrna1", 1655),
                                                      ("tr_c68_g1_i1.mrna1.173", 275),
                                                      ("tr_c58_g1_i11.mrna2", 404),
                                                      ("tr_c58_g1_i10.mrna2", 383),
                                                      ("tr_c58_g1_i8.mrna2", 383),
                                                      ("tr_c58_g1_i12.mrna1", 1725),
                                                      ("tr_c58_g1_i9.mrna1", 1988),
                                                      ("tr_c58_g1_i4.mrna1.6", 4862),
                                                      ("tr_c58_g1_i6.mrna1.3", 4945),
                                                      ("tr_c58_g1_i8.mrna1", 1819),
                                                      ("tr_c113_g1_i1.mrna1.94", 379),
                                                      ("tr_c77_g1_i1.mrna1.153", 635),
                                                      ("tr_c109_g1_i1.mrna1.102", 310),
                                                      ("tr_c60_g1_i1.mrna1.189", 1916),
                                                      ("tr_c115_g1_i1.mrna1.88", 243),
                                                      ("tr_c152_g1_i1.mrna1.66", 550),
                                                      ("tr_c21_g1_i1.mrna1.302", 475),
                                                      ("tr_c41_g1_i1.mrna1.224", 608),
                                                      ("tr_c6_g1_i1.mrna1.412", 1295),
                                                      ("tr_c74_g1_i1.mrna1.154", 320),
                                                      ("tr_c71_g1_i1.mrna1.167", 203),
                                                      ("tr_c114_g1_i1.mrna1.89", 206),
                                                      ("tr_c137_g1_i1.mrna1.77", 743),
                                                      ("tr_c37_g1_i2.mrna1.66", 417),
                                                      ("tr_c37_g1_i1.mrna1.234", 449),
                                                      ("tr_c120_g1_i1.mrna1.89", 269)])
        # cls.trinity_res = sorted(cls.trinity_res)

        cls.cuff_results = {"cl_cufflinks_star_at.23553.1": 1735,
                            "cl_cufflinks_star_at.23551.1": 851,
                            "cl_cufflinks_star_at.23551.2": 608,
                            "cl_cufflinks_star_at.23555.1": 1990,
                            "cl_cufflinks_star_at.23555.2": 1902,
                            "cl_cufflinks_star_at.23555.3": 1798,
                            "cl_cufflinks_star_at.23555.4": 688,
                            "cl_cufflinks_star_at.23563.3": 2326,
                            "cl_cufflinks_star_at.23563.2": 2423,
                            "cl_cufflinks_star_at.23563.1": 2418,
                            "cl_cufflinks_star_at.23563.4": 2285,
                            "cl_cufflinks_star_at.23556.1": 1669,
                            "cl_cufflinks_star_at.23557.1": 1410,
                            "cl_cufflinks_star_at.23558.1": 1114,
                            "cl_cufflinks_star_at.23559.1": 323,
                            "cl_cufflinks_star_at.23560.1": 1178,
                            "cl_cufflinks_star_at.23561.1": 504,
                            "cl_cufflinks_star_at.23562.1": 1302,
                            "cl_cufflinks_star_at.23562.2": 1045}

        cls.maxDiff = None

    @classmethod
    def tearDownClass(cls):
        """"""

        cls.__genomefile__.close()
        os.remove(cls.__genomefile__.name)
        if os.path.exists("{}.fai".format(cls.__genomefile__.name)):
            os.remove("{}.fai".format(cls.__genomefile__.name))

    def setUp(self):

        self.conf = configurator.to_json(None)
        self.conf["reference"]["genome"] = self.__genomefile__.name
        self.logger = create_null_logger("prepare")
        self.conf["prepare"]["keep_redundant"] = True

    def tearDown(self):
        logging.shutdown()
        for fname in ["mikado_prepared.fasta", "mikado_prepared.gtf"]:
            fname = os.path.join(tempfile.gettempdir(), fname)
            if os.path.exists(fname):
                os.remove(fname)

    def test_prepare_trinity_gff(self):

        self.conf["prepare"]["files"]["labels"].append("tr")
        self.conf["prepare"]["files"]["output_dir"] = tempfile.gettempdir()
        args = Namespace()
        args.json_conf = self.conf

        for test_file in ("trinity.gff3",
                          "trinity.match_matchpart.gff3",
                          "trinity.cDNA_match.gff3",
                          "trinity.gtf"):
            with self.subTest(test_file=test_file):
                self.conf["prepare"]["files"]["gff"] = [pkg_resources.resource_filename("Mikado.tests",
                                                                                        test_file)]

                prepare.prepare(args, self.logger)

                # Now that the program has run, let's check the output
                fasta = os.path.join(self.conf["prepare"]["files"]["output_dir"],
                                                "mikado_prepared.fasta")
                self.assertGreater(os.stat(fasta).st_size, 0, test_file)
                fa = pyfaidx.Fasta(fasta)
                res = dict((_, len(fa[_])) for _ in fa.keys())
                fa.close()
                self.assertEqual(res, self.trinity_res)
                os.remove(os.path.join(self.conf["prepare"]["files"]["output_dir"],
                                       "mikado_prepared.fasta.fai"))

    def test_prepare_trinity_and_cufflinks(self):

        self.conf["prepare"]["files"]["labels"] = ["cl", "tr"]
        self.conf["prepare"]["files"]["gff"].append(pkg_resources.resource_filename("Mikado.tests",
                                                                                    "cufflinks.gtf"))
        self.conf["prepare"]["files"]["gff"].append("")
        self.conf["prepare"]["files"]["output_dir"] = tempfile.gettempdir()
        self.conf["prepare"]["files"]["out_fasta"] = "mikado_prepared.fasta"
        self.conf["prepare"]["files"]["out"] = "mikado_prepared.gtf"
        args = Namespace()

        for test_file in ("trinity.gff3",
                          "trinity.match_matchpart.gff3",
                          "trinity.cDNA_match.gff3",
                          "trinity.gtf"):
            with self.subTest(test_file=test_file):
                self.conf["prepare"]["files"]["gff"][1] = pkg_resources.resource_filename("Mikado.tests",
                                                                                          test_file)
                self.conf["prepare"]["files"]["out_fasta"] = "mikado_prepared.fasta"
                self.conf["prepare"]["files"]["out"] = "mikado_prepared.gtf"
                args.strip_cds = True

                args.json_conf = self.conf
                prepare.prepare(args, self.logger)

                # Now that the program has run, let's check the output
                self.assertTrue(os.path.exists(os.path.join(self.conf["prepare"]["files"]["output_dir"],
                                                            "mikado_prepared.fasta")))
                self.assertGreater(os.stat(os.path.join(self.conf["prepare"]["files"]["output_dir"],
                                                        "mikado_prepared.fasta")).st_size, 0)

                fa = pyfaidx.Fasta(os.path.join(self.conf["prepare"]["files"]["output_dir"],
                                                "mikado_prepared.fasta"))
                res = dict((_, len(fa[_])) for _ in fa.keys())
                fa.close()
                precal = self.trinity_res.copy()
                precal.update(self.cuff_results)
                self.assertEqual(res, precal)
                os.remove(os.path.join(self.conf["prepare"]["files"]["output_dir"],
                                       "mikado_prepared.fasta.fai"))

    def test_prepare_with_cds(self):

        rev_strand = {"+": "-", "-": "+"}

        self.conf["prepare"]["files"]["labels"] = ["ann"]
        ann_gff3 = pkg_resources.resource_filename("Mikado.tests", "annotation.gff3")
        rev_ann_gff3 = tempfile.NamedTemporaryFile(suffix=".gff3", mode="wt")
        with open(ann_gff3) as ann:
            for line in ann:
                line = GffLine(line)
                if line.header is True:
                    continue
                line.strand = rev_strand[line.strand]  # Invert strand.
                print(line, file=rev_ann_gff3)
        rev_ann_gff3.flush()

        self.conf["prepare"]["files"]["gff"] = []
        self.conf["prepare"]["files"]["output_dir"] = tempfile.gettempdir()
        self.conf["prepare"]["files"]["out_fasta"] = "mikado_prepared.fasta"
        self.conf["prepare"]["files"]["out"] = "mikado_prepared.gtf"
        args = Namespace()
        args.json_conf = self.conf

        for fname in [ann_gff3, rev_ann_gff3.name]:
            for strip in (True, False):
                with self.subTest(fname=fname, strip=strip):
                    self.conf["prepare"]["files"]["gff"] = [fname]
                    args.json_conf["prepare"]["strip_cds"] = strip
                    prepare.prepare(args, logger=self.logger)
                    fasta = os.path.join(self.conf["prepare"]["files"]["output_dir"], "mikado_prepared.fasta")
                    self.assertTrue(os.path.exists(fasta))
                    if strip is True or (strip is False and fname == ann_gff3):
                        self.assertGreater(os.stat(fasta).st_size, 0)
                        fa = pyfaidx.Fasta(fasta)
                        self.assertEqual(len(fa.keys()), 2)
                        fa.close()
                    else:
                        self.assertEqual(os.stat(fasta).st_size, 0)
                    # Now verify that no model has CDS
                    gtf = os.path.join(self.conf["prepare"]["files"]["output_dir"], "mikado_prepared.gtf")
                    models = dict()
                    with to_gff(gtf) as file_gtf:
                        for line in file_gtf:
                            if line.header:
                                continue
                            elif line.is_transcript:
                                models[line.id] = Transcript(line)
                            else:
                                models[line.parent[0]].add_exon(line)
                    [models[model].finalize() for model in models]
                    for model in models:
                        if strip is False:
                            self.assertTrue(models[model].is_coding, models[model].format("gtf"))
                        else:
                            self.assertFalse(models[model].is_coding, models[model].format("gtf"))

    def test_cdna_redundant_cds_not(self):
        """This test will verify whether the new behaviour of not considering redundant two models with same
        exon structure but different CDS does function properly."""

        gtf = pkg_resources.resource_filename("Mikado.tests", "cds_test_1.gtf")
        self.conf["prepare"]["files"]["gff"] = [gtf]
        self.conf["prepare"]["files"]["labels"] = [""]
        self.conf["prepare"]["files"]["out_fasta"] = "mikado_prepared.fasta"
        self.conf["prepare"]["files"]["out"] = "mikado_prepared.gtf"
        self.conf["prepare"]["files"]["log"] = "prepare.log"
        self.conf["prepare"]["strip_cds"] = False

        args = Namespace(default=None)
        args.strip_cds = False
        args.json_conf = self.conf.copy()
        del args.json_conf["prepare"]["files"]["output_dir"]
        args.log = None
        for b in (False, True):
            with self.subTest(b=b):
                folder = tempfile.mktemp()
                os.makedirs(folder)
                # _ = open(os.path.join(folder, args.json_conf["prepare"]["files"]["log"]), "wt")
                args.output_dir = folder
                args.list = None
                args.gffs = None
                args.strand_specific_assemblies = None
                args.labels = None
                args.json_conf = self.conf
                args.keep_redundant = b
                args.out, args.out_fasta = None, None
                args.json_conf["prepare"]["files"]["log"] = "prepare.log"
                args.log = None
                self.logger.setLevel("DEBUG")
                args, _ = prepare_setup(args)
                self.assertEqual(args.output_dir, folder)
                self.assertEqual(args.json_conf["prepare"]["files"]["output_dir"], folder)
                self.assertIn(os.path.dirname(args.json_conf["prepare"]["files"]["out_fasta"]),
                              (folder, ""), args.json_conf)
                self.assertIn(os.path.dirname(args.json_conf["prepare"]["files"]["out"]),
                              (folder, ""), args.json_conf)

                with self.assertRaises(SystemExit) as exi:
                    prepare_launcher(args)
                self.assertTrue(os.path.exists(folder))
                self.assertTrue(os.path.isdir(folder))
                self.assertEqual(exi.exception.code, 0)
                self.assertTrue(os.path.exists(os.path.join(folder,
                                                            "mikado_prepared.fasta")),
                                open(os.path.join(folder, "prepare.log")).read())
                fa = pyfaidx.Fasta(os.path.join(folder,
                                                "mikado_prepared.fasta"))

                if b is True:
                    self.assertEqual(len(fa.keys()), 5)
                    self.assertEqual(sorted(fa.keys()), sorted(["AT5G01530."+str(_) for _ in range(5)]))
                else:
                    self.assertEqual(len(fa.keys()), 4)
                    self.assertIn("AT5G01530.0", fa.keys())
                    self.assertTrue("AT5G01530.1" in fa.keys() or "AT5G01530.2" in fa.keys())
                    self.assertIn("AT5G01530.3", fa.keys())
                    self.assertIn("AT5G01530.4", fa.keys())
                gtf_file = os.path.join(folder, "mikado_prepared.gtf")
                fa.close()
                coding_count = 0
                with to_gff(gtf_file) as gtf:
                    lines = [line for line in gtf]
                    transcripts = dict()
                    for line in lines:
                        if line.is_transcript:
                            transcript = Transcript(line)
                            transcripts[transcript.id] = transcript
                        elif line.is_exon:
                            transcripts[line.transcript].add_exon(line)
                    [transcripts[_].finalize() for _ in transcripts]
                    for transcript in transcripts.values():
                        if transcript.is_coding:
                            coding_count += 1
                            self.assertIn("has_start_codon", transcript.attributes, str(transcript.format("gtf")))
                            self.assertIn("has_stop_codon", transcript.attributes, str(transcript.format("gtf")))
                            self.assertEqual(transcript.attributes["has_start_codon"],
                                             transcript.has_start_codon,
                                             (transcript.id,
                                              transcript.attributes["has_start_codon"],
                                              transcript.has_start_codon))
                            self.assertEqual(transcript.attributes["has_stop_codon"],
                                             transcript.has_stop_codon,
                                             (transcript.id, transcript.attributes["has_stop_codon"],
                                             transcript.has_stop_codon))
                            self.assertEqual(transcript.is_complete,
                                             transcript.has_start_codon and transcript.has_stop_codon)
                    self.assertIn("AT5G01530.3", transcripts)
                    a5 = transcripts["AT5G01530.3"]
                    self.assertTrue(a5.is_coding)
                    self.assertIn("has_start_codon", a5.attributes)
                    self.assertIn("has_stop_codon", a5.attributes)
                    self.assertTrue(a5.has_start_codon)
                    self.assertTrue(a5.has_stop_codon)
                    self.assertTrue(a5.is_complete)

                self.assertGreater(coding_count, 0)
                shutil.rmtree(folder)

    def test_negative_cdna_redundant_cds_not(self):
        """This test will verify whether the new behaviour of not considering redundant two models with same
        exon structure but different CDS does function properly."""

        gtf = pkg_resources.resource_filename("Mikado.tests", "cds_test_2.gtf")
        self.conf["prepare"]["files"]["gff"] = [gtf]
        self.conf["prepare"]["files"]["labels"] = [""]
        self.conf["prepare"]["files"]["out_fasta"] = "mikado_prepared.fasta"
        self.conf["prepare"]["files"]["out"] = "mikado_prepared.gtf"
        self.conf["prepare"]["files"]["log"] = "prepare.log"
        self.conf["prepare"]["strip_cds"] = False
        self.conf["prepare"]["minimum_length"] = 150  # Necessary for testing A5

        args = Namespace()
        args.strip_cds = False
        args.json_conf = self.conf
        for b in (False, ):
            with self.subTest(b=b):
                folder = tempfile.mktemp()
                os.makedirs(folder)
                args.json_conf = self.conf
                args.keep_redundant = b
                args.output_dir = folder
                args.log = None
                args.gff = None
                args.list = None
                args.strand_specific_assemblies = None
                args, _ = prepare_setup(args)
                prepare.prepare(args, self.logger)
                self.assertTrue(os.path.exists(os.path.join(self.conf["prepare"]["files"]["output_dir"],
                                                            "mikado_prepared.fasta")))
                fa = pyfaidx.Fasta(os.path.join(self.conf["prepare"]["files"]["output_dir"],
                                                "mikado_prepared.fasta"))
                if b is True:
                    self.assertEqual(len(fa.keys()), 6)
                    self.assertEqual(sorted(fa.keys()), sorted(["AT5G01015." + str(_) for _ in range(6)]))
                else:
                    self.assertEqual(len(fa.keys()), 5)
                    self.assertIn("AT5G01015.0", fa.keys())
                    self.assertTrue("AT5G01015.1" in fa.keys() or "AT5G01015.2" in fa.keys())
                    self.assertIn("AT5G01015.3", fa.keys())
                    self.assertIn("AT5G01015.4", fa.keys())

                gtf_file = os.path.join(self.conf["prepare"]["files"]["output_dir"], "mikado_prepared.gtf")

                coding_count = 0
                with to_gff(gtf_file) as gtf:
                    lines = [line for line in gtf]
                    transcripts = dict()
                    for line in lines:
                        if line.is_transcript:
                            transcript = Transcript(line)
                            transcripts[transcript.id] = transcript
                        elif line.is_exon:
                            transcripts[line.transcript].add_exon(line)
                    [transcripts[_].finalize() for _ in transcripts]
                    for transcript in transcripts.values():
                        if transcript.is_coding:
                            coding_count += 1
                            self.assertIn("has_start_codon", transcript.attributes, str(transcript.format("gtf")))
                            self.assertIn("has_stop_codon", transcript.attributes, str(transcript.format("gtf")))
                            self.assertEqual(transcript.attributes["has_start_codon"],
                                             transcript.has_start_codon,
                                             (transcript.id,
                                              transcript.attributes["has_start_codon"],
                                              transcript.has_start_codon))
                            self.assertEqual(transcript.attributes["has_stop_codon"],
                                             transcript.has_stop_codon,
                                             (transcript.id, transcript.attributes["has_stop_codon"],
                                              transcript.has_stop_codon))
                            self.assertEqual(transcript.is_complete,
                                             transcript.has_start_codon and transcript.has_stop_codon)

                    if "AT5G01015.1" in transcripts:
                        a_first = transcripts["AT5G01015.1"]
                    else:
                        a_first = transcripts["AT5G01015.2"]
                    self.assertTrue(a_first.is_coding)
                    self.assertIn("has_start_codon", a_first.attributes)
                    self.assertIn("has_stop_codon", a_first.attributes)
                    self.assertTrue(a_first.has_start_codon)
                    self.assertTrue(a_first.has_stop_codon)
                    self.assertTrue(a_first.is_complete)

                    a0 = transcripts["AT5G01015.0"]
                    self.assertFalse(a0.has_start_codon, a0.combined_cds_start)
                    self.assertFalse(a0.has_stop_codon)

                    a4 = transcripts["AT5G01015.4"]
                    self.assertEqual(a4.canonical_intron_proportion, 1)
                    self.assertTrue(a4.has_start_codon,
                                    "\n".join([str(line) for line in lines if line.transcript == "AT5G01015.4"]))
                    self.assertFalse(a4.attributes["has_stop_codon"],
                                     a4.attributes)
                    self.assertFalse(a4.has_stop_codon,
                                     (a4.attributes["has_stop_codon"],
                                      "\n".join([str(line) for line in lines if line.transcript == "AT5G01015.4"])))
                                     # "\n".join([str(line) for line in lines if line.transcript == "AT5G01015.4"]))
                    
                    a5 = transcripts["AT5G01015.5"]
                    self.assertFalse(a5.has_start_codon,
                                     "\n".join([str(line) for line in lines if line.transcript == "AT5G01015.5"]))
                    self.assertTrue(a5.has_stop_codon,
                                    "\n".join([str(line) for line in lines if line.transcript == "AT5G01015.5"]))

                self.assertGreater(coding_count, 0)
                fa.close()
                shutil.rmtree(folder)

    def test_truncated_cds(self):
        files = ["test_truncated_cds.gff3"]
        files = [pkg_resources.resource_filename("Mikado.tests", filename) for filename in files]
        self.conf["prepare"]["files"]["gff"] = files
        self.conf["prepare"]["files"]["labels"] = [""]
        self.conf["prepare"]["files"]["output_dir"] = tempfile.gettempdir()
        self.conf["prepare"]["files"]["out_fasta"] = "mikado_prepared.fasta"
        self.conf["prepare"]["files"]["out"] = "mikado_prepared.gtf"
        self.conf["prepare"]["strip_cds"] = False

        self.conf["reference"]["genome"] = pkg_resources.resource_filename("Mikado.tests",
                                                                           "test_truncated_cds.fa")
        args = Namespace()
        args.strip_cds = False
        args.json_conf = self.conf
        prepare.prepare(args, self.logger)
        self.assertTrue(os.path.exists(os.path.join(self.conf["prepare"]["files"]["output_dir"],
                                                    "mikado_prepared.fasta")))
        fa = pyfaidx.Fasta(os.path.join(self.conf["prepare"]["files"]["output_dir"],
                                        "mikado_prepared.fasta"))
        self.assertEqual(len(fa.keys()), 1)
        gtf_file = os.path.join(self.conf["prepare"]["files"]["output_dir"], "mikado_prepared.gtf")
        with to_gff(gtf_file) as gtf_handle:
            lines = [line for line in gtf_handle]
        cds = [_ for _ in lines if _.feature == "CDS"]
        self.assertEqual(len(cds), 1)
        self.assertEqual(cds[0].frame, 1)
        self.assertEqual(cds[0].phase, 2)
        fa.close()

    def test_source_selection(self):

        # Chr5	TAIR10	mRNA	208937	210445	.	+	.	gene_id "AT5G01530"; transcript_id "AT5G01530.0";
        # Chr5	TAIR10	exon	208937	209593	.	+	.	gene_id "AT5G01530"; transcript_id "AT5G01530.0";
        # Chr5	TAIR10	exon	209881	210445	.	+	.	gene_id "AT5G01530"; transcript_id "AT5G01530.0";

        self.conf["prepare"]["files"]["output_dir"] = tempfile.gettempdir()
        self.conf["prepare"]["files"]["out_fasta"] = "mikado_prepared.fasta"
        self.conf["prepare"]["files"]["out"] = "mikado_prepared.gtf"
        self.conf["prepare"]["strip_cds"] = False
        self.conf["prepare"]["keep_redundant"] = False

        self.conf["reference"]["genome"] = self.__genomefile__.name

        t = Transcript()
        t.chrom, t.start, t.end, t.strand = "Chr5", 208937, 210445, "+"
        t.add_exons([(208937, 209593), (209881, 210445)])
        t.id = "file1.1"
        t.parent = "file1"
        t.finalize()

        t2 = t.deepcopy()
        t2.id = "file2.1"
        t2.parent = "file2"

        t_file = tempfile.NamedTemporaryFile(mode="wt", suffix=".gtf")
        t2_file = tempfile.NamedTemporaryFile(mode="wt", suffix=".gtf")
        print(t.format("gtf"), file=t_file)
        t_file.flush()
        print(t2.format("gtf"), file=t2_file)
        t2_file.flush()

        self.conf["prepare"]["files"]["gff"] = [t_file.name, t2_file.name]
        self.conf["prepare"]["files"]["labels"] = ["T1", "T2"]

        rounds = {
            0: [0, 0, "rand"],
            1: [-2, -2, "rand"],
            2: [10, 10, "rand"],
            3: [0, 0, "rand"],
            4: [1, 0, "T1_file1.1"],
            5: [0, 1, "T2_file2.1"],
            6: [0, 1, "T2_file2.1"],
            7: [0, 1, "T2_file2.1"],
            8: [1, 0, "T1_file1.1"],
            9: [10, 9, "T1_file1.1"],
            10: [9, 10, "T2_file2.1"]
        }

        for round in rounds:
            with self.subTest(round=round, msg="Starting round {} ({})".format(round, rounds[round])):
                t1, t2, res = rounds[round]
                self.conf["prepare"]["files"]["source_score"] = {"T1": t1, "T2": t2}
                args = Namespace()
                args.strip_cds = False
                args.json_conf = self.conf
                prepare.prepare(args, self.logger)
                self.assertGreater(os.stat(self.conf["prepare"]["files"]["out_fasta"]).st_size, 0)
                fa = pyfaidx.Fasta(self.conf["prepare"]["files"]["out_fasta"])
                self.assertEqual(len(fa.keys()), 1, round)
                if res != "rand":
                    key = list(fa.keys())[0]
                    self.assertEqual(key, res, round)

    def test_reference_selection(self):

        self.conf["prepare"]["files"]["output_dir"] = outdir = tempfile.gettempdir()
        self.conf["prepare"]["files"]["out_fasta"] = "mikado_prepared.fasta"
        self.conf["prepare"]["files"]["out"] = "mikado_prepared.gtf"
        self.conf["prepare"]["strip_cds"] = False
        self.conf["prepare"]["keep_redundant"] = False

        self.conf["reference"]["genome"] = self.__genomefile__.name

        t = Transcript()
        # This is *key*. Transcript T1 should never be selected, unless we are having a "lenient" analysis.
        # However, by flagging T1 as a "reference" transcript, we should retain it
        self.conf["prepare"]["lenient"] = False
        t.chrom, t.start, t.end, t.strand = "Chr5", 208937, 210445, "-"
        t.add_exons([(208937, 209593), (209881, 210445)])
        t.id = "file1.1"
        t.parent = "file1"
        t.finalize()

        # Likewise. We will also test that even when different scores are applied, we *still* retain all
        # transcripts.
        t2 = t.deepcopy()
        t2.id = "file2.1"
        t2.parent = "file2"

        rounds = {
            # Standard cases. I expect the transcripts to be reversed, and one to be picked
            0: [0, 0, [], "rand", "+"],
            1: [-2, -2, [], "rand", "+"],
            2: [10, 10, [], "rand", "+"],
            # T1 as reference. T1 should be kept, T1 should be on the negative strand.
            3: [0, 0, ["T1"], sorted(["T1_file1.1"]), "-"],
            # Both as reference. Both should be kept, both should be on the negative strand.
            4: [1, 0, ["T1", "T2"],  sorted(["T1_file1.1", "T2_file2.1"]), "-"],
            5: [0, 0, ["T2"], sorted(["T2_file2.1"]), "-"]

        }

        for fformat in ("gff3", "gtf", "bed12"):
            t_file = tempfile.NamedTemporaryFile(mode="wt", suffix=".{}".format(fformat))
            t2_file = tempfile.NamedTemporaryFile(mode="wt", suffix=".{}".format(fformat))
            print(t.format(fformat), file=t_file)
            t_file.flush()
            print(t2.format(fformat), file=t2_file)
            t2_file.flush()

            self.conf["prepare"]["files"]["gff"] = [t_file.name, t2_file.name]
            self.conf["prepare"]["files"]["strand_specific_assemblies"] = [t_file.name, t2_file.name]
            self.conf["prepare"]["files"]["labels"] = ["T1", "T2"]

            for round in rounds:
                for iteration in range(4):  # Repeat each test 4 times
                    with self.subTest(round=round, format=format, iteration=iteration,
                                      msg="Starting round {} ({})".format(round, rounds[round])):
                        t1_score, t2_score, is_ref, res, corr_strand = rounds[round]
                        self.conf["prepare"]["files"]["source_score"] = {"T1": t1_score,
                                                                         "T2": t2_score}
                        self.conf["prepare"]["files"]["reference"] = []
                        for label in is_ref:
                            if label == "T1":
                                self.conf["prepare"]["files"]["reference"].append(t_file.name)
                            elif label == "T2":
                                self.conf["prepare"]["files"]["reference"].append(t2_file.name)

                        args = Namespace()
                        args.strip_cds = False
                        args.json_conf = self.conf
                        with self.assertLogs(self.logger, "INFO") as cm:
                            prepare.prepare(args, self.logger)
                        self.assertGreater(os.stat(self.conf["prepare"]["files"]["out_fasta"]).st_size, 0,
                                           (round, fformat, iteration, "\n".join(cm.output)))
                        fa = pyfaidx.Fasta(os.path.join(outdir, os.path.basename(
                            self.conf["prepare"]["files"]["out_fasta"])))
                        if res != "rand":
                            key = sorted(list(fa.keys()))
                            self.assertEqual(key, res, round)
                        else:
                            self.assertEqual(len(fa.keys()), 1, (round, fa.keys(), res))
                        gtf = os.path.join(outdir, os.path.basename(
                            self.conf["prepare"]["files"]["out"]))
                        strand = [_.strand for _ in to_gff(gtf)]
                        self.assertEqual(len(set(strand)), 1, strand)
                        self.assertEqual(set(strand).pop(), corr_strand,
                                         (round, self.conf["prepare"]["files"]["reference"],
                                          set(strand), corr_strand))


class CompareCheck(unittest.TestCase):

    """Test to check that compare interacts correctly with match, match_part, cDNA_match"""

    def tearDown(self):
        logging.shutdown()
        for suff in ["log", "refmap", "tmap", "stats"]:
            [os.remove(_) for _ in glob.glob(os.path.join(tempfile.gettempdir(),
                                                          "compare_*.{}".format(suff)))]

    def test_index(self):

        # Create the list of files
        files = ["trinity.gtf",
                 "trinity.gff3",
                 "trinity.cDNA_match.gff3",
                 "trinity.match_matchpart.gff3",
                 "trinity.bed12"]
        # files = [pkg_resources.resource_filename("Mikado.tests", filename) for filename in files]

        namespace = Namespace(default=False)
        namespace.distance = 2000
        namespace.index = True
        namespace.prediction = None
        namespace.log = os.path.join(tempfile.gettempdir(), "index.log")
        logger = create_null_logger("null")

        for ref in files:
            with self.subTest(ref=ref):
                temp_ref = os.path.join(tempfile.gettempdir(), ref)
                with pkg_resources.resource_stream("Mikado.tests", ref) as ref_handle,\
                        open(temp_ref, "wb") as out_handle:
                    out_handle.write(ref_handle.read())
                namespace.reference = to_gff(temp_ref)
                compare(namespace)

                self.assertTrue(os.path.exists(namespace.log))
                self.assertTrue(os.path.exists("{}.midx".format(namespace.reference.name)))
                self.assertGreater(os.stat("{}.midx".format(namespace.reference.name)).st_size, 0)
                genes, positions = load_index(namespace, logger)
                self.assertIsInstance(genes, dict)
                self.assertIsInstance(positions, dict)
                self.assertEqual(len(genes), 38)
                os.remove(namespace.reference.name)
                os.remove(namespace.log)
                os.remove("{}.midx".format(namespace.reference.name))
                namespace.reference.close()

    def test_compare_trinity(self):

        # Create the list of files
        files = ["trinity.gtf",
                 "trinity.gff3",
                 "trinity.cDNA_match.gff3",
                 "trinity.match_matchpart.gff3",
                 "trinity.bed12"]
        files = [pkg_resources.resource_filename("Mikado.tests", filename) for filename in files]

        namespace = Namespace(default=False)
        namespace.distance = 2000
        namespace.no_save_index = True

        for ref, pred in itertools.permutations(files, 2):

            with self.subTest(ref=ref, pred=pred):
                namespace.reference = to_gff(ref)
                namespace.prediction = to_gff(pred)
                namespace.log = os.path.join(tempfile.gettempdir(), "compare_{}_{}.log".format(
                    files.index(ref), files.index(pred)))
                namespace.out = os.path.join(tempfile.gettempdir(), "compare_{}_{}".format(
                    files.index(ref), files.index(pred)))
                compare(namespace)
                refmap = "{}.refmap".format(namespace.out)
                tmap = "{}.tmap".format(namespace.out)
                stats = "{}.stats".format(namespace.out)

                self.assertTrue(os.path.exists(namespace.log))
                # with open(log) as log_handle:
                #     log = [_.rstrip() for _ in log_handle]
                for fname in [refmap, stats, tmap]:
                    self.assertTrue(os.path.exists(fname))
                    self.assertGreater(os.stat(fname).st_size, 0)

                with open(refmap) as _:
                    reader = csv.DictReader(_, delimiter="\t")
                    counter = 0
                    for counter, line in enumerate(reader, start=1):
                        ccode = line["ccode"]
                        self.assertIn(ccode,
                                      ("_", "=", "f,_", "f,="),
                                      (ref, pred, line))

                    self.assertEqual(counter, 38)

        for suff in ["log", "refmap", "tmap", "stats"]:
            [os.remove(_) for _ in glob.glob(os.path.join(tempfile.gettempdir(),
                                                          "compare_*.{}".format(suff)))]


class StatCheck(unittest.TestCase):

    """This unit test takes care of verifying that statistics are generated correctly when
    considering four different inputs. Output will be checked against a standard file."""

    def test_stat(self):

        files = ["trinity.gtf",
                 "trinity.gff3",
                 "trinity.cDNA_match.gff3",
                 "trinity.match_matchpart.gff3",
                 "trinity.bed12"]
        files = [pkg_resources.resource_filename("Mikado.tests", filename) for filename in files]

        std_lines = []
        with pkg_resources.resource_stream("Mikado.tests", "trinity_stats.txt") as t_stats:
            for line in t_stats:
                std_lines.append(line.decode().rstrip())

        namespace = Namespace(default=False)
        namespace.tab_stats = None
        for filename in files:
            with self.subTest(filename=filename):
                namespace.gff = to_gff(filename)
                with open(os.path.join(tempfile.gettempdir(),
                                       "{}.txt".format(os.path.basename(filename))), "w") as out:
                    namespace.out = out
                    Calculator(namespace)()
                self.assertGreater(os.stat(out.name).st_size, 0)
                with open(out.name) as out_handle:
                    lines = [_.rstrip() for _ in out_handle]
                self.assertEqual(std_lines, lines)
                os.remove(out.name)
                namespace.gff.close()


class ConfigureCheck(unittest.TestCase):

    """Test for creating configuration files"""

    __genomefile__ = None

    @classmethod
    def setUpClass(cls):
        cls.__genomefile__ = tempfile.NamedTemporaryFile(mode="wb", delete=False, suffix=".fa", prefix="configure")

        with pkg_resources.resource_stream("Mikado.tests", "chr5.fas.gz") as _:
            cls.__genomefile__.write(gzip.decompress(_.read()))
        cls.__genomefile__.flush()

    @classmethod
    def tearDownClass(cls):
        """"""

        cls.__genomefile__.close()
        os.remove(cls.__genomefile__.name)

    def test_mikado_config(self):
        namespace = Namespace(default=False)
        namespace.scoring = None
        namespace.intron_range = None
        namespace.reference = ""
        namespace.external = None
        namespace.mode = ["permissive"]
        namespace.threads = 1
        namespace.blast_targets = []
        namespace.junctions = []
        namespace.new_scoring = None
        out = os.path.join(tempfile.gettempdir(), "configuration.yaml")
        with open(out, "w") as out_handle:
            namespace.out = out_handle
            Mikado.subprograms.configure.create_config(namespace)
        self.assertGreater(os.stat(out).st_size, 0)
        conf = Mikado.configuration.configurator.to_json(out)
        conf = Mikado.configuration.configurator.check_json(conf)
        conf = Mikado.configuration.configurator.check_json(conf)
        self.assertNotIn("asm_methods", conf)
        os.remove(out)

    def test_mikado_config_full(self):
        namespace = Namespace(default=False)
        namespace.scoring = None
        namespace.intron_range = None
        namespace.reference = ""
        namespace.external = None
        namespace.mode = ["permissive"]
        namespace.threads = 1
        namespace.blast_targets = []
        namespace.junctions = []
        namespace.new_scoring = None
        namespace.full = True
        namespace.daijin = False
        out = os.path.join(tempfile.gettempdir(), "configuration.yaml")
        with open(out, "w") as out_handle:
            namespace.out = out_handle
            Mikado.subprograms.configure.create_config(namespace)
        self.assertGreater(os.stat(out).st_size, 0)
        conf = Mikado.configuration.configurator.to_json(out)
        conf = Mikado.configuration.configurator.check_json(conf)
        conf = Mikado.configuration.configurator.check_json(conf)
        self.assertNotIn("asm_methods", conf)
        os.remove(out)

    def test_mikado_config_daijin(self):
        namespace = Namespace(default=False)
        namespace.scoring = None
        namespace.intron_range = None
        namespace.reference = ""
        namespace.external = None
        namespace.mode = ["permissive"]
        namespace.threads = 1
        namespace.blast_targets = []
        namespace.junctions = []
        namespace.new_scoring = None
        namespace.full = True
        namespace.daijin = True
        out = os.path.join(tempfile.gettempdir(), "configuration.yaml")
        with open(out, "w") as out_handle:
            namespace.out = out_handle
            Mikado.subprograms.configure.create_config(namespace)
        self.assertGreater(os.stat(out).st_size, 0)
        conf = Mikado.configuration.configurator.to_json(out)
        conf = Mikado.configuration.configurator.check_json(conf)
        conf = Mikado.configuration.configurator.check_json(conf)
        self.assertIn("asm_methods", conf)
        os.remove(out)

    def test_mikado_config_daijin_set_from_mode(self):
        namespace = Namespace(default=False)
        namespace.scoring = None
        namespace.intron_range = None
        namespace.reference = ""
        namespace.external = None
        namespace.mode = ["permissive", "split"]
        namespace.threads = 1
        namespace.blast_targets = []
        namespace.junctions = []
        namespace.new_scoring = None
        namespace.full = True
        namespace.daijin = False
        out = os.path.join(tempfile.gettempdir(), "configuration.yaml")
        with open(out, "w") as out_handle:
            namespace.out = out_handle
            Mikado.subprograms.configure.create_config(namespace)
        self.assertGreater(os.stat(out).st_size, 0)
        conf = Mikado.configuration.configurator.to_json(out)
        conf = Mikado.configuration.configurator.check_json(conf)
        conf = Mikado.configuration.configurator.check_json(conf)
        self.assertIn("asm_methods", conf)
        os.remove(out)

    @unittest.skipUnless((sys.version_info.minor > 4),
                         "Due to a bug in JSONSCHEMA, Daijin configure fails with Python versions lower than 3.5.")
    def test_daijin_config(self):

        # Check the basic function actually functions
        _ = daijin_configurator.create_daijin_base_config()

        namespace = Namespace(default=False)
        namespace.r1 = []
        namespace.r2 = []
        namespace.samples = []
        namespace.strandedness = []
        namespace.asm_methods = []
        namespace.aligners = []
        namespace.modes = ["nosplit"]
        namespace.cluster_config = None
        namespace.scheduler = ""
        namespace.flank = None
        namespace.intron_range = None
        namespace.prot_db = []
        namespace.genome = self.__genomefile__.name
        namespace.transcriptome = ""
        namespace.name = "Daijin"
        namespace.out_dir = tempfile.gettempdir()
        namespace.threads = 1

        namespace.scoring = random.choice(
            pkg_resources.resource_listdir("Mikado.configuration", "scoring_files"))

        out = os.path.join(tempfile.gettempdir(), "configuration.yaml")
        with open(out, "wt") as out_handle:
            namespace.out = out_handle
            daijin_configurator.create_daijin_config(namespace, level="ERROR")
        self.assertGreater(os.stat(out).st_size, 0)

        with open(out) as out_handle:
            config = yaml.load(out_handle)

        daijin_configurator.check_config(config)


class PickTest(unittest.TestCase):

    """This unit test will check that pick functions correctly."""

    def test_single_proc(self):

        json_conf = configurator.to_json(None)
        json_conf["pick"]["run_options"]["procs"] = 1
        json_conf["db_settings"]["db"] = pkg_resources.resource_filename("Mikado.tests", "mikado.db")

        json_conf["pick"]["files"]["input"] = pkg_resources.resource_filename("Mikado.tests",
                                                                              "mikado_prepared.gtf")
        json_conf["pick"]["files"]["output_dir"] = tempfile.gettempdir()
        json_conf["pick"]["files"]["loci_out"] = "mikado.monoproc.loci.gff3"
        json_conf["pick"]["files"]["subloci_out"] = "mikado.monoproc.subloci.gff3"
        json_conf["pick"]["files"]["monoloci_out"] = "mikado.monoproc.monoloci.gff3"
        json_conf["pick"]["files"]["log"] = "mikado.monoproc.log"
        json_conf["log_settings"]["log_level"] = "WARNING"

        pick_caller = picker.Picker(json_conf=json_conf)
        with self.assertRaises(SystemExit), self.assertLogs("main_logger", "INFO") as cm:
            pick_caller()
        self.assertTrue(os.path.exists(os.path.join(tempfile.gettempdir(), "mikado.monoproc.loci.gff3")))
        with to_gff(os.path.join(tempfile.gettempdir(), "mikado.monoproc.loci.gff3")) as inp_gff:
            lines = [_ for _ in inp_gff if not _.header is True]
            self.assertGreater(len(lines), 0)
            self.assertGreater(len([_ for _ in lines if _.is_transcript is True]), 0)
            self.assertGreater(len([_ for _ in lines if _.feature == "mRNA"]), 0,
                               [_ for _ in cm.output if "WARNING" in _])
            self.assertGreater(len([_ for _ in lines if _.feature == "CDS"]), 0)

        [os.remove(_) for _ in glob.glob(os.path.join(tempfile.gettempdir(), "mikado.monoproc.") + "*")]

    def test_multi_proc(self):
        json_conf = configurator.to_json(None)
        json_conf["pick"]["run_options"]["procs"] = 2
        json_conf["pick"]["files"]["input"] = pkg_resources.resource_filename("Mikado.tests",
                                                                              "mikado_prepared.gtf")
        json_conf["pick"]["files"]["output_dir"] = tempfile.gettempdir()
        json_conf["pick"]["files"]["loci_out"] = "mikado.multiproc.loci.gff3"
        json_conf["pick"]["files"]["subloci_out"] = "mikado.multiproc.subloci.gff3"
        json_conf["pick"]["files"]["monoloci_out"] = "mikado.multiproc.monoloci.gff3"
        json_conf["pick"]["files"]["log"] = "mikado.multiproc.log"
        json_conf["db_settings"]["db"] = pkg_resources.resource_filename("Mikado.tests", "mikado.db")
        json_conf["log_settings"]["log_level"] = "WARNING"

        pick_caller = picker.Picker(json_conf=json_conf)
        with self.assertRaises(SystemExit), self.assertLogs("main_logger", "INFO"):
            pick_caller()
        self.assertTrue(os.path.exists(os.path.join(tempfile.gettempdir(), "mikado.multiproc.loci.gff3")))
        with to_gff(os.path.join(tempfile.gettempdir(), "mikado.multiproc.loci.gff3")) as inp_gff:
            lines = [_ for _ in inp_gff if not _.header is True]
            self.assertGreater(len(lines), 0)
            self.assertGreater(len([_ for _ in lines if _.is_transcript is True]), 0)
            self.assertGreater(len([_ for _ in lines if _.feature == "mRNA"]), 0)
            self.assertGreater(len([_ for _ in lines if _.feature == "CDS"]), 0)

        [os.remove(_) for _ in glob.glob(os.path.join(tempfile.gettempdir(), "mikado.multiproc.") + "*")]

    def test_subprocess(self):
        
        json_conf = configurator.to_json(None)
        
        json_conf["pick"]["files"]["input"] = pkg_resources.resource_filename("Mikado.tests",
                                                                              "mikado_prepared.gtf")
        json_conf["pick"]["files"]["output_dir"] = tempfile.gettempdir()
        json_conf["pick"]["files"]["loci_out"] = "mikado.subproc.loci.gff3"
        json_conf["pick"]["files"]["subloci_out"] = "mikado.subproc.subloci.gff3"
        json_conf["pick"]["files"]["monoloci_out"] = "mikado.subproc.monoloci.gff3"
        json_conf["pick"]["files"]["log"] = "mikado.subproc.log"
        json_conf["db_settings"]["db"] = pkg_resources.resource_filename("Mikado.tests", "mikado.db")
        json_conf["log_settings"]["log_level"] = "WARNING"
        
        for num in (1, 2):
            with self.subTest(num=num):

                json_conf["pick"]["run_options"]["procs"] = num
                json_conf["pick"]["run_options"]["single_thread"] = (num == 1)
                json_file = os.path.join(tempfile.gettempdir(), "mikado.yaml")

                with open(json_file, "wt") as json_handle:
                    Mikado.subprograms.configure.print_config(yaml.dump(json_conf, default_flow_style=False),
                                                              json_handle)

                sys.argv = ["mikado", "pick", "--json-conf", json_file]
                with self.assertRaises(SystemExit):
                    pkg_resources.load_entry_point("Mikado", "console_scripts", "mikado")()
                self.assertTrue(os.path.exists(os.path.join(tempfile.gettempdir(), "mikado.subproc.loci.gff3")))
                with to_gff(os.path.join(tempfile.gettempdir(), "mikado.subproc.loci.gff3")) as inp_gff:
                    lines = [_ for _ in inp_gff if not _.header is True]
                    self.assertGreater(len(lines), 0)
                    self.assertGreater(len([_ for _ in lines if _.is_transcript is True]), 0)
                    self.assertGreater(len([_ for _ in lines if _.feature == "mRNA"]), 0)
                    self.assertGreater(len([_ for _ in lines if _.feature == "CDS"]), 0)

                [os.remove(_) for _ in glob.glob(os.path.join(tempfile.gettempdir(), "mikado.subproc.") + "*")]

    def test_different_scoring(self):

        json_conf = configurator.to_json(None)

        json_conf["pick"]["files"]["input"] = pkg_resources.resource_filename("Mikado.tests",
                                                                              "mikado_prepared.gtf")
        json_conf["pick"]["files"]["output_dir"] = tempfile.gettempdir()
        json_conf["pick"]["files"]["loci_out"] = "mikado.test_diff.loci.gff3"
        json_conf["pick"]["files"]["subloci_out"] = "mikado.test_diff.subloci.gff3"
        json_conf["pick"]["files"]["monoloci_out"] = "mikado.test_diff.monoloci.gff3"
        json_conf["pick"]["files"]["log"] = "mikado.test_diff.log"
        json_conf["db_settings"]["db"] = pkg_resources.resource_filename("Mikado.tests", "mikado.db")
        json_conf["log_settings"]["log_level"] = "WARNING"
        self.assertEqual(os.path.basename(json_conf["pick"]["scoring_file"]),
                         "plants.yaml")
        json_file = os.path.join(tempfile.gettempdir(), "mikado.yaml")
        with open(json_file, "wt") as json_handle:
            Mikado.subprograms.configure.print_config(yaml.dump(json_conf, default_flow_style=False),
                                                      json_handle)

        sys.argv = ["mikado", "pick", "--json-conf", json_file]
        with self.assertRaises(SystemExit):
            pkg_resources.load_entry_point("Mikado", "console_scripts", "mikado")()

        import csv
        with open(os.path.join(json_conf["pick"]["files"]["output_dir"], "mikado.test_diff.loci.scores.tsv")) as tsv:
            reader = csv.DictReader(tsv, delimiter="\t")
            score_names = [_ for _ in json_conf["scoring"]]
            score_header = [_ for _ in reader.fieldnames if _ not in
                            ("tid", "alias", "parent", "score", "source_score")]
            self.assertEqual(sorted(score_names), sorted(score_header))

        scoring_file = pkg_resources.resource_filename("Mikado.tests", "scoring_only_cds.yaml")
        sys.argv = ["mikado", "pick", "--json-conf", json_file, "--scoring-file", scoring_file]
        with self.assertRaises(SystemExit):
            pkg_resources.load_entry_point("Mikado", "console_scripts", "mikado")()
        with open(os.path.join(json_conf["pick"]["files"]["output_dir"], "mikado.test_diff.loci.scores.tsv")) as tsv:
            reader = csv.DictReader(tsv, delimiter="\t")
            score_header = [_ for _ in reader.fieldnames if _ not in
                            ("tid", "alias", "parent", "score", "source_score")]
            self.assertEqual(score_header, ["selected_cds_length"])

        [os.remove(_) for _ in glob.glob(os.path.join(tempfile.gettempdir(), "mikado.test_diff.") + "*")]

    def test_purging(self):

        gtf = """Chr1	foo	transcript	100	1000	.	+	.	gene_id "foo1"; transcript_id "foo1.1"
Chr1	foo	exon	100	1000	.	+	.	gene_id "foo1"; transcript_id "foo1.1"
Chr1	foo	transcript	100	2000	.	+	.	gene_id "foo1"; transcript_id "foo1.2"
Chr1	foo	exon	100	800	.	+	.	gene_id "foo1"; transcript_id "foo1.2"
Chr1	foo	exon	1900	2000	.	+	.	gene_id "foo1"; transcript_id "foo1.2"
Chr1	foo	transcript	10000	20000	.	+	.	gene_id "foo2"; transcript_id "foo2.1"
Chr1	foo	exon	10000	13000	.	+	.	gene_id "foo2; transcript_id "foo2.1"
Chr1	foo	exon	19000	20000	.	+	.	gene_id "foo"; transcript_id "foo2.1"""

        temp_gtf = tempfile.NamedTemporaryFile(mode="wt", suffix=".gtf", delete=True)

        temp_gtf.write(gtf)
        temp_gtf.flush()

        json_conf = configurator.to_json(None)

        json_conf["pick"]["files"]["input"] = temp_gtf.name
        json_conf["db_settings"]["db"] = os.path.join(tempfile.gettempdir(), "mikado.db")
        json_conf["pick"]["files"]["output_dir"] = tempfile.gettempdir()
        json_conf["log_settings"]["log_level"] = "WARNING"

        # Now the scoring
        scoring = dict()
        
        scoring["requirements"] = dict()
        scoring["requirements"]["expression"] = ["exon_num"]
        scoring["requirements"]["parameters"] = dict()
        scoring["requirements"]["parameters"]["exon_num"] = dict()
        scoring["requirements"]["parameters"]["exon_num"]["name"] = "exon_num"
        scoring["requirements"]["parameters"]["exon_num"]["operator"] = "gt"
        scoring["requirements"]["parameters"]["exon_num"]["value"] = 1

        import copy
        scoring["as_requirements"] = copy.deepcopy(scoring["requirements"])
        scoring["not_fragmentary"] = copy.deepcopy(scoring["requirements"].copy())

        scoring["scoring"] = dict()
        scoring["scoring"]["cdna_length"] = dict()
        scoring["scoring"]["cdna_length"]["rescaling"] = "max"
        scoring["scoring"]["cdna_length"]["filter"] = dict()
        scoring["scoring"]["cdna_length"]["filter"]["operator"] = "gt"
        scoring["scoring"]["cdna_length"]["filter"]["value"] = 2000

        scoring_file = tempfile.NamedTemporaryFile(suffix=".yaml", delete=True, mode="wt")
        yaml.dump(scoring, scoring_file)
        scoring_file.flush()
        json_conf["pick"]["scoring_file"] = scoring_file.name
        del json_conf["scoring"]
        del json_conf["requirements"]
        del json_conf["as_requirements"]
        del json_conf["not_fragmentary"]

        for purging in (False, True):
            with self.subTest(purging=purging):
                json_conf["pick"]["files"]["loci_out"] = "mikado.purging_{}.loci.gff3".format(purging)
                json_conf["pick"]["files"]["log"] = "mikado.purging_{}.log".format(purging)
                json_conf["pick"]["clustering"]["purge"] = purging
                json_conf["pick"]["scoring_file"] = scoring_file.name
                json_conf = configurator.check_json(json_conf)
                self.assertEqual(len(json_conf["scoring"].keys()), 1, json_conf["scoring"].keys())

                pick_caller = picker.Picker(json_conf=json_conf)
                with self.assertRaises(SystemExit), self.assertLogs("main_logger", "INFO"):
                    pick_caller()

                with to_gff(os.path.join(tempfile.gettempdir(),
                                                              json_conf["pick"]["files"]["loci_out"])) as gff:

                    lines = [line for line in gff if line.header is False]
                self.assertGreater(len(lines), 0)
                self.assertTrue(any([_ for _ in lines if _.attributes.get("alias", "") == "foo2.1"]))
                if purging is True:
                    self.assertFalse(any([_ for _ in lines if _.attributes.get("alias", "") in ("foo1.2", "foo1.1")]))
                else:
                    found_line = [_ for _ in lines if _.attributes.get("alias", "") in ("foo1.2", "foo1.1")]
                    self.assertTrue(any(found_line))
                    self.assertTrue(any([_ for _ in found_line if _.score == 0]))

            # Clean up
            for fname in ["mikado.db", "mikado.purging_{}.*".format(purging)]:
                [os.remove(_) for _ in glob.glob(os.path.join(tempfile.gettempdir(), fname))]

        scoring_file.close()
        # Now let us test with a scoring which will create transcripts with negative scores
        scoring["scoring"] = dict()
        scoring["scoring"]["cdna_length"] = dict()
        scoring["scoring"]["cdna_length"]["rescaling"] = "min"
        scoring["scoring"]["cdna_length"]["multiplier"] = -10
        scoring["scoring"]["cdna_length"]["filter"] = dict()
        scoring["scoring"]["cdna_length"]["filter"]["operator"] = "lt"
        scoring["scoring"]["cdna_length"]["filter"]["value"] = 1000

        scoring["scoring"]["exon_num"] = dict()
        scoring["scoring"]["exon_num"]["rescaling"] = "max"

        scoring_file = tempfile.NamedTemporaryFile(suffix=".yaml", delete=True, mode="wt")
        yaml.dump(scoring, scoring_file)
        scoring_file.flush()
        json_conf["pick"]["scoring_file"] = scoring_file.name

        for purging in (False, True):
            with self.subTest(purging=purging):
                json_conf["pick"]["files"]["loci_out"] = "mikado.purging_{}.loci.gff3".format(purging)
                json_conf["pick"]["files"]["subloci_out"] = "mikado.purging_{}.subloci.gff3".format(purging)
                json_conf["pick"]["files"]["log"] = os.path.join(
                    tempfile.gettempdir(),
                    "mikado.purging_{}.log".format(purging))
                json_conf["pick"]["clustering"]["purge"] = purging
                json_conf["pick"]["scoring_file"] = scoring_file.name
                json_conf = configurator.check_json(json_conf)
                self.assertEqual(len(json_conf["scoring"].keys()), 2, json_conf["scoring"].keys())

                pick_caller = picker.Picker(json_conf=json_conf)
                with self.assertRaises(SystemExit), self.assertLogs("main_logger", "INFO"):
                    pick_caller()

                with to_gff(os.path.join(tempfile.gettempdir(),
                                                              json_conf["pick"]["files"]["loci_out"])) as gff:
                    lines = [line for line in gff if line.header is False]
                self.assertGreater(len(lines), 0)
                self.assertTrue(any([_ for _ in lines if _.attributes.get("alias", "") == "foo2.1"]))
                if purging is True:
                    self.assertFalse(any([_ for _ in lines if _.attributes.get("alias", "") in ("foo1.2", "foo1.1")]))
                else:
                    found_line = [_ for _ in lines if _.attributes.get("alias", "") in ("foo1.2", "foo1.1")]
                    self.assertTrue(any(found_line))
                    self.assertTrue(any([_ for _ in found_line if _.score <= 0]))

            # Clean up
            for fname in ["mikado.db", "mikado.purging_{}.*".format(purging)]:
                [os.remove(_) for _ in glob.glob(os.path.join(tempfile.gettempdir(), fname))]

        temp_gtf.close()

        temp_gtf = tempfile.NamedTemporaryFile(mode="wt", suffix=".gtf", delete=True)

        gtf = "\n".join([_ for _ in gtf.split("\n") if "foo1.1" not in _])

        temp_gtf.write(gtf)
        temp_gtf.flush()
        json_conf["pick"]["files"]["input"] = temp_gtf.name

        for purging in (False, True):
            with self.subTest(purging=purging):
                json_conf["pick"]["files"]["loci_out"] = "mikado.purging_{}.loci.gff3".format(purging)
                json_conf["pick"]["files"]["subloci_out"] = "mikado.purging_{}.subloci.gff3".format(purging)
                json_conf["pick"]["files"]["log"] = "mikado.purging_{}.log".format(purging)
                json_conf["pick"]["clustering"]["purge"] = purging
                json_conf["pick"]["scoring_file"] = scoring_file.name
                json_conf = configurator.check_json(json_conf)
                self.assertEqual(len(json_conf["scoring"].keys()), 2, json_conf["scoring"].keys())

                pick_caller = picker.Picker(json_conf=json_conf)
                with self.assertRaises(SystemExit), self.assertLogs("main_logger", "INFO"):
                    pick_caller()

                with to_gff(os.path.join(tempfile.gettempdir(),
                                                              json_conf["pick"]["files"]["loci_out"])) as gff:
                    lines = [line for line in gff if line.header is False]
                self.assertGreater(len(lines), 0)
                self.assertTrue(any([_ for _ in lines if _.attributes.get("alias", "") == "foo2.1"]))
                if purging is True:
                    self.assertFalse(any([_ for _ in lines if _.attributes.get("alias", "") == "foo1.2"]))
                else:
                    found_line = [_ for _ in lines if _.attributes.get("alias", "") == "foo1.2"]
                    self.assertTrue(any(found_line))
                    self.assertTrue(any([_ for _ in found_line if _.score <= 0]),
                                    "\n".join([str(_) for _ in found_line]))

            # Clean up
            for fname in ["mikado.db", "mikado.purging_{}.*".format(purging)]:
                [os.remove(_) for _ in glob.glob(os.path.join(tempfile.gettempdir(), fname))]


if __name__ == "__main__":
    unittest.main()
