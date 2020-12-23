import csv
import glob
import gzip
import itertools
import logging
import os
import numpy as np
import pandas as pd
import sys
import tempfile
import unittest
import pkg_resources
import pyfaidx
import yaml
from ..configuration import print_config
try:
    from yaml import CSafeLoader as yLoader
except ImportError:
    from yaml import SafeLoader as yLoader
from pytest import mark
from .. import configuration
from ..subprograms import configure as sub_configure
from ..configuration import configurator, daijin_configurator
from ..picking import picker
from ..preparation import prepare
from ..scales.compare import compare
from ..scales.reference_preparation.indexing import load_index
from ..scales.calculator import Calculator
from ..subprograms.prepare import prepare_launcher
from ..subprograms.prepare import setup as prepare_setup
from ..transcripts.transcript import Namespace
from ..utilities.log_utils import create_null_logger, create_default_logger
from ..parsers.GFF import GffLine
import sqlite3
import shutil
from ..parsers import to_gff
from ..transcripts import Transcript
import threading
from time import sleep
import pysam
import io


class ConvertCheck(unittest.TestCase):

    @mark.slow
    def test_convert_from_bam(self):

        bam_inp = pkg_resources.resource_filename("Mikado.tests", "test_mRNA.bam")
        for outp in ("gff3", "gtf", "bed12"):
            with self.subTest(outp=outp):
                outfile = tempfile.NamedTemporaryFile(mode="wt")
                outfile.close()
                sys.argv = ["", "util", "convert", "-of", outp, bam_inp, outfile.name]
                # with self.assertRaises(SystemExit):
                pkg_resources.load_entry_point("Mikado", "console_scripts", "mikado")()
                self.assertGreater(os.stat(outfile.name).st_size, 0)
                lines = [_ for _ in open(outfile.name)]
                self.assertTrue(any(["TraesCS2B02G055500.1" in line for line in lines]))

    @mark.slow
    def test_convert_from_problematic(self):
        probl = pkg_resources.resource_filename("Mikado.tests", "Chrysemys_picta_bellii_problematic.gff3")
        for outp in ("gtf", "bed12"):
            with self.subTest(outp=outp):
                outfile = tempfile.NamedTemporaryFile(mode="wt")
                outfile.close()
                sys.argv = ["", "util", "convert", "-of", outp, probl, outfile.name]
                # with self.assertRaises(SystemExit):
                pkg_resources.load_entry_point("Mikado", "console_scripts", "mikado")()
                self.assertGreater(os.stat(outfile.name).st_size, 0)
                lines = [_ for _ in open(outfile.name)]
                self.assertTrue(any(["rna-NC_023890.1:71..1039" in line for line in lines]))
                self.assertTrue(any(["rna-NC_023890.1:1040..1107" in line for line in lines]))
                self.assertTrue(any(["gene-LOC112059550" in line for line in lines]))
                self.assertTrue(any(["id-LOC112059311" in line for line in lines]))


# @mark.slow
class PrepareCheck(unittest.TestCase):

    __genomefile__ = None

    @classmethod
    def setUpClass(cls):
        # cls.__genomefile__ = tempfile.NamedTemporaryFile(mode="wb", delete=False, suffix=".fa.gz", prefix="prepare")
        # cls.__genomefile__.write(pkg_resources.resource_stream("Mikado.tests", "chr5.fas.gz").read())
        # cls.__genomefile__.flush()
        # cls.fai = pyfaidx.Fasta(cls.__genomefile__.name)
        cls.fai = pysam.FastaFile(pkg_resources.resource_filename("Mikado.tests", "chr5.fas.gz"))

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

        cls._trinity_redundant = ["c58_g1_i10.mrna2", "c58_g1_i8.mrna2"]

        cls.maxDiff = None

    def setUp(self):

        self.conf = configurator.to_json(None)
        self.conf["seed"] = 1066
        self.conf["reference"]["genome"] = self.fai.filename.decode()
        assert isinstance(self.conf["reference"]["genome"], str)
        self.logger = create_null_logger("prepare")
        self.conf["prepare"]["exclude_redundant"] = False

    def tearDown(self):
        logging.shutdown()

    @mark.slow
    def test_varying_max_intron(self):

        self.conf["prepare"]["files"]["labels"].append("tr")
        dir = tempfile.TemporaryDirectory(prefix="test_varying_max_intron")
        self.conf["prepare"]["files"]["output_dir"] = dir.name
        args = Namespace()
        args.json_conf = self.conf
        test_file = "trinity.gtf"
        self.conf["prepare"]["files"]["gff"] = [pkg_resources.resource_filename("Mikado.tests",
                                                                                test_file)]
        self.conf["prepare"]["files"]["strip_cds"] = [False]

        for max_intron in (20, 200, 1000, 5000):
            with self.subTest(max_intron=max_intron):
                self.conf["prepare"]["max_intron_length"] = max_intron
                prepare.prepare(args, self.logger)
                gtf = os.path.join(self.conf["prepare"]["files"]["output_dir"], "mikado_prepared.gtf")
                self.assertGreater(os.stat(gtf).st_size, 0, test_file)
                transcripts = dict()
                for row in to_gff(gtf):
                    if row.is_transcript:
                        transcripts[row.transcript] = Transcript(row)
                    else:
                        transcripts[row.transcript].add_exon(row)
                self.assertGreater(len(transcripts), 0)
                [_.finalize() for _ in transcripts.values()]
                self.assertLessEqual(max([_.max_intron_length for _ in transcripts.values()]),
                                     max_intron)
                os.remove(gtf)

    @mark.slow
    def test_prepare_trinity_gff(self):

        self.conf["prepare"]["files"]["labels"].append("tr")
        dir = tempfile.TemporaryDirectory(prefix="test_prepare_trinity_gff")
        self.conf["prepare"]["files"]["output_dir"] = dir.name
        args = Namespace()
        args.json_conf = self.conf
        args.procs = 1
        args.single_thread = True
        args.seed = 10

        for test_file in ("trinity.gff3",
                          "trinity.match_matchpart.gff3",
                          "trinity.cDNA_match.gff3",
                          "trinity.gtf",
                          "trinity.no_transcript_feature.gtf"):
            with self.subTest(test_file=test_file):
                self.conf["prepare"]["files"]["strip_cds"] = [False]
                self.conf["prepare"]["files"]["gff"] = [pkg_resources.resource_filename("Mikado.tests",
                                                                                        test_file)]

                with self.assertLogs(self.logger) as cm:
                    try:
                        prepare.prepare(args, self.logger)
                    except OSError:
                        raise OSError(cm.output)

                # Now that the program has run, let's check the output
                fasta = os.path.join(self.conf["prepare"]["files"]["output_dir"],
                                                "mikado_prepared.fasta")
                self.assertGreater(os.stat(fasta).st_size, 0, (test_file, cm.output))
                fa = pyfaidx.Fasta(fasta)
                res = dict((_, len(fa[_])) for _ in fa.keys())
                fa.close()
                check = dict(_ for _ in self.trinity_res.items())
                check.pop(self.conf["prepare"]["files"]["labels"][0] + "_" + self._trinity_redundant[0])
                self.assertEqual(res, check)
                os.remove(os.path.join(self.conf["prepare"]["files"]["output_dir"],
                                       "mikado_prepared.fasta.fai"))
        dir.cleanup()

    @mark.slow
    def test_prepare_trinity_and_cufflinks(self):

        self.conf["prepare"]["files"]["labels"] = ["cl", "tr"]

        self.conf["prepare"]["files"]["gff"] = [None, None]
        self.conf["prepare"]["files"]["strip_cds"] = [False, True]
        dir = tempfile.TemporaryDirectory(prefix="test_prepare_trinity_and_cufflinks")
        self.conf["prepare"]["files"]["output_dir"] = dir.name
        self.conf["prepare"]["files"]["out_fasta"] = "mikado_prepared.fasta"
        self.conf["prepare"]["files"]["out"] = "mikado_prepared.gtf"
        args = Namespace()

        for cuff_file, test_file in itertools.product(
                ("cufflinks.gtf", "cufflinks.no_transcript.gtf"),
                (("trinity.gff3", "trinity.match_matchpart.gff3", "trinity.cDNA_match.gff3", "trinity.gtf",
                  "trinity.no_transcript_feature.gtf"))):
            for proc in (1, 3):
                with self.subTest(test_file=test_file, cuff_file=cuff_file, proc=proc):
                    self.conf["prepare"]["files"]["gff"][0] = pkg_resources.resource_filename("Mikado.tests",
                                                                                              cuff_file)
                    self.conf["prepare"]["files"]["gff"][1] = pkg_resources.resource_filename("Mikado.tests",
                                                                                              test_file)
                    self.conf["prepare"]["files"]["out_fasta"] = "mikado_prepared.fasta"
                    self.conf["prepare"]["files"]["out"] = "mikado_prepared.gtf"
                    args.strip_cds = True
                    args.json_conf = self.conf
                    args.seed = 10
                    args.json_conf["threads"] = proc
                    args.exclude_redundant = False
                    prepare.prepare(args, self.logger)

                    # Now that the program has run, let's check the output
                    self.assertTrue(os.path.exists(os.path.join(self.conf["prepare"]["files"]["output_dir"],
                                                                "mikado_prepared.fasta")))
                    self.assertGreater(os.stat(os.path.join(self.conf["prepare"]["files"]["output_dir"],
                                                            "mikado_prepared.fasta")).st_size, 0)

                    fa = pysam.FastaFile(os.path.join(self.conf["prepare"]["files"]["output_dir"],
                                                      "mikado_prepared.fasta"))
                    res = dict((name, length) for name, length in zip(fa.references, fa.lengths))
                    fa.close()
                    precal = self.trinity_res.copy()
                    precal.update(self.cuff_results)
                    precal.pop("tr_" + self._trinity_redundant[0])
                    self.assertEqual(res, precal)
                    os.remove(os.path.join(self.conf["prepare"]["files"]["output_dir"],
                                           "mikado_prepared.fasta.fai"))
        dir.cleanup()

    @mark.slow
    def test_prepare_with_cds(self):

        rev_strand = {"+": "-", "-": "+"}

        self.conf["prepare"]["files"]["labels"] = ["ann"]
        folder = tempfile.TemporaryDirectory(prefix="test_prepare_with_cds")
        ann_gff3 = pkg_resources.resource_filename("Mikado.tests", "annotation.gff3")
        rev_ann_gff3 = tempfile.NamedTemporaryFile(suffix=".gff3", mode="wt", dir=folder.name)
        with open(ann_gff3) as ann:
            for line in ann:
                line = GffLine(line)
                if line.header is True:
                    continue
                line.strand = rev_strand[line.strand]  # Invert strand.
                print(line, file=rev_ann_gff3)
        rev_ann_gff3.flush()

        self.conf["prepare"]["files"]["gff"] = []
        self.conf["prepare"]["files"]["output_dir"] = folder.name
        self.conf["prepare"]["files"]["out_fasta"] = "mikado_prepared.fasta"
        self.conf["prepare"]["files"]["out"] = "mikado_prepared.gtf"
        args = Namespace()
        args.json_conf = self.conf

        for fname in [ann_gff3, rev_ann_gff3.name]:
            for strip in (True, False):
                for proc in (1, 3):
                    with self.subTest(fname=fname, strip=strip, proc=proc):
                        self.conf["prepare"]["files"]["gff"] = [fname]
                        args.json_conf["prepare"]["strip_cds"] = False
                        self.conf["prepare"]["files"]["strip_cds"] = [strip]
                        args.json_conf["threads"] = proc
                        with self.assertLogs(self.logger, "INFO") as cm:
                            try:
                                prepare.prepare(args, logger=self.logger)
                            except SystemExit:
                                raise SystemExit("\n".join(cm.output))
                        fasta = os.path.join(self.conf["prepare"]["files"]["output_dir"], "mikado_prepared.fasta")
                        self.assertTrue(os.path.exists(fasta), "\n".join(cm.output))
                        if strip is True or (strip is False and fname == ann_gff3):
                            self.assertGreater(os.stat(fasta).st_size, 0, "\n".join(cm.output))
                            fa = pyfaidx.Fasta(fasta)
                            self.assertEqual(len(fa.keys()), 2, "\n".join(cm.output))
                            fa.close()
                        else:
                            self.assertEqual(os.stat(fasta).st_size, 0,
                                             str(strip) + " " + str(fname) + "\n" + "\n".join(cm.output))

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
        rev_ann_gff3.close()

    @mark.slow
    def test_cdna_redundant_cds_not(self):
        """This test will verify whether the new behaviour of not considering redundant two models with same
        exon structure but different CDS does function properly."""

        gtf = pkg_resources.resource_filename("Mikado.tests", "cds_test_1.gtf")
        self.conf["prepare"]["files"]["gff"] = [gtf]
        self.conf["prepare"]["files"]["labels"] = [""]
        self.conf["prepare"]["files"]["strip_cds"] = [False]
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
                self.conf["prepare"]["files"]["exclude_redundant"] = [b]
                folder = tempfile.TemporaryDirectory()
                args.output_dir = folder.name
                args.seed = 10
                args.list = None
                args.gffs = None
                args.strand_specific_assemblies = None
                args.labels = None
                args.json_conf = self.conf
                args.exclude_redundant = b
                args.out, args.out_fasta = None, None
                args.json_conf["prepare"]["files"]["log"] = "prepare.log"
                if isinstance(args.json_conf["reference"]["genome"], bytes):
                    args.json_conf["reference"]["genome"] = args.json_conf["reference"]["genome"].decode()
                args.log = open(os.path.join(args.output_dir, "prepare.log"), "wt")
                self.logger.setLevel("DEBUG")
                args, _ = prepare_setup(args)
                self.assertEqual(args.output_dir, folder.name)
                self.assertEqual(args.json_conf["prepare"]["files"]["output_dir"], folder.name)
                self.assertIn(os.path.dirname(args.json_conf["prepare"]["files"]["out_fasta"]),
                              (folder.name, ""), args.json_conf)
                self.assertIn(os.path.dirname(args.json_conf["prepare"]["files"]["out"]),
                              (folder.name, ""), args.json_conf)

                with self.assertRaises(SystemExit) as exi:
                    prepare_launcher(args)
                self.assertTrue(os.path.exists(folder.name))
                self.assertTrue(os.path.isdir(folder.name))
                self.assertEqual(exi.exception.code, 0)
                self.assertTrue(os.path.exists(os.path.join(folder.name,
                                                            "mikado_prepared.fasta")),
                                open(os.path.join(folder.name,
                                                  "prepare.log")).read())
                fa = pyfaidx.Fasta(os.path.join(folder.name,
                                                "mikado_prepared.fasta"))
                logged = [_ for _ in open(args.json_conf["prepare"]["files"]["log"])]
                self.assertFalse("AT5G01530.1" in fa.keys(), (b, sorted(list(fa.keys()))))
                self.assertTrue("AT5G01530.2" in fa.keys())
                if b is False:
                    self.assertEqual(len(fa.keys()), 4)
                    self.assertEqual(sorted(fa.keys()), sorted(["AT5G01530."+str(_) for _ in [0, 2, 3, 4]]))
                else:
                    self.assertEqual(len(fa.keys()), 3, (fa.keys(), logged))
                    self.assertIn("AT5G01530.0", fa.keys())
                    self.assertIn("AT5G01530.2", fa.keys())
                    self.assertNotIn("AT5G01530.3", fa.keys())
                    self.assertIn("AT5G01530.4", fa.keys())
                gtf_file = os.path.join(folder.name, "mikado_prepared.gtf")
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
                    # self.assertIn("AT5G01530.3", transcripts)
                    a5 = transcripts["AT5G01530.2"]
                    self.assertTrue(a5.is_coding)
                    self.assertIn("has_start_codon", a5.attributes)
                    self.assertIn("has_stop_codon", a5.attributes)
                    self.assertTrue(a5.has_start_codon)
                    self.assertTrue(a5.has_stop_codon)
                    self.assertTrue(a5.is_complete)

                self.assertGreater(coding_count, 0)

    @mark.slow
    def test_negative_cdna_redundant_cds_not(self):
        """This test will verify whether the new behaviour of not considering redundant two models with same
        exon structure but different CDS does function properly."""

        gtf = pkg_resources.resource_filename("Mikado.tests", "cds_test_2.gtf")
        self.conf["prepare"]["files"]["gff"] = [gtf]
        dir = tempfile.TemporaryDirectory(prefix="test_negative_cdna_redundant_cds_not")
        self.conf["prepare"]["files"]["output_dir"] = dir.name
        self.conf["prepare"]["files"]["labels"] = [""]
        self.conf["prepare"]["files"]["out_fasta"] = "mikado_prepared.fasta"
        self.conf["prepare"]["files"]["strip_cds"] = [False]
        self.conf["prepare"]["files"]["out"] = "mikado_prepared.gtf"
        self.conf["prepare"]["files"]["log"] = "prepare.log"
        self.conf["prepare"]["strip_cds"] = False
        self.conf["prepare"]["minimum_cdna_length"] = 150  # Necessary for testing A5

        args = Namespace()
        args.strip_cds = False
        args.json_conf = self.conf
        for b in (False, True):
            with self.subTest(b=b):
                self.conf["prepare"]["files"]["exclude_redundant"] = [b]
                folder = tempfile.TemporaryDirectory()
                args.json_conf = self.conf
                args.json_conf["seed"] = 10
                args.exclude_redundant = b
                args.output_dir = folder.name
                args.log = None
                args.gff = None
                args.list = None
                args.strand_specific_assemblies = None
                if isinstance(args.json_conf["reference"]["genome"], bytes):
                    args.json_conf["reference"]["genome"] = args.json_conf["reference"]["genome"].decode()
                args, _ = prepare_setup(args)
                prepare.prepare(args, self.logger)
                self.assertTrue(os.path.exists(os.path.join(self.conf["prepare"]["files"]["output_dir"],
                                                            "mikado_prepared.fasta")))
                fa = pyfaidx.Fasta(os.path.join(self.conf["prepare"]["files"]["output_dir"],
                                                "mikado_prepared.fasta"))
                if b is False:
                    self.assertEqual(len(fa.keys()), 5)
                    self.assertEqual(sorted(fa.keys()), sorted(["AT5G01015." + str(_) for _ in
                                                                [0, 1, 3, 4, 5]]))
                else:
                    self.assertEqual(len(fa.keys()), 4, "\n".join(list(fa.keys())))
                    self.assertIn("AT5G01015.0", fa.keys())
                    self.assertTrue("AT5G01015.1" in fa.keys())
                    self.assertNotIn("AT5G01015.3", fa.keys())
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

                    a_first = transcripts["AT5G01015.1"]
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
                folder.cleanup()

    @mark.slow
    def test_truncated_cds(self):
        files = ["test_truncated_cds.gff3"]
        files = [pkg_resources.resource_filename("Mikado.tests", filename) for filename in files]
        self.conf["prepare"]["files"]["gff"] = files
        self.conf["prepare"]["files"]["strip_cds"] = [False] * len(files)
        self.conf["prepare"]["files"]["labels"] = [""]
        dir = tempfile.TemporaryDirectory(prefix="test_truncated_cds")
        self.conf["prepare"]["files"]["output_dir"] = dir.name
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

    @mark.slow
    def test_source_selection(self):

        # Chr5	TAIR10	mRNA	208937	210445	.	+	.	gene_id "AT5G01530"; transcript_id "AT5G01530.0";
        # Chr5	TAIR10	exon	208937	209593	.	+	.	gene_id "AT5G01530"; transcript_id "AT5G01530.0";
        # Chr5	TAIR10	exon	209881	210445	.	+	.	gene_id "AT5G01530"; transcript_id "AT5G01530.0";

        dir = tempfile.TemporaryDirectory(prefix="test_source_selection")
        self.conf["prepare"]["files"]["output_dir"] = dir.name
        self.conf["prepare"]["files"]["out_fasta"] = "mikado_prepared.fasta"
        self.conf["prepare"]["files"]["out"] = "mikado_prepared.gtf"
        self.conf["prepare"]["strip_cds"] = False
        self.conf["prepare"]["exclude_redundant"] = True
        self.conf["threads"] = 1

        self.conf["reference"]["genome"] = self.fai.filename.decode()

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
        self.conf["prepare"]["files"]["strip_cds"] = [True, True]

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

    @mark.slow
    def test_reference_selection(self):

        dir = tempfile.TemporaryDirectory(prefix="test_reference_selection")
        self.conf["prepare"]["files"]["output_dir"] = outdir = dir.name
        self.conf["prepare"]["files"]["out_fasta"] = "mikado_prepared.fasta"
        self.conf["prepare"]["files"]["out"] = "mikado_prepared.gtf"
        self.conf["prepare"]["files"]["strip_cds"] = [False, False]
        self.conf["prepare"]["files"]["exclude_redundant"] = [False, False]
        self.conf["prepare"]["strip_cds"] = False
        self.conf["prepare"]["exclude_redundant"] = True

        self.conf["reference"]["genome"] = self.fai.filename.decode()

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
                for iteration in range(2):  # Repeat each test 2 times, not more for time length reasons
                    with self.subTest(round=round, format=format, iteration=iteration,
                                      msg="Starting round {} ({})".format(round, rounds[round])):
                        t1_score, t2_score, is_ref, res, corr_strand = rounds[round]
                        self.conf["threads"] = 1
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

    @mark.slow
    def test_reference_cds_kept(self):

        t = Transcript()
        # This is *key*. Transcript T1 should never be selected, unless we are having a "lenient" analysis.
        # However, by flagging T1 as a "reference" transcript, we should retain it
        self.conf["prepare"]["lenient"] = False
        t.chrom, t.start, t.end, t.strand = "Chr5", 208937, 210445, "+"
        t.add_exons([(208937, 209593), (209881, 210445)])
        t.add_exons([(209084, 209593), (209881, 210243)], features=["CDS"] * 2)
        t.id = "file1.1"
        t.parent = "file1"
        t.finalize()

        # Likewise. We will also test that even when different scores are applied, we *still* retain all
        # transcripts.
        t2 = t.deepcopy()
        t2.id = "file2.1"
        t2.parent = "file2"

        dir = tempfile.TemporaryDirectory(prefix="test_reference_cds_kept")
        self.conf["prepare"]["files"]["output_dir"] = outdir = dir.name
        self.conf["prepare"]["files"]["out_fasta"] = "mikado_prepared.fasta"
        self.conf["prepare"]["files"]["out"] = "mikado_prepared.gtf"
        self.conf["prepare"]["strip_cds"] = True
        self.conf["prepare"]["exclude_redundant"] = True
        self.conf["reference"]["genome"] = self.fai.filename.decode()

        rounds = {
            # Standard cases. I expect the transcripts to be reversed, and one to be picked
            0: [0, 0, [], "rand", set()],
            1: [-2, -2, [], "rand", set()],
            2: [10, 10, [], "rand", set()],
            # T1 as reference. T1 should be kept as coding, T2 should be kept as non-coding.
            3: [0, 0, ["T1"], sorted(["T1_file1.1", "T2_file2.1"]), {"T1_file1.1"}],
            # Both as reference. Both should be kept, both should be on the negative strand.
            4: [1, 0, ["T1", "T2"], sorted(["T1_file1.1", "T2_file2.1"]), {"T1_file1.1", "T2_file2.1"}],
            5: [0, 0, ["T2"], sorted(["T1_file1.1", "T2_file2.1"]), {"T2_file2.1"}]

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
            self.conf["prepare"]["files"]["exclude_redundant"] = [True, True]
            self.conf["prepare"]["files"]["strip_cds"] = [True, True]

            for round in rounds:
                for iteration in range(2):  # Repeat each test 2 times for time length reasons
                    with self.subTest(round=round, format=format, iteration=iteration,
                                      msg="Starting round {} ({})".format(round, rounds[round])):
                        t1_score, t2_score, is_ref, res, coding = rounds[round]
                        self.conf["prepare"]["files"]["source_score"] = {"T1": t1_score,
                                                                         "T2": t2_score}
                        self.conf["prepare"]["files"]["reference"] = []
                        self.conf["threads"] = 1
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
                            self.assertEqual(key, res, (round, rounds[round], cm.output))
                        else:
                            self.assertEqual(len(fa.keys()), 1, (round, fa.keys(), res))
                        gtf = os.path.join(outdir, os.path.basename(
                            self.conf["prepare"]["files"]["out"]))
                        with_cds = set()
                        for line in to_gff(gtf):
                            if line.feature == "CDS":
                                with_cds.add(line.transcript)
                            else:
                                continue
                        self.assertEqual(coding, with_cds)


@mark.slow
class CompareCheck(unittest.TestCase):

    """Test to check that compare interacts correctly with match, match_part, cDNA_match"""

    def tearDown(self):
        logging.shutdown()

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
        dir = tempfile.mkdtemp(prefix="test_index")
        namespace.log = None
        logger = create_null_logger("null")

        for ref in files:
            with self.subTest(ref=ref):
                temp_ref = os.path.join(dir, ref)
                with pkg_resources.resource_stream("Mikado.tests", ref) as ref_handle,\
                        open(temp_ref, "wb") as out_handle:
                    out_handle.write(ref_handle.read())
                namespace.reference = to_gff(temp_ref)
                with self.assertLogs("main_compare") as ctx:
                    compare(namespace)

                self.assertTrue(os.path.exists("{}.midx".format(namespace.reference.name)))
                self.assertGreater(os.stat("{}.midx".format(namespace.reference.name)).st_size, 0)
                genes, positions = load_index(namespace, logger)
                self.assertIsInstance(genes, dict)
                self.assertIsInstance(positions, dict)
                self.assertEqual(len(genes), 38)
                os.remove(namespace.reference.name)
                os.remove("{}.midx".format(namespace.reference.name))
                namespace.reference.close()
        shutil.rmtree(dir)

    @mark.slow
    def test_compare_trinity(self):

        # Create the list of files
        files = ["trinity.gtf",
                 "trinity.gff3",
                 "trinity.cDNA_match.gff3",
                 "trinity.match_matchpart.gff3",
                 "trinity.bed12"]
        files = [pkg_resources.resource_filename("Mikado.tests", filename) for filename in files]
        bam = pkg_resources.resource_filename("Mikado.tests", "trinity.minimap2.bam")

        namespace = Namespace(default=False)
        namespace.distance = 2000
        namespace.no_save_index = True
        namespace.gzip = False

        for ref, pred in itertools.chain(itertools.permutations(files, 2),
                                         [(ref, bam) for ref in files]):
            with self.subTest(ref=ref, pred=pred):
                namespace.reference = to_gff(ref)
                namespace.prediction = to_gff(pred)
                namespace.processes = 2
                folder = tempfile.mkdtemp(prefix="test_compare_trinity_{}_{}".format(
                    os.path.splitext(ref)[-1], os.path.splitext(pred)[-1]))
                namespace.log = None
                if pred != bam:
                    namespace.out = os.path.join(folder,
                                                 "compare_{}_{}".format(
                        files.index(ref), files.index(pred)))
                else:
                    namespace.out = os.path.join(folder,
                                                 "compare_{}_{}".format(
                        files.index(ref), len(files) + 1))

                with self.assertLogs("main_compare") as ctx:
                    try:
                        compare(namespace)
                    except (ValueError, TypeError) as exc:
                        self.assertTrue(False, (ref, pred, exc))
                sleep(0.1)
                refmap = "{}.refmap".format(namespace.out)
                tmap = "{}.tmap".format(namespace.out)
                stats = "{}.stats".format(namespace.out)

                self.assertTrue(len(ctx.records) > 0)
                log = [str(_) for _ in ctx.records]
                for fname in [refmap, stats, tmap]:
                    self.assertTrue(os.path.exists(fname),
                                    (fname, ref, pred, glob.glob(namespace.out + "*"), "\n".join(log)))
                    self.assertGreater(os.stat(fname).st_size, 0, (fname, ref, pred, "\n".join(log)))

                with open(refmap) as _:
                    reader = csv.DictReader(_, delimiter="\t")
                    for counter, line in enumerate(reader, start=1):
                        ccode = line["ccode"]
                        if pred != bam:
                            self.assertIn(ccode,
                                          ("_", "=", "f,_", "f,="),
                                          (ref, pred, line))
                    self.assertEqual(counter, 38)
                if pred == bam:
                    with open(tmap) as _:
                        reader = csv.DictReader(_, delimiter="\t")
                        for counter, line in enumerate(reader, start=1):
                            pass
                    self.assertEqual(counter, 38)
                shutil.rmtree(folder)

    def test_compare_problematic(self):

        problematic = pkg_resources.resource_filename("Mikado.tests", "Chrysemys_picta_bellii_problematic.gff3")
        namespace = Namespace(default=False)
        namespace.distance = 2000
        namespace.no_save_index = True
        namespace.protein_coding = False
        namespace.exclude_utr = False
        namespace.self = False
        namespace.gzip = False

        folders = []
        for proc in (1, 3):
            namespace.reference = to_gff(problematic)
            namespace.prediction = to_gff(problematic)
            namespace.processes = proc
            folder = tempfile.mkdtemp(prefix="test_compare_problematic_{}".format(proc))
            namespace.log = os.path.join(folder, "compare_problematic_{proc}.log".format(proc=proc))
            namespace.out = os.path.join(folder, "compare_problematic_{proc}".format(proc=proc))
            folders.append(folder)
            compare(namespace)
            sleep(1)
            refmap = "{}.refmap".format(namespace.out)
            tmap = "{}.tmap".format(namespace.out)
            stats = "{}.stats".format(namespace.out)

            self.assertTrue(os.path.exists(namespace.log))
            with open(namespace.log) as log_handle:
                log = [_.rstrip() for _ in log_handle]
            for fname in [refmap, stats, tmap]:
                self.assertTrue(os.path.exists(fname), (glob.glob(namespace.out + "*"), "\n".join(log)))
                self.assertGreater(os.stat(fname).st_size, 0,
                                   "\n".join(log))
            with open(refmap) as _:
                reader = csv.DictReader(_, delimiter="\t")
                for counter, line in enumerate(reader, start=1):
                    ccode = line["ccode"]
                    self.assertIn(ccode, ("_", "=", "f,_", "f,="), line)
                self.assertEqual(counter, 4)
            with open(tmap) as _:
                reader = csv.DictReader(_, delimiter="\t")
                for counter, line in enumerate(reader, start=1):
                    pass
            self.assertEqual(counter, 4)
        [shutil.rmtree(folder) for folder in folders]


class ConfigureCheck(unittest.TestCase):

    """Test for creating configuration files"""

    __genomefile__ = None

    @classmethod
    def setUpClass(cls):
        cls.fai = pysam.FastaFile(pkg_resources.resource_filename("Mikado.tests", "chr5.fas.gz"))

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
        namespace.seed = None
        namespace.min_clustering_cds_overlap = 0.2
        namespace.min_clustering_cdna_overlap = 0.2
        dir = tempfile.TemporaryDirectory()
        out = os.path.join(dir.name, "configuration.yaml")
        with open(out, "w") as out_handle:
            namespace.out = out_handle
            sub_configure.create_config(namespace)
        self.assertGreater(os.stat(out).st_size, 0)
        conf = configuration.configurator.to_json(out)
        conf = configuration.configurator.check_json(conf)
        conf = configuration.configurator.check_json(conf)
        self.assertNotIn("asm_methods", conf)
        dir.cleanup()

    def test_seed(self):
        namespace = Namespace(default=False)
        namespace.scoring = None
        namespace.intron_range = None
        namespace.reference = ""
        namespace.external = None
        namespace.threads = 1
        namespace.blast_targets = []
        namespace.junctions = []
        namespace.new_scoring = None
        namespace.min_clustering_cds_overlap = 0.2
        namespace.min_clustering_cdna_overlap = 0.2
        dir = tempfile.TemporaryDirectory()
        out = os.path.join(dir.name, "configuration.yaml")
        for trial in (None, 1066, 175108):
            with self.subTest(trial=trial):
                namespace.mode = ["permissive"]
                namespace.seed = trial
                with open(out, "w") as out_handle:
                    namespace.out = out_handle
                    sub_configure.create_config(namespace)
                self.assertGreater(os.stat(out).st_size, 0)
                conf = configuration.configurator.to_json(out)
                conf = configuration.configurator.check_json(conf)
                conf = configuration.configurator.check_json(conf)
                self.assertNotIn("asm_methods", conf)
                if trial is not None:
                    self.assertEqual(conf["seed"], trial)
                else:
                    self.assertNotEqual(conf["seed"], trial)
                    self.assertIsInstance(conf["seed"], int)

        for mistake in (False, "hello", 10.5, b"890"):
            with self.subTest(mistake=mistake):
                namespace.mode = ["permissive"]
                with self.assertRaises(OSError):
                    namespace.seed = mistake
                    with open(out, "w") as out_handle:
                        namespace.out = out_handle
                        sub_configure.create_config(namespace)

        dir.cleanup()

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
        namespace.seed = None
        namespace.min_clustering_cds_overlap = 0.2
        namespace.min_clustering_cdna_overlap = 0.2
        dir = tempfile.TemporaryDirectory()
        out = os.path.join(dir.name, "configuration.yaml")
        with open(out, "w") as out_handle:
            namespace.out = out_handle
            sub_configure.create_config(namespace)
        self.assertGreater(os.stat(out).st_size, 0)
        conf = configuration.configurator.to_json(out)
        conf = configuration.configurator.check_json(conf)
        conf = configuration.configurator.check_json(conf)
        self.assertNotIn("asm_methods", conf)
        dir.cleanup()

    def test_mikado_config_daijin(self):
        namespace = Namespace(default=False)
        namespace.scoring = None
        namespace.intron_range = None
        namespace.reference = pkg_resources.resource_filename("Mikado.tests", "chr5.fas.gz")
        namespace.external = None
        namespace.mode = ["permissive"]
        namespace.threads = 1
        namespace.blast_targets = []
        namespace.junctions = []
        namespace.new_scoring = None
        namespace.full = True
        namespace.daijin = True
        namespace.seed = None
        namespace.min_clustering_cds_overlap = 0.2
        namespace.min_clustering_cdna_overlap = 0.2
        dir = tempfile.TemporaryDirectory()
        out = os.path.join(dir.name, "configuration.yaml")
        with open(out, "w") as out_handle:
            namespace.out = out_handle
            sub_configure.create_config(namespace)
        self.assertGreater(os.stat(out).st_size, 0)
        conf = configuration.configurator.to_json(out)
        conf = configuration.configurator.check_json(conf)
        conf = configuration.configurator.check_json(conf)
        dir.cleanup()

    def test_mikado_config_daijin_set_from_mode(self):
        namespace = Namespace(default=False)
        namespace.scoring = None
        namespace.intron_range = None
        namespace.reference = pkg_resources.resource_filename("Mikado.tests", "chr5.fas.gz")
        namespace.external = None
        namespace.mode = ["permissive", "split"]
        # namespace.mode = ["permissive"]
        namespace.threads = 1
        namespace.blast_targets = []
        namespace.junctions = []
        namespace.new_scoring = None
        namespace.full = True
        namespace.daijin = False
        namespace.seed = None
        namespace.min_clustering_cds_overlap = 0.2
        namespace.min_clustering_cdna_overlap = 0.2
        dir = tempfile.TemporaryDirectory()
        out = os.path.join(dir.name, "configuration.yaml")
        with open(out, "w") as out_handle:
            namespace.out = out_handle
            sub_configure.create_config(namespace)
        self.assertGreater(os.stat(out).st_size, 0)
        conf = configuration.configurator.to_json(out)
        conf = configuration.configurator.check_json(conf)
        conf = configuration.configurator.check_json(conf)
        dir.cleanup()

    @unittest.skipUnless((sys.version_info.minor > 4),
                         "Due to a bug in JSONSCHEMA, Daijin configure fails with Python versions lower than 3.5.")
    @mark.slow
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
        namespace.genome = self.fai.filename.decode()
        namespace.transcriptome = ""
        namespace.name = "Daijin"
        namespace.threads = 1
        namespace.full = False
        namespace.seed = None
        namespace.long_aln_methods = []

        for iteration in range(20):
            with self.subTest(iteration=iteration):
                dir = tempfile.TemporaryDirectory()
                namespace.out_dir = dir.name
                scorers = sorted(pkg_resources.resource_listdir("Mikado.configuration", "scoring_files"))

                namespace.scoring = scorers[np.random.choice(len(scorers))]

                out = os.path.join(dir.name, "configuration.yaml")
                with open(out, "wt") as out_handle:
                    namespace.out = out_handle
                    daijin_configurator.create_daijin_config(namespace, level="ERROR")
                self.assertGreater(os.stat(out).st_size, 0)

                with open(out) as out_handle:
                    config = yaml.load(out_handle, Loader=yLoader)

                daijin_configurator.check_config(config)
                dir.cleanup()


@mark.slow
class PickTest(unittest.TestCase):

    """This unit test will check that pick functions correctly."""

    def setUp(self):
        
        self.json_conf = configurator.to_json(None)
        self.json_conf["reference"]["genome"] = self.fai.filename.decode()

    @classmethod
    def setUpClass(cls):
        cls.fai = pysam.FastaFile(pkg_resources.resource_filename("Mikado.tests", "chr5.fas.gz"))

    def tearDown(self):
        all_child_threads = [thread for thread in threading.enumerate() if thread != threading.main_thread()]
        for thread in all_child_threads:
            lock = thread._tstate_lock
            if lock is not None:
                lock.acquire(False)
                lock.release()
            thread._wait_for_tstate_lock(block=True, timeout=0.00001)
            thread._stop()

    @mark.slow
    def test_single_proc(self):

        self.json_conf["threads"] = 1
        self.json_conf["db_settings"]["db"] = pkg_resources.resource_filename("Mikado.tests", "mikado.db")

        self.json_conf["pick"]["files"]["input"] = pkg_resources.resource_filename("Mikado.tests",
                                                                                   "mikado_prepared.gtf")
        dir = tempfile.TemporaryDirectory()
        self.json_conf["pick"]["files"]["output_dir"] = dir.name
        self.json_conf["pick"]["files"]["loci_out"] = "mikado.monoproc.loci.gff3"
        self.json_conf["pick"]["files"]["subloci_out"] = "mikado.monoproc.subloci.gff3"
        self.json_conf["pick"]["files"]["monoloci_out"] = "mikado.monoproc.monoloci.gff3"
        self.json_conf["pick"]["files"]["log"] = "mikado.monoproc.log"
        self.json_conf["pick"]["alternative_splicing"]["pad"] = False
        self.json_conf["log_settings"]["log_level"] = "WARNING"

        pick_caller = picker.Picker(json_conf=self.json_conf)
        with self.assertRaises(SystemExit), self.assertLogs("main_logger", "INFO") as cm:
            pick_caller()
        self.assertTrue(os.path.exists(os.path.join(dir.name, "mikado.monoproc.loci.gff3")))
        with to_gff(os.path.join(dir.name, "mikado.monoproc.loci.gff3")) as inp_gff:
            lines = [_ for _ in inp_gff if not _.header is True]
            self.assertGreater(len(lines), 0)
            self.assertGreater(len([_ for _ in lines if _.is_transcript is True]), 0)
            self.assertGreater(len([_ for _ in lines if _.feature == "mRNA"]), 0,
                               [_ for _ in cm.output if "WARNING" in _])
            self.assertGreater(len([_ for _ in lines if _.feature == "CDS"]), 0)

        dir.cleanup()

    @mark.slow
    def test_multi_proc(self):
        self.json_conf["threads"] = 2
        self.json_conf["pick"]["files"]["input"] = pkg_resources.resource_filename("Mikado.tests",
                                                                              "mikado_prepared.gtf")
        dir = tempfile.TemporaryDirectory()
        self.json_conf["pick"]["files"]["output_dir"] = dir.name
        self.json_conf["pick"]["files"]["loci_out"] = "mikado.multiproc.loci.gff3"
        self.json_conf["pick"]["files"]["subloci_out"] = "mikado.multiproc.subloci.gff3"
        self.json_conf["pick"]["files"]["monoloci_out"] = "mikado.multiproc.monoloci.gff3"
        self.json_conf["pick"]["files"]["log"] = "mikado.multiproc.log"
        self.json_conf["db_settings"]["db"] = pkg_resources.resource_filename("Mikado.tests", "mikado.db")
        self.json_conf["pick"]["alternative_splicing"]["pad"] = False
        self.json_conf["log_settings"]["log_level"] = "WARNING"

        pick_caller = picker.Picker(json_conf=self.json_conf)
        with self.assertRaises(SystemExit), self.assertLogs("main_logger", "INFO"):
            pick_caller()
        self.assertTrue(os.path.exists(os.path.join(dir.name, "mikado.multiproc.loci.gff3")))
        with to_gff(os.path.join(dir.name, "mikado.multiproc.loci.gff3")) as inp_gff:
            lines = [_ for _ in inp_gff if not _.header is True]
            self.assertGreater(len(lines), 0)
            self.assertGreater(len([_ for _ in lines if _.is_transcript is True]), 0)
            self.assertGreater(len([_ for _ in lines if _.feature == "mRNA"]), 0)
            self.assertGreater(len([_ for _ in lines if _.feature == "CDS"]), 0)

        dir.cleanup()

    @mark.slow
    def test_subprocess(self):
                
        self.json_conf["pick"]["files"]["input"] = pkg_resources.resource_filename("Mikado.tests",
                                                                                   "mikado_prepared.gtf")
        dir = tempfile.TemporaryDirectory()
        self.json_conf["pick"]["files"]["output_dir"] = dir.name
        self.json_conf["pick"]["files"]["loci_out"] = "mikado.subproc.loci.gff3"
        self.json_conf["pick"]["files"]["subloci_out"] = "mikado.subproc.subloci.gff3"
        self.json_conf["pick"]["files"]["monoloci_out"] = "mikado.subproc.monoloci.gff3"
        self.json_conf["pick"]["alternative_splicing"]["pad"] = False
        self.json_conf["pick"]["files"]["log"] = "mikado.subproc.log"
        self.json_conf["db_settings"]["db"] = str(pkg_resources.resource_filename("Mikado.tests", "mikado.db"))
        self.json_conf["log_settings"]["log_level"] = "WARNING"

        for num in (1, 2):
            with self.subTest(num=num):
                self.json_conf["pick"]["threads"] = num
                self.json_conf["pick"]["run_options"]["single_thread"] = (num == 1)
                json_file = os.path.join(dir.name, "mikado.yaml")

                # Printing out would crash without removing these compiled bits
                self.json_conf["requirements"].pop("compiled", None)
                self.json_conf["as_requirements"].pop("compiled", None)
                self.json_conf["not_fragmentary"].pop("compiled", None)

                with open(json_file, "wt") as json_handle:
                    print_config(yaml.dump(self.json_conf, default_flow_style=False), json_handle)

                sys.argv = ["mikado", "pick", "--json-conf", json_file, "--seed", "1078"]
                with self.assertRaises(SystemExit):
                    pkg_resources.load_entry_point("Mikado", "console_scripts", "mikado")()

                self.assertTrue(os.path.exists(os.path.join(dir.name, "mikado.subproc.loci.gff3")))
                with to_gff(os.path.join(dir.name, "mikado.subproc.loci.gff3")) as inp_gff:
                    lines = [_ for _ in inp_gff if not _.header is True]
                    self.assertGreater(len(lines), 0)
                    self.assertGreater(len([_ for _ in lines if _.is_transcript is True]), 0)
                    self.assertGreater(len([_ for _ in lines if _.feature == "mRNA"]), 0)
                    self.assertGreater(len([_ for _ in lines if _.feature == "CDS"]), 0)
                [os.remove(_) for _ in glob.glob(os.path.join(dir.name, "mikado.subproc.") + "*")]

        dir.cleanup()

    @mark.slow
    @unittest.skipUnless(os.path.exists("/dev/shm") and os.access("/dev/shm", os.W_OK),
                         "/dev/shm not present or not writeable")
    def test_subprocess_shm(self):
        self.json_conf["pick"]["files"]["input"] = pkg_resources.resource_filename("Mikado.tests",
                                                                                   "mikado_prepared.gtf")
        dir = tempfile.TemporaryDirectory()
        self.json_conf["pick"]["files"]["output_dir"] = dir.name
        self.json_conf["pick"]["files"]["loci_out"] = "mikado.subproc.loci.gff3"
        self.json_conf["pick"]["files"]["subloci_out"] = "mikado.subproc.subloci.gff3"
        self.json_conf["pick"]["files"]["monoloci_out"] = "mikado.subproc.monoloci.gff3"
        self.json_conf["pick"]["alternative_splicing"]["pad"] = False
        self.json_conf["pick"]["files"]["log"] = "mikado.subproc.log"
        self.json_conf["db_settings"]["db"] = str(pkg_resources.resource_filename("Mikado.tests", "mikado.db"))
        self.json_conf["log_settings"]["log_level"] = "WARNING"

        for num, shm in itertools.product((1, 2), (True,)):
            with self.subTest(num=num, shm=shm):

                self.json_conf["pick"]["threads"] = num
                self.json_conf["pick"]["run_options"]["single_thread"] = (num == 1)
                json_file = os.path.join(dir.name, "mikado.yaml")

                # Printing out would crash without removing these compiled bits
                self.json_conf["requirements"].pop("compiled", None)
                self.json_conf["as_requirements"].pop("compiled", None)
                self.json_conf["not_fragmentary"].pop("compiled", None)

                with open(json_file, "wt") as json_handle:
                    print_config(yaml.dump(self.json_conf, default_flow_style=False), json_handle)

                log = "pick.log"
                if os.path.exists(os.path.join(dir.name, log)):
                    os.remove(os.path.join(dir.name, log))
                sys.argv = ["mikado", "pick", "--json-conf", json_file, "--seed", "1078", "--log", log]
                if shm is True:
                    sys.argv.append("--shm")
                with self.assertRaises(SystemExit):
                    pkg_resources.load_entry_point("Mikado", "console_scripts", "mikado")()

                self.assertTrue(os.path.exists(os.path.join(dir.name, "mikado.subproc.loci.gff3")))
                with to_gff(os.path.join(dir.name, "mikado.subproc.loci.gff3")) as inp_gff:
                    lines = [_ for _ in inp_gff if not _.header is True]
                    self.assertGreater(len(lines), 0)
                    self.assertGreater(len([_ for _ in lines if _.is_transcript is True]), 0)
                    self.assertGreater(len([_ for _ in lines if _.feature == "mRNA"]), 0)
                    self.assertGreater(len([_ for _ in lines if _.feature == "CDS"]), 0)
                with open(os.path.join(dir.name, log)) as hlog:
                    log_lines = [_.rstrip() for _ in hlog]
                if shm is True:
                    self.assertTrue(any("Copying Mikado database into a SHM db" in _ for _ in log_lines))

                [os.remove(_) for _ in glob.glob(os.path.join(dir.name, "mikado.subproc.") + "*")]

        dir.cleanup()


    @mark.slow
    def test_different_scoring(self):

        dir = tempfile.TemporaryDirectory()
        self.json_conf["pick"]["files"]["output_dir"] = os.path.abspath(dir.name)
        self.json_conf["pick"]["files"]["input"] = pkg_resources.resource_filename("Mikado.tests",
                                                                                   "mikado_prepared.gtf")

        self.json_conf["pick"]["files"]["loci_out"] = "mikado.test_diff.loci.gff3"
        self.json_conf["pick"]["files"]["subloci_out"] = "mikado.test_diff.subloci.gff3"
        self.json_conf["pick"]["files"]["monoloci_out"] = "mikado.test_diff.monoloci.gff3"
        self.json_conf["pick"]["files"]["log"] = "mikado.test_diff.log"
        self.json_conf["pick"]["alternative_splicing"]["pad"] = False
        self.json_conf["log_settings"]["log_level"] = "DEBUG"

        self.assertEqual(os.path.basename(self.json_conf["pick"]["scoring_file"]),
                         "plant.yaml")
        shutil.copy(pkg_resources.resource_filename("Mikado.tests", "mikado.db"),
                    os.path.join(self.json_conf["pick"]["files"]["output_dir"], "mikado.db"))
        self.json_conf["db_settings"]["db"] = os.path.join(self.json_conf["pick"]["files"]["output_dir"],
                                                           "mikado.db")
        json_file = os.path.join(self.json_conf["pick"]["files"]["output_dir"], "mikado.yaml")
        with open(json_file, "wt") as json_handle:
            print_config(yaml.dump(self.json_conf, default_flow_style=False), json_handle)
        sys.argv = ["mikado", "pick", "--json-conf", json_file, "--single", "--seed", "1078"]
        with self.assertRaises(SystemExit):
            pkg_resources.load_entry_point("Mikado", "console_scripts", "mikado")()

        import csv
        with open(os.path.join(self.json_conf["pick"]["files"]["output_dir"], "mikado.test_diff.loci.scores.tsv")) as tsv:
            reader = csv.DictReader(tsv, delimiter="\t")
            score_names = [_ for _ in self.json_conf["scoring"]]
            score_header = [_ for _ in reader.fieldnames if _ not in
                            ("tid", "alias", "parent", "score", "source_score")]
            self.assertEqual(sorted(score_names), sorted(score_header))
        dir.cleanup()

    @mark.slow
    def test_different_scoring_2(self):

        self.json_conf["pick"]["files"]["input"] = pkg_resources.resource_filename("Mikado.tests",
                                                                                   "mikado_prepared.gtf")

        self.json_conf["pick"]["files"]["loci_out"] = "mikado.test_diff.loci.gff3"
        self.json_conf["pick"]["files"]["subloci_out"] = "mikado.test_diff.subloci.gff3"
        self.json_conf["pick"]["files"]["monoloci_out"] = "mikado.test_diff.monoloci.gff3"
        self.json_conf["pick"]["files"]["log"] = "mikado.test_diff.log"
        self.json_conf["pick"]["alternative_splicing"]["pad"] = False
        self.json_conf["log_settings"]["log_level"] = "DEBUG"

        self.assertEqual(os.path.basename(self.json_conf["pick"]["scoring_file"]),
                         "plant.yaml")

        outdir = tempfile.TemporaryDirectory()
        shutil.copy(pkg_resources.resource_filename("Mikado.tests", "mikado.db"),
                    os.path.join(outdir.name, "mikado.db"))
        self.json_conf["db_settings"]["db"] = os.path.join(outdir.name, "mikado.db")
        self.json_conf["pick"]["files"]["output_dir"] = os.path.join(outdir.name)
        json_file = os.path.join(outdir.name, "mikado.yaml")
        with open(json_file, "wt") as json_handle:
            print_config(yaml.dump(self.json_conf, default_flow_style=False),
                                                      json_handle)
        self.json_conf["pick"]["files"]["output_dir"] = os.path.join(outdir.name)
        scoring_file = pkg_resources.resource_filename("Mikado.tests", "scoring_only_cds.yaml")
        sys.argv = ["mikado", "pick", "--json-conf", json_file, "--single",
                    "--scoring-file", scoring_file, "--seed", "1078"]

        with self.assertRaises(SystemExit):
            pkg_resources.load_entry_point("Mikado", "console_scripts", "mikado")()

        import csv
        with open(os.path.join(self.json_conf["pick"]["files"]["output_dir"],
                               "mikado.test_diff.loci.scores.tsv")) as tsv:
            reader = csv.DictReader(tsv, delimiter="\t")
            score_header = [_ for _ in reader.fieldnames if _ not in
                            ("tid", "alias", "parent", "score", "source_score")]
            self.assertEqual(score_header, ["selected_cds_length"])

        outdir.cleanup()

    def __get_purgeable_gff(self):

        gtf = """Chr1	foo	transcript	100	1000	.	+	.	gene_id "foo1"; transcript_id "foo1.1"
        Chr1	foo	exon	100	1000	.	+	.	gene_id "foo1"; transcript_id "foo1.1"
        Chr1	foo	transcript	100	2000	.	+	.	gene_id "foo1"; transcript_id "foo1.2"
        Chr1	foo	exon	100	800	.	+	.	gene_id "foo1"; transcript_id "foo1.2"
        Chr1	foo	exon	1900	2000	.	+	.	gene_id "foo1"; transcript_id "foo1.2"
        Chr1	foo	transcript	10000	20000	.	+	.	gene_id "foo2"; transcript_id "foo2.1"
        Chr1	foo	exon	10000	13000	.	+	.	gene_id "foo2"; transcript_id "foo2.1"
        Chr1	foo	exon	19000	20000	.	+	.	gene_id "foo"; transcript_id "foo2.1\""""

        dir = tempfile.TemporaryDirectory()
        temp_gtf = tempfile.NamedTemporaryFile(mode="wt", suffix=".gtf", dir=dir.name, delete=True)

        temp_gtf.write(gtf)
        temp_gtf.flush()

        self.json_conf["pick"]["files"]["input"] = temp_gtf.name
        self.json_conf["db_settings"]["db"] = os.path.join(dir.name, "mikado.db")
        self.json_conf["pick"]["files"]["output_dir"] = dir.name
        self.json_conf["log_settings"]["log_level"] = "WARNING"
        self.json_conf["pick"]["alternative_splicing"]["pad"] = False  # Necessary!
        del self.json_conf["scoring"]
        del self.json_conf["requirements"]
        del self.json_conf["as_requirements"]
        del self.json_conf["not_fragmentary"]
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

        return gtf, dir, temp_gtf, scoring

    @mark.slow
    def test_purging1(self):

        # Now the scoring
        gtf, dir, temp_gtf, scoring = self.__get_purgeable_gff()

        scoring["scoring"] = dict()
        scoring["scoring"]["cdna_length"] = dict()
        scoring["scoring"]["cdna_length"]["rescaling"] = "max"
        scoring["scoring"]["cdna_length"]["filter"] = dict()
        scoring["scoring"]["cdna_length"]["filter"]["operator"] = "gt"
        scoring["scoring"]["cdna_length"]["filter"]["value"] = 2000

        scoring_file = tempfile.NamedTemporaryFile(suffix=".yaml", delete=True, mode="wt", dir=dir.name)
        yaml.dump(scoring, scoring_file)
        scoring_file.flush()
        self.json_conf["pick"]["scoring_file"] = scoring_file.name

        for purging in (False, True):
            with self.subTest(purging=purging):
                self.json_conf["pick"]["files"]["loci_out"] = "mikado.purging_{}.loci.gff3".format(purging)
                self.json_conf["pick"]["files"]["log"] = os.path.join(
                    dir.name,
                    "mikado.purging_{}.log".format(purging))
                self.json_conf["pick"]["clustering"]["purge"] = purging
                self.json_conf["pick"]["scoring_file"] = scoring_file.name
                self.json_conf = configurator.check_json(self.json_conf)
                self.assertEqual(len(self.json_conf["scoring"].keys()), 1, self.json_conf["scoring"].keys())

                pick_caller = picker.Picker(json_conf=self.json_conf)
                with self.assertRaises(SystemExit), self.assertLogs("main_logger", "INFO"):
                    pick_caller()

                with to_gff(os.path.join(dir.name,
                                         self.json_conf["pick"]["files"]["loci_out"])) as gff:
                    lines = [line for line in gff if line.header is False]
                self.assertGreater(len(lines), 0)
                self.assertTrue(any([_ for _ in lines if _.attributes.get("alias", "") == "foo2.1"]),
                                "\n".join([str(_) for _ in lines]))
                if purging is True:
                    self.assertFalse(any([_ for _ in lines if _.attributes.get("alias", "") in ("foo1.2", "foo1.1")]))
                else:
                    found_line = [_ for _ in lines if _.attributes.get("alias", "") in ("foo1.2", "foo1.1")]
                    self.assertTrue(any(found_line))
                    self.assertTrue(any([_ for _ in found_line if _.score == 0]))

            # Clean up
            for fname in ["mikado.db", "mikado.purging_{}.*".format(purging)]:
                [os.remove(_) for _ in glob.glob(os.path.join(dir.name, fname))]

        scoring_file.close()
        temp_gtf.close()
        dir.cleanup()

    @mark.slow
    def test_purging2(self):

        gtf, folder, temp_gtf, scoring = self.__get_purgeable_gff()

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

        scoring_file = tempfile.NamedTemporaryFile(suffix=".yaml", delete=True, mode="wt", dir=folder.name)
        yaml.dump(scoring, scoring_file)
        scoring_file.flush()
        self.json_conf["pick"]["scoring_file"] = scoring_file.name

        for purging in (False, True):
            with self.subTest(purging=purging):
                self.json_conf["pick"]["files"]["loci_out"] = "mikado.purging_{}.loci.gff3".format(purging)
                self.json_conf["pick"]["files"]["subloci_out"] = "mikado.purging_{}.subloci.gff3".format(purging)
                self.json_conf["pick"]["files"]["log"] = os.path.join(
                    folder.name,
                    "mikado.purging_{}.log".format(purging))
                self.json_conf["pick"]["clustering"]["purge"] = purging
                self.json_conf["pick"]["scoring_file"] = scoring_file.name
                self.json_conf = configurator.check_json(self.json_conf)
                self.assertEqual(len(self.json_conf["scoring"].keys()), 2,
                                 self.json_conf["scoring"].keys())

                pick_caller = picker.Picker(json_conf=self.json_conf)
                with self.assertRaises(SystemExit), self.assertLogs("main_logger", "INFO"):
                    pick_caller()

                with to_gff(os.path.join(folder.name,
                                         self.json_conf["pick"]["files"]["loci_out"])) as gff:
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
                [os.remove(_) for _ in glob.glob(os.path.join(folder.name, fname))]

        temp_gtf.close()
        scoring_file.close()
        folder.cleanup()

    @mark.slow
    def test_purging3(self):

        gtf, dir, temp_gtf, scoring = self.__get_purgeable_gff()
        temp_gtf.close()  # We are going to redo this

        scoring["scoring"] = dict()
        scoring["scoring"]["cdna_length"] = dict()
        scoring["scoring"]["cdna_length"]["rescaling"] = "min"
        scoring["scoring"]["cdna_length"]["multiplier"] = -10
        scoring["scoring"]["cdna_length"]["filter"] = dict()
        scoring["scoring"]["cdna_length"]["filter"]["operator"] = "lt"
        scoring["scoring"]["cdna_length"]["filter"]["value"] = 1000

        scoring["scoring"]["exon_num"] = dict()
        scoring["scoring"]["exon_num"]["rescaling"] = "max"

        scoring_file = tempfile.NamedTemporaryFile(suffix=".yaml", delete=True, mode="wt", dir=dir.name)
        yaml.dump(scoring, scoring_file)
        scoring_file.flush()
        self.json_conf["pick"]["scoring_file"] = scoring_file.name

        temp_gtf = tempfile.NamedTemporaryFile(mode="wt", suffix=".gtf", delete=True, dir=dir.name)

        gtf = "\n".join([_ for _ in gtf.split("\n") if "foo1.1" not in _])

        temp_gtf.write(gtf)
        temp_gtf.flush()
        self.json_conf["pick"]["files"]["input"] = temp_gtf.name

        for purging in (False, True):
            with self.subTest(purging=purging):
                self.json_conf["pick"]["files"]["loci_out"] = "mikado.purging_{}.loci.gff3".format(purging)
                self.json_conf["pick"]["files"]["subloci_out"] = "mikado.purging_{}.subloci.gff3".format(purging)
                self.json_conf["pick"]["files"]["log"] = "mikado.purging_{}.log".format(purging)
                self.json_conf["pick"]["clustering"]["purge"] = purging
                self.json_conf["pick"]["scoring_file"] = scoring_file.name
                self.json_conf = configurator.check_json(self.json_conf)
                self.assertEqual(len(self.json_conf["scoring"].keys()), 2, self.json_conf["scoring"].keys())

                pick_caller = picker.Picker(json_conf=self.json_conf)
                with self.assertRaises(SystemExit), self.assertLogs("main_logger", "INFO"):
                    pick_caller()

                with to_gff(os.path.join(dir.name,
                                         self.json_conf["pick"]["files"]["loci_out"])) as gff:
                    lines = [line for line in gff if line.header is False]
                self.assertGreater(len(lines), 0)
                self.assertTrue(any([_ for _ in lines if _.attributes.get("alias", "") == "foo2.1"]),
                                "\n".join([str(_) for _ in lines]))
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
        temp_gtf.close()
        scoring_file.close()
        dir.cleanup()


@mark.slow
class SerialiseChecker(unittest.TestCase):

    def setUp(self):
        self.json_conf = configurator.to_json(None)
        self.json_conf["reference"]["genome"] = self.fai.filename.decode()

    @classmethod
    def setUpClass(cls):
        cls.fai = pysam.FastaFile(pkg_resources.resource_filename("Mikado.tests", "chr5.fas.gz"))

    def test_subprocess_single(self):

        xml = pkg_resources.resource_filename("Mikado.tests", "chunk-001-proteins.xml.gz")
        transcripts = pkg_resources.resource_filename("Mikado.tests", "mikado_prepared.fasta")
        junctions = pkg_resources.resource_filename("Mikado.tests", "junctions.bed")
        orfs = pkg_resources.resource_filename("Mikado.tests", "transcripts.fasta.prodigal.gff3")
        uniprot = pkg_resources.resource_filename("Mikado.tests", "uniprot_sprot_plants.fasta.gz")
        mobjects = 300  # Let's test properly the serialisation for BLAST

        dir = tempfile.TemporaryDirectory()
        json_file = os.path.join(dir.name, "mikado.yaml")
        db = os.path.join(dir.name, "mikado.db")
        log = os.path.join(dir.name, "serialise.log")
        uni_out = os.path.join(dir.name, "uniprot_sprot_plants.fasta")
        with gzip.open(uniprot, "rb") as uni, open(uni_out, "wb") as uni_out_handle:
            uni_out_handle.write(uni.read())

        with open(json_file, "wt") as json_handle:
            print_config(yaml.dump(self.json_conf, default_flow_style=False),
                                                      json_handle)
        # Set up the command arguments
        for procs in (1,):
            with self.subTest(proc=procs):
                sys.argv = [str(_) for _ in ["mikado", "serialise", "--json-conf", json_file,
                            "--transcripts", transcripts, "--blast_targets", uni_out,
                            "--orfs", orfs, "--junctions", junctions, "--xml", xml, "-od", dir.name,
                            "-p", procs, "-mo", mobjects, "--log", os.path.basename(log), "--seed", "1078",
                                             os.path.basename(db)]]
                pkg_resources.load_entry_point("Mikado", "console_scripts", "mikado")()
                logged = [_.rstrip() for _ in open(log)]

                self.assertTrue(os.path.exists(db))
                conn = sqlite3.connect(db)
                cursor = conn.cursor()
                self.assertEqual(cursor.execute("select count(*) from hit").fetchall()[0][0], 562, logged)
                self.assertEqual(cursor.execute("select count(*) from hsp").fetchall()[0][0], 669)
                self.assertEqual(cursor.execute("select count(distinct(query_id)) from hsp").fetchall()[0][0], 71)
                self.assertEqual(cursor.execute("select count(distinct(query_id)) from hit").fetchall()[0][0], 71)
                self.assertEqual(cursor.execute("select count(distinct(target_id)) from hsp").fetchall()[0][0], 32)
                self.assertEqual(cursor.execute("select count(distinct(target_id)) from hit").fetchall()[0][0], 32)
                self.assertEqual(cursor.execute("select count(*) from junctions").fetchall()[0][0], 372)
                self.assertEqual(cursor.execute("select count(distinct(chrom_id)) from junctions").fetchall()[0][0], 2)
                self.assertEqual(cursor.execute("select count(*) from orf").fetchall()[0][0], 168,
                                 "\n".join(logged))
                self.assertEqual(cursor.execute("select count(distinct(query_id)) from orf").fetchall()[0][0], 81)
                os.remove(db)
        dir.cleanup()

    def test_subprocess_multi(self):

        xml = pkg_resources.resource_filename("Mikado.tests", "chunk-001-proteins.xml.gz")
        transcripts = pkg_resources.resource_filename("Mikado.tests", "mikado_prepared.fasta")
        junctions = pkg_resources.resource_filename("Mikado.tests", "junctions.bed")
        orfs = pkg_resources.resource_filename("Mikado.tests", "transcripts.fasta.prodigal.gff3")
        uniprot = pkg_resources.resource_filename("Mikado.tests", "uniprot_sprot_plants.fasta.gz")
        mobjects = 300  # Let's test properly the serialisation for BLAST

        # Set up the command arguments
        for procs in (1, 3,):
            with self.subTest(proc=procs):
                dir = tempfile.TemporaryDirectory(suffix="test_subprocess_multi_{}".format(procs))
                json_file = os.path.join(dir.name, "mikado.yaml")
                db = os.path.join(dir.name, "mikado.db")
                log = os.path.join(dir.name, "serialise.log")
                uni_out = os.path.join(dir.name, "uniprot_sprot_plants.fasta")
                with gzip.open(uniprot, "rb") as uni, open(uni_out, "wb") as uni_out_handle:
                    uni_out_handle.write(uni.read())

                with open(json_file, "wt") as json_handle:
                    print_config(yaml.dump(self.json_conf, default_flow_style=False),
                                               json_handle)
                sys.argv = [str(_) for _ in ["mikado", "serialise", "--json-conf", json_file,
                            "--transcripts", transcripts, "--blast_targets", uni_out,
                            "--orfs", orfs, "--junctions", junctions, "--xml", xml, "-od", dir.name,
                            "-p", procs, "-mo", mobjects, "--log", os.path.basename(log),
                                             "--seed", "1078", os.path.basename(db)]]
                pkg_resources.load_entry_point("Mikado", "console_scripts", "mikado")()
                logged = [_.rstrip() for _ in open(log)]

                self.assertTrue(os.path.exists(db))
                conn = sqlite3.connect(db)
                cursor = conn.cursor()
                self.assertEqual(cursor.execute("select count(*) from hit").fetchall()[0][0], 562, logged)
                self.assertEqual(cursor.execute("select count(*) from hsp").fetchall()[0][0], 669)
                self.assertEqual(cursor.execute("select count(distinct(query_id)) from hsp").fetchall()[0][0], 71)
                self.assertEqual(cursor.execute("select count(distinct(query_id)) from hit").fetchall()[0][0], 71)
                self.assertEqual(cursor.execute("select count(distinct(target_id)) from hsp").fetchall()[0][0], 32)
                self.assertEqual(cursor.execute("select count(distinct(target_id)) from hit").fetchall()[0][0], 32)
                self.assertEqual(cursor.execute("select count(*) from junctions").fetchall()[0][0], 372)
                self.assertEqual(cursor.execute("select count(distinct(chrom_id)) from junctions").fetchall()[0][0], 2)
                self.assertEqual(cursor.execute("select count(*) from orf").fetchall()[0][0], 168,
                                 "\n".join(logged))
                self.assertEqual(cursor.execute("select count(distinct(query_id)) from orf").fetchall()[0][0], 81)
                os.remove(db)
                dir.cleanup()

    @mark.slow
    def test_xml_vs_tsv(self):
        xml = pkg_resources.resource_filename("Mikado.tests", os.path.join("blast_data", "diamond.0.9.30.xml.gz"))
        tsv = pkg_resources.resource_filename("Mikado.tests", os.path.join("blast_data", "diamond.0.9.30.tsv.gz"))
        queries = pkg_resources.resource_filename("Mikado.tests", os.path.join("blast_data", "transcripts.fasta"))
        prots = pkg_resources.resource_filename("Mikado.tests", "uniprot_sprot_plants.fasta.gz")
        logs = dict()
        dbs = dict()
        base = tempfile.TemporaryDirectory()
        for name, blast in zip(["xml", "tsv"], [xml, tsv]):
            db = "{}.db".format(name)
            log = "{}.log".format(name)
            sys.argv = [str(_) for _ in ["mikado", "serialise", "-od", base.name,
                                         "--transcripts", queries, "--blast_targets", prots,
                                         "--xml", xml, "-mo", 1000, "--log", log, "--seed", "1078",
                                         db]]
            pkg_resources.load_entry_point("Mikado", "console_scripts", "mikado")()
            dbs[name] = os.path.join(base.name, db)
            logged = [_.rstrip() for _ in open(os.path.join(base.name, log))]
            logs[name] = logged

        def prep_dbs(name):
            hsp, hit, query, target = [pd.read_sql(table, name) for table in ["hsp", "hit", "query", "target"]]
            hit = hit.join(target.set_index("target_id"), on=["target_id"], how="inner").join(
                query.set_index("query_id"), on=["query_id"], how="inner")
            hsp = hsp.join(target.set_index("target_id"), on=["target_id"], how="inner").join(
                query.set_index("query_id"), on=["query_id"], how="inner")
            hsp.set_index(["query_name", "target_name", "counter"], inplace=True)
            hit.set_index(["query_name", "target_name"], inplace=True)
            return hit, hsp

        xml_hit, xml_hsp = prep_dbs("sqlite:///" + dbs["xml"])
        tsv_hit, tsv_hsp = prep_dbs("sqlite:///" + dbs["tsv"])
        hit = pd.merge(xml_hit, tsv_hit, left_index=True, right_index=True, suffixes=("_xml", "_tsv"))
        hsp = pd.merge(xml_hsp, tsv_hsp, left_index=True, right_index=True, suffixes=("_xml", "_tsv"))
        self.assertTrue(hit.shape[0] == xml_hit.shape[0] == tsv_hit.shape[0] > 0)
        self.assertTrue(hsp.shape[0] == xml_hsp.shape[0] == tsv_hsp.shape[0] > 0)
        # Get the columns
        hitcols, hspcols = dict(), dict()
        for d, df in zip([hitcols, hspcols], [hit, hsp]):
            for col in df.columns:
                name = col[:-4]
                if name not in d:
                    d[name] = []
                d[name].append(col)
            failed = []
            for col in d:
                if col in ("query_id", "target_id"):
                    continue
                catch = df[d[col]].apply(lambda row: row[0] == row[1] or
                                                     np.isclose(row[0], row[1], atol=.01, rtol=.01), axis=1)
                if not (catch).all():
                    failed.append(col)
            self.assertEqual(len(failed), 0, failed)

    def test_subprocess_multi_empty_orfs(self):

        xml = pkg_resources.resource_filename("Mikado.tests", "chunk-001-proteins.xml.gz")
        transcripts = pkg_resources.resource_filename("Mikado.tests", "mikado_prepared.fasta")
        junctions = pkg_resources.resource_filename("Mikado.tests", "junctions.bed")
        # orfs = pkg_resources.resource_filename("Mikado.tests", "transcripts.fasta.prodigal.gff3")
        tmp_orf = tempfile.NamedTemporaryFile(suffix=".bed12")
        tmp_orf.write(b"#track\n")
        tmp_orf.write(
            b"cufflinks_star_at.23553.1\t0\t1733\tID=1_1;partial=01;start_type=ATG\t0\t+\t312\t1733\t0,0,0\t1\t1733\t0\n")
        tmp_orf.flush()
        uniprot = pkg_resources.resource_filename("Mikado.tests", "uniprot_sprot_plants.fasta.gz")
        mobjects = 300  # Let's test properly the serialisation for BLAST

        # Set up the command arguments
        for procs in (3, 1):
            with self.subTest(procs=procs):
                dir = tempfile.TemporaryDirectory(prefix="has_to_fail")
                json_file = os.path.join(dir.name, "mikado.yaml")
                db = os.path.join(dir.name, "mikado.db")
                log = "failed_serialise.log"
                uni_out = os.path.join(dir.name, "uniprot_sprot_plants.fasta")
                self.json_conf["serialise"]["files"]["log"] = os.path.basename(log)
                self.json_conf["multiprocessing_method"] = "fork"
                with gzip.open(uniprot, "rb") as uni, open(uni_out, "wb") as uni_out_handle:
                    uni_out_handle.write(uni.read())

                with open(json_file, "wt") as json_handle:
                    print_config(yaml.dump(self.json_conf, default_flow_style=False), json_handle)
                with self.subTest(proc=procs):
                    sys.argv = [str(_) for _ in ["mikado", "serialise", "--json-conf", json_file,
                                                 "--transcripts", transcripts, "--blast_targets", uni_out,
                                                 "--log", log,
                                                 "-od", dir.name,
                                                 "--orfs", tmp_orf.name, "--junctions", junctions, "--xml", xml,
                                                 "-p", procs, "-mo", mobjects, db,
                                                 "--seed", "1078"]]
                    log = os.path.join(dir.name, log)
                    with self.assertRaises(SystemExit):
                        pkg_resources.load_entry_point("Mikado", "console_scripts", "mikado")()
                    self.assertTrue("failed" in log)
                    self.assertTrue(os.path.exists(log), log)
                    self.assertTrue(os.stat(log).st_size > 0, log)
                    logged = [_.rstrip() for _ in open(log)]
                    self.assertGreater(len(logged), 0)
                    self.assertFalse(os.path.exists(db), logged)
                    self.assertTrue(any(
                        "Mikado serialise failed due to problems with the input data. Please check the logs." in line
                        for line in logged))
                    self.assertTrue(any(
                        "The provided ORFs do not match the transcripts provided and already present in the database."
                        in line for line in logged),
                    print("\n".join(logged)))
                dir.cleanup()

    def test_serialise_external(self):

        base = pkg_resources.resource_filename("Mikado.tests", "test_external")
        external_conf = os.path.join(base, "mikado.configuration.testds.yaml")
        external_scores = os.path.join(base, "annotation_run1.metrics.testds.txt")
        fasta = os.path.join(base, "mikado_prepared.testds.fasta")

        for procs in (1, 3):
            with self.subTest(procs=procs):
                dir = tempfile.TemporaryDirectory(suffix="test_serialise_external")
                log = "serialise.log"
                sys.argv = [str(_) for _ in ["mikado", "serialise", "--json-conf", external_conf,
                                             "--transcripts", fasta, "-od", dir.name,
                                             "-l", log,
                                             "--external-scores", external_scores,
                                             "--seed", 10, "--procs", procs, "mikado.db"]]
                log = os.path.join(dir.name, log)
                pkg_resources.load_entry_point("Mikado", "console_scripts", "mikado")()
                conn = sqlite3.connect(os.path.join(dir.name, "mikado.db"))
                total = conn.execute("SELECT count(*) FROM external").fetchone()[0]
                self.assertEqual(total, 190)
                tot_sources = conn.execute("SELECT count(*) FROM external_sources").fetchone()[0]
                self.assertEqual(tot_sources, 95)


class StatsTest(unittest.TestCase):

    def test_annotation_stats(self):

        annotation_file = pkg_resources.resource_filename("Mikado.tests", "annotation.gff3")
        annotation_check = pkg_resources.resource_filename("Mikado.tests", "annotation.gff3.stats")

        dir = tempfile.TemporaryDirectory(prefix="test_annotation_stats")
        out = os.path.join(dir.name, "annotation.gff3.stats")
        sys.argv = [str(_) for _ in ["mikado", "util", "stats", annotation_file, out]]
        # with self.assertRaises(SystemExit):
        #     pkg_resources.load_entry_point("Mikado", "console_scripts", "mikado")()
        pkg_resources.load_entry_point("Mikado", "console_scripts", "mikado")()

        self.assertTrue(os.path.exists(out))

        with open(annotation_check) as check:
            check_lines = [_.strip() for _ in check]

        with open(out) as out_handle:
            out_lines = [_.strip() for _ in out_handle]

        self.assertEqual(check_lines, out_lines)

    def test_stat(self):

        """This unit test takes care of verifying that statistics are generated correctly when
            considering four different inputs. Output will be checked against a standard file."""

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
                dir = tempfile.TemporaryDirectory(prefix="test_stat")
                with open(os.path.join(dir.name,
                                       "{}.txt".format(os.path.basename(filename))), "w") as out:
                    namespace.out = out
                    Calculator(namespace)()
                self.assertGreater(os.stat(out.name).st_size, 0)
                with open(out.name) as out_handle:
                    lines = [_.rstrip() for _ in out_handle]
                self.assertEqual(std_lines, lines)
                os.remove(out.name)
                namespace.gff.close()
                dir.cleanup()

    def test_problematic(self):

        files = [pkg_resources.resource_filename("Mikado.tests", "Chrysemys_picta_bellii_problematic.gff3")]

        std_lines = []
        with pkg_resources.resource_stream("Mikado.tests", "Chrysemys_picta_bellii_problematic.txt") as t_stats:
            for line in t_stats:
                std_lines.append(line.decode().rstrip())

        namespace = Namespace(default=False)
        namespace.tab_stats = None
        for filename in files:
            with self.subTest(filename=filename):
                namespace.gff = to_gff(filename)
                dir = tempfile.TemporaryDirectory(prefix="test_problematic")
                with open(os.path.join(dir.name,
                                       "{}.txt".format(os.path.basename(filename))), "w") as out:
                    namespace.out = out
                    Calculator(namespace)()
                self.assertGreater(os.stat(out.name).st_size, 0)
                with open(out.name) as out_handle:
                    lines = [_.rstrip() for _ in out_handle]
                self.assertEqual(std_lines, lines)
                os.remove(out.name)
                namespace.gff.close()
                dir.cleanup()


class GrepTest(unittest.TestCase):

    @mark.slow
    def test_grep(self):
        files = [pkg_resources.resource_filename("Mikado.tests", fname) for fname in ["trinity.gtf", "trinity.gff3"]]
        with io.TextIOWrapper(pkg_resources.resource_stream("Mikado.tests", "trinity.ids")) as id_file:
            ids = [tuple(line.rstrip().split("\t")) for line in id_file]

        id_temp_file = tempfile.NamedTemporaryFile("wt", suffix=".txt")
        to_write = [ids[_] for _ in np.random.choice(len(ids), 10, replace=False)]
        [print(*idline, sep="\t", file=id_temp_file) for idline in to_write]
        id_temp_file.flush()

        for fname in files:
            with self.subTest(fname=fname):
                form = os.path.splitext(fname)[1]
                outfile = tempfile.NamedTemporaryFile("wt", suffix=form)
                outfile.close()
                self.assertFalse(os.path.exists(outfile.name))
                sys.argv = ["mikado", "util", "grep", id_temp_file.name, fname, outfile.name]
                pkg_resources.load_entry_point("Mikado", "console_scripts", "mikado")()
                self.assertTrue(os.path.exists(outfile.name))
                found = set()
                with to_gff(outfile.name, input_format=form[1:]) as stream:
                    for record in stream:
                        if record.is_transcript:
                            found.add(record.transcript)
                self.assertEqual(len(found), 10)
                self.assertEqual(found, set(_[0] for _ in to_write))

    @mark.slow
    def test_v_grep(self):
        files = [pkg_resources.resource_filename("Mikado.tests", fname) for fname in ["trinity.gtf", "trinity.gff3"]]
        with io.TextIOWrapper(pkg_resources.resource_stream("Mikado.tests", "trinity.ids")) as id_file:
            ids = [tuple(line.rstrip().split("\t")) for line in id_file]

        id_temp_file = tempfile.NamedTemporaryFile("wt", suffix=".txt")
        to_write = [ids[_] for _ in np.random.choice(len(ids), 10, replace=False)]
        others = [_ for _ in ids if _ not in to_write]
        [print(*idline, sep="\t", file=id_temp_file) for idline in to_write]
        id_temp_file.flush()

        for fname in files:
            with self.subTest(fname=fname):
                form = os.path.splitext(fname)[1]
                outfile = tempfile.NamedTemporaryFile("wt", suffix=form)
                outfile.close()
                self.assertFalse(os.path.exists(outfile.name))
                sys.argv = ["mikado", "util", "grep", "-v", id_temp_file.name, fname, outfile.name]
                pkg_resources.load_entry_point("Mikado", "console_scripts", "mikado")()
                self.assertTrue(os.path.exists(outfile.name))
                found = set()
                with to_gff(outfile.name, input_format=form[1:]) as stream:
                    for record in stream:
                        if record.is_transcript:
                            found.add(record.transcript)
                self.assertEqual(len(found), len(others))
                self.assertEqual(found, set(_[0] for _ in others))

    @mark.slow
    def test_problem_grep(self):
        fname = pkg_resources.resource_filename("Mikado.tests", "Chrysemys_picta_bellii_problematic.gff3")
        flist = pkg_resources.resource_filename("Mikado.tests", "Chrysemys_picta_bellii_problematic.list.txt")

        for flag in ("", "-v"):
            with self.subTest(flag=flag):
                form = os.path.splitext(fname)[1]
                outfile = tempfile.NamedTemporaryFile("wt", suffix=form)
                outfile.close()
                self.assertFalse(os.path.exists(outfile.name))
                if flag:
                    sys.argv = ["mikado", "util", "grep", flag, flist, fname, outfile.name]
                else:
                    sys.argv = ["mikado", "util", "grep", flist, fname, outfile.name]
                print(*sys.argv)
                pkg_resources.load_entry_point("Mikado", "console_scripts", "mikado")()
                self.assertTrue(os.path.exists(outfile.name))
                found = set()

                others = ["NC_023890.1:1..16875"]
                if flag != "-v":
                    for line in pkg_resources.resource_stream("Mikado.tests",
                                                              "Chrysemys_picta_bellii_problematic.list.txt"):
                        rec = line.decode().rstrip().split()[0]
                        print(line, rec)
                        others.append(rec)

                with to_gff(outfile.name, input_format=form[1:]) as stream:
                    for record in stream:
                        if record.feature in ("exon", "CDS"):
                            continue
                        if record.is_transcript:
                            found.add(record.transcript)
                        elif record.feature in ("pseudogene", "region"):
                            found.add(record.id)

                self.assertEqual(len(found), len(others))
                self.assertEqual(found, set(others))



if __name__ == "__main__":
    unittest.main()
