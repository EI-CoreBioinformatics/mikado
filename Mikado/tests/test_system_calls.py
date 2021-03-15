import collections
import csv
import glob
import gzip
import itertools
import copy
import logging
import os
import re
import numpy as np
import pandas as pd
import sys
import tempfile
import unittest
import pkg_resources
import pyfaidx
import toml
import yaml
from Mikado.configuration.configurator import load_and_validate_config
from Mikado.exceptions import InvalidConfiguration, InvalidParsingFormat
from Mikado.daijin.mikado import mikado_pipeline
from Mikado.daijin.assemble import assemble_transcripts_pipeline
from Mikado.configuration import print_config, DaijinConfiguration, MikadoConfiguration
import rapidjson as json
from Mikado.subprograms._utils import _set_pick_mode, check_log_settings_and_create_logger
from Mikado.subprograms.pick import _parse_regions
from Mikado.subprograms.serialise import serialise
from pytest import mark
from Mikado import configuration
from Mikado.subprograms import configure as sub_configure
from Mikado.subprograms.util.convert import launch as convert_launch, convert_parser
from Mikado.configuration import configurator, daijin_configurator
from Mikado.picking import picker
from Mikado.preparation import prepare
from Mikado.scales.compare import compare
from Mikado.scales.reference_preparation.indexing import load_index
from Mikado.scales.calculator import Calculator
from Mikado.subprograms.prepare import prepare_launcher, parse_gff_args
from Mikado.subprograms.prepare import setup as prepare_setup
from Mikado.utilities import to_region, IntervalTree
from Mikado.utilities.namespace import Namespace
from Mikado.utilities.log_utils import create_null_logger
from Mikado.parsers.GFF import GffLine
import sqlite3
import shutil
from Mikado.parsers import parser_factory
from Mikado.transcripts import Transcript
import threading
from time import sleep
import pysam
import io
try:
    from yaml import CSafeLoader as yLoader
except ImportError:
    from yaml import SafeLoader as yLoader


class ConvertCheck(unittest.TestCase):

    @mark.slow
    def test_convert_from_bam(self):

        bam_inp = pkg_resources.resource_filename("Mikado.tests", "test_mRNA.bam")
        for outp in ("gff3", "gtf", "bed12"):
            with self.subTest(outp=outp), tempfile.NamedTemporaryFile(mode="wt") as outfile:
                # sys.argv = ["", "util", "convert", "-of", outp, bam_inp, outfile.name]
                argv = ["-of", outp, bam_inp, outfile.name]
                parser = convert_parser()
                args = parser.parse_args(argv)
                convert_launch(args)
                # pkg_resources.load_entry_point("Mikado", "console_scripts", "mikado")()
                outfile.flush()
                self.assertGreater(os.stat(outfile.name).st_size, 0)
                lines = [_ for _ in open(outfile.name)]
                if outp == "gff3":
                    self.assertEqual(len(lines), 1826)
                elif outp == "bed12":
                    self.assertEqual(len(lines), 270)
                print(os.stat(outfile.name).st_size, len(lines))
                self.assertTrue(any(["TraesCS2B02G055500.1" in line for line in lines]))

    @mark.slow
    def test_convert_from_problematic(self):
        probl = pkg_resources.resource_filename("Mikado.tests", "Chrysemys_picta_bellii_problematic.gff3")
        for outp in ("gtf", "bed12"):
            with self.subTest(outp=outp), tempfile.NamedTemporaryFile(mode="wt") as outfile:
                sys.argv = ["", "util", "convert", "-of", outp, probl, outfile.name]
                pkg_resources.load_entry_point("Mikado", "console_scripts", "mikado")()
                self.assertGreater(os.stat(outfile.name).st_size, 0)
                lines = [_ for _ in open(outfile.name)]
                self.assertTrue(any(["rna-NC_023890.1:71..1039" in line for line in lines]))
                self.assertTrue(any(["rna-NC_023890.1:1040..1107" in line for line in lines]))
                self.assertTrue(any(["gene-LOC112059550" in line for line in lines]))
                self.assertTrue(any(["id-LOC112059311" in line for line in lines]))

    @mark.slow
    def test_generic_convert_check(self):
        # Test all possible combinations
        gtf = pkg_resources.resource_filename("Mikado.tests", "trinity.gtf")
        gff3 = pkg_resources.resource_filename("Mikado.tests", "trinity.gff3")
        bed12 = pkg_resources.resource_filename("Mikado.tests", "trinity.bed12")
        bam = pkg_resources.resource_filename("Mikado.tests", "trinity.bam")
        for inp_file, inp_format in [(gtf, "gtf"), (gff3, "gff3"), (bed12, "bed12"), (bam, "bam")]:
            for outp in ["gff3", "bed12", "gtf"]:
                outfile = tempfile.NamedTemporaryFile(suffix=f".{outp}")
                sys.argv = ["-of", outp, inp_file, outfile.name]
                args = convert_parser().parse_args(sys.argv)
                if outp == inp_format:
                    with self.assertRaises(SystemExit) as exit:
                        convert_launch(args)
                    self.assertEqual(exit.exception.code, 1)
                else:
                    convert_launch(args)
                    self.assertGreater(os.stat(outfile.name).st_size, 0)

    @mark.slow
    def test_unsorted_inputs(self):
        gtf = pkg_resources.resource_filename("Mikado.tests", "trinity_unsorted.gtf")
        gff3 = pkg_resources.resource_filename("Mikado.tests", "trinity_unsorted.gff3")
        args = Namespace(default=None)

        for inp_file, assume_sorted in itertools.product([gtf, gff3], [False, True]):
            if inp_file == gtf:
                outf = "gff3"
            else:
                outf = "gtf"
            with self.subTest(assume_sorted=assume_sorted), \
                    tempfile.NamedTemporaryFile(suffix=f".{outf}", mode="wt", delete=False) as out:
                args.assume_sorted = assume_sorted
                args.gf = inp_file
                args.out = out
                args.out_format = outf
                if assume_sorted:
                    with self.assertRaises(InvalidParsingFormat):
                        convert_launch(args=args)
                else:
                    convert_launch(args)
                    self.assertGreater(os.stat(out.name).st_size, 0)
                os.remove(out.name)


class PrepareCheck(unittest.TestCase):

    __genomefile__ = None

    @classmethod
    def setUpClass(cls):
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
                                                      ("tr_c58_g1_i8.mrna2", 383),
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

        cls._trinity_redundant = ["c58_g1_i8.mrna2", "c58_g1_i10.mrna2"]

        cls.maxDiff = None

    def setUp(self):

        self.conf = configurator.load_and_validate_config(None)
        self.conf.filename = "nofile"  # FIXME: Filename is added dynamically to this obj and not validated
        self.conf.seed = 1066
        self.conf.reference.genome = self.fai.filename.decode()
        assert isinstance(self.conf.reference.genome, str)
        self.logger = create_null_logger("prepare")
        self.conf.prepare.exclude_redundant = False

    def tearDown(self):
        logging.shutdown()

    def test_parsing_gffs_from_cli(self):
        """"""
        # parse_gff_args

        args = Namespace(default=None)
        args.gff = ["foo.gff", "foo.gff", "bar.gff"]
        with self.assertRaises(InvalidConfiguration) as exc:
            parse_gff_args(self.conf, args)
        self.assertTrue(re.search(r"Repeated elements among the input GFFs", str(exc.exception)))
        # Now test that
        args.gff = ["foo.gff", "bar.gff"]
        args.strand_specific_assemblies = "foo.gff,bar.gff,absent.gff"
        with self.assertRaises(InvalidConfiguration) as exc:
            parse_gff_args(self.conf, args)
        self.assertTrue(re.search(r"Incorrect number of strand-specific assemblies specified!", str(exc.exception)),
                        str(exc.exception))
        args.strand_specific_assemblies = "foo.gff,absent.gff"
        with self.assertRaises(InvalidConfiguration) as exc:
            parse_gff_args(self.conf, args)
        self.assertTrue(re.search(r"Incorrect assembly file specified as strand-specific", str(exc.exception)),
                        str(exc.exception))
        args.strand_specific_assemblies = "foo.gff"
        args.labels = "foo,"
        with self.assertRaises(InvalidConfiguration) as exc:
            parse_gff_args(self.conf, args)
        self.assertTrue(re.search(r"Empty labels provided!", str(exc.exception)), str(exc.exception))
        args.labels = "foo,foo"
        with self.assertRaises(InvalidConfiguration) as exc:
            parse_gff_args(self.conf, args)
        self.assertTrue(re.search(r"Duplicated labels detected", str(exc.exception)), str(exc.exception))
        args.labels = "foo"
        with self.assertRaises(InvalidConfiguration) as exc:
            parse_gff_args(self.conf, args)
        self.assertTrue(re.search(r"Incorrect number of labels specified", str(exc.exception)), str(exc.exception))
        args.labels = None
        conf = parse_gff_args(self.conf, args)
        self.assertEqual(conf.prepare.files.labels, ["1", "2"])
        self.assertEqual(conf.prepare.files.exclude_redundant, [False, False])
        self.assertEqual(conf.prepare.files.reference, [False, False])
        self.assertEqual(conf.prepare.files.strand_specific_assemblies, ["foo.gff"])

    def test_parse_prepare_options(self):
        """Tests for setting the values correctly from the command line argument parser."""

    @mark.slow
    def test_varying_max_intron(self):

        self.conf.prepare.files.labels.append("tr")
        dir = tempfile.TemporaryDirectory(prefix="test_varying_max_intron")
        self.conf.prepare.files.output_dir = dir.name
        args = Namespace()
        args.configuration = self.conf
        test_file = "trinity.gtf"
        self.conf.prepare.files.gff = [pkg_resources.resource_filename("Mikado.tests",
                                                                                test_file)]
        self.conf.prepare.files.strip_cds = [False]

        for max_intron in (20, 200, 1000, 5000):
            with self.subTest(max_intron=max_intron):
                self.conf.prepare.max_intron_length = max_intron
                prepare.prepare(args.configuration, self.logger)
                gtf = os.path.join(self.conf.prepare.files.output_dir, "mikado_prepared.gtf")
                self.assertGreater(os.stat(gtf).st_size, 0, test_file)
                transcripts = dict()
                for row in parser_factory(gtf):
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

        self.conf.prepare.files.labels.append("tr")
        dir = tempfile.TemporaryDirectory(prefix="test_prepare_trinity_gff")
        self.conf.prepare.files.output_dir = dir.name
        args = Namespace()
        args.configuration = self.conf
        args.procs = 1
        args.single_thread = True
        args.seed = 10
        self.conf.seed = 10

        for test_file in ("trinity.gff3",
                          "trinity.match_matchpart.gff3",
                          "trinity.cDNA_match.gff3",
                          "trinity.gtf",
                          "trinity.no_transcript_feature.gtf"):
            with self.subTest(test_file=test_file):
                self.conf.prepare.files.strip_cds = [False]
                self.conf.prepare.files.gff = [pkg_resources.resource_filename("Mikado.tests",
                                                                                        test_file)]

                with self.assertLogs(self.logger) as cm:
                    try:
                        prepare.prepare(args.configuration, self.logger)
                    except OSError:
                        raise OSError(cm.output)

                # Now that the program has run, let's check the output
                fasta = os.path.join(self.conf.prepare.files.output_dir,
                                                "mikado_prepared.fasta")
                self.assertGreater(os.stat(fasta).st_size, 0, (test_file, cm.output))
                fa = pyfaidx.Fasta(fasta)
                res = dict((_, len(fa[_])) for _ in fa.keys())
                fa.close()
                check = dict(_ for _ in self.trinity_res.items())
                check.pop(self.conf.prepare.files.labels[0] + "_" + self._trinity_redundant[0])
                self.assertEqual(res, check)
                os.remove(os.path.join(self.conf.prepare.files.output_dir,
                                       "mikado_prepared.fasta.fai"))
        dir.cleanup()

    @mark.slow
    def test_prepare_trinity_and_cufflinks(self):

        self.conf.prepare.files.labels = ["cl", "tr"]

        self.conf.prepare.files.gff = [None, None]
        self.conf.prepare.files.strip_cds = [False, True]
        self.conf.prepare.files.out_fasta = "mikado_prepared.fasta"
        self.conf.prepare.files.out = "mikado_prepared.gtf"
        self.conf.multiprocessing_method = "spawn"
        args = Namespace()

        with tempfile.TemporaryDirectory(prefix="test_prepare_trinity_and_cufflinks") as folder:
            self.conf.prepare.files.output_dir = folder
            for cuff_file, test_file in itertools.product(
                    ("cufflinks.gtf", "cufflinks.no_transcript.gtf"),
                    (("trinity.gff3", "trinity.match_matchpart.gff3", "trinity.cDNA_match.gff3", "trinity.gtf",
                      "trinity.no_transcript_feature.gtf", "trinity.bam"))):
                for proc in (1, 3):
                    with self.subTest(test_file=test_file, cuff_file=cuff_file, proc=proc):
                        self.conf.prepare.files.gff[0] = pkg_resources.resource_filename("Mikado.tests", cuff_file)
                        self.conf.prepare.files.gff[1] = pkg_resources.resource_filename("Mikado.tests", test_file)
                        self.conf.prepare.files.out_fasta = "mikado_prepared.fasta"
                        self.conf.prepare.files.out = "mikado_prepared.gtf"
                        args.configuration = self.conf
                        args.configuration.seed = 10
                        args.configuration.threads = 1  # proc
                        args.configuration.prepare.exclude_redundant = False
                        args.configuration.prepare.strip_cds = True
                        prepare.prepare(args.configuration, self.logger)

                        # Now that the program has run, let's check the output
                        self.assertTrue(os.path.exists(os.path.join(self.conf.prepare.files.output_dir,
                                                                    "mikado_prepared.fasta")))
                        self.assertGreater(os.stat(os.path.join(self.conf.prepare.files.output_dir,
                                                                "mikado_prepared.fasta")).st_size, 0)

                        fa = pysam.FastaFile(os.path.join(self.conf.prepare.files.output_dir, "mikado_prepared.fasta"))
                        res = dict((name, length) for name, length in zip(fa.references, fa.lengths))
                        fa.close()
                        precal = self.trinity_res.copy()
                        precal.update(self.cuff_results)
                        precal.pop("tr_" + self._trinity_redundant[0])
                        self.assertEqual(res, precal)

    @mark.slow
    def test_prepare_with_bam(self):
        self.conf.prepare.files.labels = ["cl", "pb"]

        self.conf.prepare.files.gff = [pkg_resources.resource_filename("Mikado.tests", "cufflinks.gtf"),
                                       pkg_resources.resource_filename("Mikado.tests", "pacbio.bam")]
        self.conf.prepare.files.strip_cds = [False, False]
        self.conf.prepare.files.out_fasta = "mikado_prepared.fasta"
        self.conf.prepare.files.out = "mikado_prepared.gtf"
        self.conf.prepare.files.source_score["cl"] = 0
        self.conf.prepare.files.source_score["pb"] = 5
        args = Namespace()

        with tempfile.TemporaryDirectory(prefix="test_prepare_cufflinks_and_bam") as folder:
            for proc in (1, 2):
                self.conf.prepare.files.out_fasta = "mikado_prepared.fasta"
                self.conf.prepare.files.out = "mikado_prepared.gtf"
                args.configuration = self.conf
                args.configuration.seed = 10
                args.configuration.threads = proc
                args.configuration.prepare.exclude_redundant = False
                args.configuration.prepare.strip_cds = True
                args.configuration.prepare.files.output_dir = folder
                prepare.prepare(args.configuration, self.logger)
                # Now that the program has run, let's check the output
                self.assertTrue(os.path.exists(os.path.join(self.conf.prepare.files.output_dir,
                                                            "mikado_prepared.fasta")))
                self.assertGreater(os.stat(os.path.join(self.conf.prepare.files.output_dir,
                                                        "mikado_prepared.fasta")).st_size, 0)

                fa = pysam.FastaFile(os.path.join(self.conf.prepare.files.output_dir, "mikado_prepared.fasta"))
                self.assertTrue(any([_.startswith("pb") for _ in fa.references]))
                self.assertTrue(any([_.startswith("cl") for _ in fa.references]))
                fa.close()

    @mark.slow
    def test_prepare_with_cds(self):

        rev_strand = {"+": "-", "-": "+"}

        self.conf.prepare.files.labels = ["ann"]
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

        self.conf.prepare.files.gff = []
        self.conf.prepare.files.output_dir = folder.name
        self.conf.prepare.files.out_fasta = "mikado_prepared.fasta"
        self.conf.prepare.files.out = "mikado_prepared.gtf"
        args = Namespace()
        args.configuration = self.conf

        for fname in [ann_gff3, rev_ann_gff3.name]:
            for strip, strip_faulty_cds in itertools.product((True, False), (True, False)):
                with self.subTest(fname=fname, strip=strip, strip_faulty_cds=strip_faulty_cds):
                    self.conf.prepare.files.gff = [fname]
                    args.configuration.prepare.strip_cds = False
                    args.configuration.prepare.single = False
                    self.conf.prepare.files.strip_cds = [strip]
                    self.conf.prepare.strip_faulty_cds = strip_faulty_cds
                    with self.assertLogs(self.logger, "INFO") as cm:
                        try:
                            prepare.prepare(args.configuration, logger=self.logger)
                        except SystemExit:
                            raise SystemExit("\n".join(cm.output))
                    fasta = os.path.join(self.conf.prepare.files.output_dir, "mikado_prepared.fasta")
                    self.assertTrue(os.path.exists(fasta), "\n".join(cm.output))
                    if fname == ann_gff3 or strip is True or strip_faulty_cds is True:
                        self.assertGreater(os.stat(fasta).st_size, 0, "\n".join(cm.output))
                        fa = pyfaidx.Fasta(fasta)
                        self.assertEqual(len(fa.keys()), 2, "\n".join(cm.output))
                        fa.close()
                    else:
                        self.assertEqual(os.stat(fasta).st_size, 0,
                                         str(strip) + " " + str(strip_faulty_cds) + " " + \
                                         str(rev_ann_gff3.name == fname) + "\n" + "\n".join(cm.output))

                    # Now verify that no model has CDS
                    gtf = os.path.join(self.conf.prepare.files.output_dir, "mikado_prepared.gtf")
                    models = dict()
                    with parser_factory(gtf) as file_gtf:
                        for line in file_gtf:
                            if line.header:
                                continue
                            elif line.is_transcript:
                                models[line.id] = Transcript(line)
                            else:
                                models[line.parent[0]].add_exon(line)
                    [models[model].finalize() for model in models]
                    for model in models:
                        if strip is False and rev_ann_gff3.name != fname and strip_faulty_cds is True:
                            self.assertTrue(models[model].is_coding,
                                                (fname, strip_faulty_cds, strip,
                                                models[model].format("gtf"))
                                            )
                        elif strip is True:
                            self.assertFalse(models[model].is_coding, (
                                fname, strip_faulty_cds, strip, models[model].format("gtf")))
        rev_ann_gff3.close()

    @mark.slow
    def test_cdna_redundant_cds_not(self):
        """This test will verify whether the new behaviour of not considering redundant two models with same
        exon structure but different CDS does function properly."""

        gtf = pkg_resources.resource_filename("Mikado.tests", "cds_test_1.gtf")
        self.conf.prepare.files.gff = [gtf]
        self.conf.prepare.files.labels = [""]
        self.conf.prepare.files.strip_cds = [False]
        self.conf.prepare.files.out_fasta = "mikado_prepared.fasta"
        self.conf.prepare.files.out = "mikado_prepared.gtf"
        self.conf.prepare.files.log = "prepare.log"
        self.conf.prepare.strip_cds = False

        args = Namespace(default=None)
        args.strip_cds = False
        args.configuration = self.conf.copy()
        del args.configuration.prepare.files.output_dir
        args.log = None
        # Keep them *OUTSIDE* otherwise weird disk errors might happen
        with tempfile.TemporaryDirectory() as folder_a, tempfile.TemporaryDirectory() as folder_b:
            for b, folder in [(True, folder_a), (False, folder_b)]:
                os.makedirs(folder, exist_ok=True)
                with self.subTest(b=b):
                    os.makedirs(folder, exist_ok=True)
                    self.conf.prepare.files.exclude_redundant = [b]
                    args.output_dir = folder
                    args.seed = 10
                    args.procs = 1
                    args.list = None
                    args.gffs = None
                    args.strand_specific_assemblies = []
                    args.labels = None
                    args.configuration = self.conf
                    args.exclude_redundant = b
                    args.out, args.out_fasta = None, None
                    args.configuration.prepare.files.log = "prepare.log"
                    if isinstance(args.configuration.reference.genome, bytes):
                        args.configuration.reference.genome = args.configuration.reference.genome.decode()
                    args.log = "prepare.log"
                    self.logger.setLevel("DEBUG")
                    assert os.path.exists(folder)
                    self.assertEqual(args.strand_specific_assemblies, [])
                    self.assertEqual(args.configuration.prepare.files.strand_specific_assemblies, [])
                    args, mikado_configuration, _logger = prepare_setup(args)
                    self.assertIsNotNone(mikado_configuration)
                    # self.assertEqual(args.output_dir, folder)
                    self.assertEqual(mikado_configuration.prepare.files.output_dir, folder)
                    self.assertIn(os.path.dirname(mikado_configuration.prepare.files.out_fasta),
                                  (folder, ""), mikado_configuration)
                    self.assertIn(os.path.dirname(mikado_configuration.prepare.files.out),
                                  (folder, ""), mikado_configuration)

                    prepare.prepare(mikado_configuration, _logger)
                    self.assertTrue(os.path.exists(folder))
                    self.assertTrue(os.path.isdir(folder))
                    self.assertTrue(os.path.exists(os.path.join(folder,
                                                                "mikado_prepared.fasta")),
                                    open(os.path.join(folder,
                                                      "prepare.log")).read())
                    fa = pyfaidx.Fasta(os.path.join(folder,
                                                    "mikado_prepared.fasta"))
                    logged = [_ for _ in open(os.path.join(mikado_configuration.prepare.files.output_dir,
                                                           mikado_configuration.prepare.files.log))]
                    self.assertFalse("AT5G01530.1" in fa.keys(), (b, sorted(list(fa.keys())),
                                                                  print(*logged, sep="\n", file=open("/tmp/log", "wt"))))
                    self.assertTrue("AT5G01530.2" in fa.keys(), print(*logged, sep="\n", file=open("/tmp/log", "wt")))
                    if b is False:
                        self.assertEqual(len(fa.keys()), 4)
                        self.assertEqual(sorted(fa.keys()), sorted(["AT5G01530."+str(_) for _ in [0, 2, 3, 4]]))
                    else:
                        self.assertEqual(len(fa.keys()), 3, (fa.keys(), logged))
                        self.assertIn("AT5G01530.0", fa.keys())
                        self.assertIn("AT5G01530.2", fa.keys())
                        self.assertNotIn("AT5G01530.3", fa.keys())
                        self.assertIn("AT5G01530.4", fa.keys())
                    gtf_file = os.path.join(folder, "mikado_prepared.gtf")
                    fa.close()
                    coding_count = 0
                    with parser_factory(gtf_file) as gtf:
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
        self.conf.prepare.files.gff = [gtf]
        self.conf.prepare.files.labels = [""]
        self.conf.prepare.files.out_fasta = "mikado_prepared.fasta"
        self.conf.prepare.files.strip_cds = [False]
        self.conf.prepare.files.out = "mikado_prepared.gtf"
        self.conf.prepare.files.log = "prepare.log"
        self.conf.prepare.strip_cds = False
        self.conf.prepare.minimum_cdna_length = 150  # Necessary for testing A5

        with tempfile.TemporaryDirectory(prefix="test_negative_cdna_redundant_cds_not") as folder:
            self.conf.prepare.files.output_dir = folder
            args = Namespace()
            args.log_level = None
            args.strip_cds = False
            args.configuration = self.conf
            args.reference = None
            for b in (False, True):
                with self.subTest(b=b):
                    self.conf.prepare.files.exclude_redundant = [b]
                    args.configuration = self.conf
                    args.configuration.seed = 10
                    args.exclude_redundant = b
                    args.output_dir = folder
                    args.out_fasta = None
                    args.log = None
                    args.gff = None
                    args.list = None
                    args.strand_specific_assemblies = None
                    if isinstance(args.configuration.reference.genome, bytes):
                        args.configuration.reference.genome = args.configuration.reference.genome.decode()
                    args, mikado_config, _ = prepare_setup(args)
                    assert hasattr(mikado_config.reference, "genome"), args.configuration.reference
                    prepare.prepare(mikado_config, self.logger)
                    self.assertTrue(os.path.exists(os.path.join(self.conf.prepare.files.output_dir,
                                                                "mikado_prepared.fasta")))
                    fa = pyfaidx.Fasta(os.path.join(self.conf.prepare.files.output_dir,
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

                    gtf_file = os.path.join(self.conf.prepare.files.output_dir, "mikado_prepared.gtf")

                    coding_count = 0
                    with parser_factory(gtf_file) as gtf:
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

    @mark.slow
    def test_truncated_cds(self):
        files = ["test_truncated_cds.gff3"]
        files = [pkg_resources.resource_filename("Mikado.tests", filename) for filename in files]
        self.conf.prepare.files.gff = files
        self.conf.prepare.files.strip_cds = [False] * len(files)
        self.conf.prepare.files.labels = [""]
        dir = tempfile.TemporaryDirectory(prefix="test_truncated_cds")
        self.conf.prepare.files.output_dir = dir.name
        self.conf.prepare.files.out_fasta = "mikado_prepared.fasta"
        self.conf.prepare.files.out = "mikado_prepared.gtf"
        self.conf.prepare.strip_cds = False

        self.conf.reference.genome = pkg_resources.resource_filename("Mikado.tests",
                                                                           "test_truncated_cds.fa")
        args = Namespace()
        args.strip_cds = False
        args.configuration = self.conf
        prepare.prepare(args.configuration, self.logger)
        self.assertTrue(os.path.exists(os.path.join(self.conf.prepare.files.output_dir,
                                                    "mikado_prepared.fasta")))
        fa = pyfaidx.Fasta(os.path.join(self.conf.prepare.files.output_dir,
                                        "mikado_prepared.fasta"))
        self.assertEqual(len(fa.keys()), 1)
        gtf_file = os.path.join(self.conf.prepare.files.output_dir, "mikado_prepared.gtf")
        with parser_factory(gtf_file) as gtf_handle:
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
        self.conf.prepare.files.output_dir = dir.name
        self.conf.prepare.files.out_fasta = "mikado_prepared.fasta"
        self.conf.prepare.files.out = "mikado_prepared.gtf"
        self.conf.prepare.strip_cds = False
        self.conf.prepare.exclude_redundant = True
        self.conf.threads = 1

        self.conf.reference.genome = self.fai.filename.decode()

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

        self.conf.prepare.files.gff = [t_file.name, t2_file.name]
        self.conf.prepare.files.labels = ["T1", "T2"]
        self.conf.prepare.files.strip_cds = [True, True]

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
                self.conf.prepare.files.source_score = {"T1": t1, "T2": t2}
                args = Namespace()
                args.strip_cds = False
                args.configuration = self.conf
                prepare.prepare(args.configuration, self.logger)
                self.assertGreater(os.stat(self.conf.prepare.files.out_fasta).st_size, 0)
                fa = pyfaidx.Fasta(self.conf.prepare.files.out_fasta)
                self.assertEqual(len(fa.keys()), 1, round)
                if res != "rand":
                    key = list(fa.keys())[0]
                    self.assertEqual(key, res, round)

    @mark.slow
    def test_reference_selection(self):

        dir = tempfile.TemporaryDirectory(prefix="test_reference_selection")
        self.conf.prepare.files.output_dir = outdir = dir.name
        self.conf.prepare.files.out_fasta = "mikado_prepared.fasta"
        self.conf.prepare.files.out = "mikado_prepared.gtf"
        self.conf.prepare.files.strip_cds = [False, False]
        self.conf.prepare.files.exclude_redundant = [False, False]
        self.conf.prepare.strip_cds = False
        self.conf.prepare.exclude_redundant = True

        self.conf.reference.genome = self.fai.filename.decode()

        t = Transcript()
        # This is *key*. Transcript T1 should never be selected, unless we are having a "lenient" analysis.
        # However, by flagging T1 as a "reference" transcript, we should retain it
        self.conf.prepare.lenient = False
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
            print(t.format(fformat, transcriptomic=False), file=t_file)
            t_file.flush()
            print(t2.format(fformat, transcriptomic=False), file=t2_file)
            t2_file.flush()

            self.conf.prepare.files.gff = [t_file.name, t2_file.name]
            self.conf.prepare.files.strand_specific_assemblies = [t_file.name, t2_file.name]
            self.conf.prepare.files.labels = ["T1", "T2"]

            for round in rounds:
                for iteration in range(2):  # Repeat each test 2 times, not more for time length reasons
                    with self.subTest(round=round, format=format, iteration=iteration,
                                      msg="Starting round {} ({})".format(round, rounds[round])):
                        t1_score, t2_score, is_ref, res, corr_strand = rounds[round]
                        self.conf.threads = 1
                        self.conf.prepare.files.source_score = {"T1": t1_score,
                                                                         "T2": t2_score}
                        self.conf.prepare.files.reference = []
                        for label in is_ref:
                            if label == "T1":
                                self.conf.prepare.files.reference.append(t_file.name)
                            elif label == "T2":
                                self.conf.prepare.files.reference.append(t2_file.name)

                        args = Namespace()
                        args.strip_cds = False
                        args.configuration = self.conf
                        with self.assertLogs(self.logger, "INFO") as cm:
                            try:
                                prepare.prepare(args.configuration, self.logger)
                            except SystemExit:
                                print(t.format(fformat))
                                print()
                                print(t2.format(fformat))
                                print()
                                print(*cm.output, sep="\n")
                                raise
                        self.assertGreater(os.stat(self.conf.prepare.files.out_fasta).st_size, 0,
                                           (round, fformat, iteration, "\n".join(cm.output)))
                        fa = pyfaidx.Fasta(os.path.join(outdir, os.path.basename(
                            self.conf.prepare.files.out_fasta)))
                        if res != "rand":
                            key = sorted(list(fa.keys()))
                            self.assertEqual(key, res, round)
                        else:
                            self.assertEqual(len(fa.keys()), 1, (round, fa.keys(), res))
                        gtf = os.path.join(outdir, os.path.basename(
                            self.conf.prepare.files.out))
                        strand = [_.strand for _ in parser_factory(gtf)]
                        self.assertEqual(len(set(strand)), 1, strand)
                        self.assertEqual(set(strand).pop(), corr_strand,
                                         (round, self.conf.prepare.files.reference,
                                          set(strand), corr_strand))

    @mark.slow
    def test_reference_cds_kept(self):

        t = Transcript()
        # This is *key*. Transcript T1 should never be selected, unless we are having a "lenient" analysis.
        # However, by flagging T1 as a "reference" transcript, we should retain it
        self.conf.prepare.lenient = False
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
        self.conf.prepare.files.output_dir = outdir = dir.name
        self.conf.prepare.files.out_fasta = "mikado_prepared.fasta"
        self.conf.prepare.files.out = "mikado_prepared.gtf"
        self.conf.prepare.strip_cds = True
        self.conf.prepare.exclude_redundant = True
        self.conf.reference.genome = self.fai.filename.decode()

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

            self.conf.prepare.files.gff = [t_file.name, t2_file.name]
            self.conf.prepare.files.strand_specific_assemblies = [t_file.name, t2_file.name]
            self.conf.prepare.files.labels = ["T1", "T2"]
            self.conf.prepare.files.exclude_redundant = [True, True]
            self.conf.prepare.files.strip_cds = [True, True]

            for round in rounds:
                for iteration in range(2):  # Repeat each test 2 times for time length reasons
                    with self.subTest(round=round, format=format, iteration=iteration,
                                      msg="Starting round {} ({})".format(round, rounds[round])):
                        t1_score, t2_score, is_ref, res, coding = rounds[round]
                        self.conf.prepare.files.source_score = {"T1": t1_score,
                                                                         "T2": t2_score}
                        self.conf.prepare.files.reference = []
                        self.conf.threads = 1
                        for label in is_ref:
                            if label == "T1":
                                self.conf.prepare.files.reference.append(t_file.name)
                            elif label == "T2":
                                self.conf.prepare.files.reference.append(t2_file.name)

                        args = Namespace()
                        args.strip_cds = False
                        args.configuration = self.conf
                        with self.assertLogs(self.logger, "INFO") as cm:
                            prepare.prepare(args.configuration, self.logger)
                        self.assertGreater(os.stat(self.conf.prepare.files.out_fasta).st_size, 0,
                                           (round, fformat, iteration, "\n".join(cm.output)))
                        fa = pyfaidx.Fasta(os.path.join(outdir, os.path.basename(
                            self.conf.prepare.files.out_fasta)))
                        if res != "rand":
                            key = sorted(list(fa.keys()))
                            self.assertEqual(key, res, (round, rounds[round], cm.output))
                        else:
                            self.assertEqual(len(fa.keys()), 1, (round, fa.keys(), res))
                        gtf = os.path.join(outdir, os.path.basename(
                            self.conf.prepare.files.out))
                        with_cds = set()
                        for line in parser_factory(gtf):
                            if line.feature == "CDS":
                                with_cds.add(line.transcript)
                            else:
                                continue
                        self.assertEqual(coding, with_cds)


@mark.slow
class CompareCheck(unittest.TestCase):

    """Test to check that compare interacts correctly with match, match_part, cDNA_match"""

    def test_index(self):

        # Create the list of files
        files = ["trinity.gtf",
                 "trinity.gff3",
                 "trinity.cDNA_match.gff3",
                 "trinity.match_matchpart.gff3",
                 "trinity.bed12"]

        namespace = Namespace(default=False)
        namespace.distance = 2000
        namespace.index = True
        namespace.prediction = None
        namespace.log = None
        logger = create_null_logger("null")
        for ref in files:
            with self.subTest(ref=ref), tempfile.TemporaryDirectory() as folder:
                temp_ref = os.path.join(folder, ref)
                os.symlink(pkg_resources.resource_filename("Mikado.tests", ref), temp_ref)
                namespace.reference = parser_factory(temp_ref)
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
            with self.subTest(ref=ref, pred=pred), tempfile.TemporaryDirectory(
                    prefix="test_compare_trinity_{}_{}".format(os.path.splitext(ref)[-1],
                                                               os.path.splitext(pred)[-1])) as folder:
                namespace.reference = parser_factory(ref)
                namespace.prediction = parser_factory(pred)
                namespace.processes = 1
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

    def test_compare_problematic(self):

        problematic = pkg_resources.resource_filename("Mikado.tests", "Chrysemys_picta_bellii_problematic.gff3")
        namespace = Namespace(default=False)
        namespace.distance = 2000
        namespace.no_save_index = True
        namespace.protein_coding = False
        namespace.exclude_utr = False
        namespace.self = False
        namespace.gzip = False

        for proc in (1, 3):
            with tempfile.TemporaryDirectory(prefix="test_compare_problematic_{}_".format(proc)) as folder:
                namespace.reference = parser_factory(problematic)
                namespace.prediction = parser_factory(problematic)
                namespace.processes = proc
                namespace.log = None
                namespace.out = os.path.join(folder, "compare_problematic")
                self.assertTrue(os.path.exists(folder))
                self.assertIsInstance(folder, str)
                with self.assertLogs("main_compare", level="INFO") as cm:
                    compare(namespace)
                refmap = "{}.refmap".format(namespace.out)
                tmap = "{}.tmap".format(namespace.out)
                stats = "{}.stats".format(namespace.out)

                log = [str(_).rstrip() for _ in cm.records]
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

    @mark.slow
    def test_internal(self):
        ref = pkg_resources.resource_filename("Mikado.tests", "reference.gff3")
        namespace = Namespace(default=False)
        namespace.distance = 2000
        namespace.no_save_index = True
        namespace.protein_coding = False
        namespace.exclude_utr = False
        namespace.self = False
        namespace.gzip = False
        namespace.internal = True
        for proc in (1, 3):
            with tempfile.TemporaryDirectory(prefix="test_internal_{}_".format(proc)) as folder, self.subTest(
                    proc=proc):
                namespace.reference = parser_factory(ref)
                namespace.prediction = parser_factory(ref)
                namespace.processes = proc
                namespace.log = None
                namespace.out = os.path.join(folder, "compare_internal")
                self.assertTrue(os.path.exists(folder))
                self.assertIsInstance(folder, str)
                with self.assertLogs("main_compare", level="INFO") as cm:
                    compare(namespace)
                refmap = "{}.refmap".format(namespace.out)
                tmap = "{}.tmap".format(namespace.out)
                stats = "{}.stats".format(namespace.out)
                # Stats and refmap should not exist
                log = [str(_).rstrip() for _ in cm.records]
                self.assertFalse(os.path.exists(refmap))
                self.assertFalse(os.path.exists(stats))
                self.assertTrue(os.path.exists(tmap))
                self.assertGreater(os.stat(tmap).st_size, 0, "\n".join(log))

                with open(tmap) as _:
                    reader = csv.DictReader(_, delimiter="\t")
                    for counter, line in enumerate(reader, start=1):
                        self.assertNotIn(line["ccode"], ("=", "_"))
                        self.assertEqual(line["ref_gene"], line["gid"])

    def test_self(self):
        ref = pkg_resources.resource_filename("Mikado.tests", "reference.gff3")
        namespace = Namespace(default=False)
        namespace.distance = 2000
        namespace.no_save_index = True
        namespace.protein_coding = False
        namespace.exclude_utr = False
        namespace.self = True
        namespace.gzip = False
        namespace.internal = False
        transcripts = [_ for _ in parser_factory(ref) if _.is_transcript is True]
        for proc in (1, 3):
            with tempfile.TemporaryDirectory(prefix="test_self_{}_".format(proc)) as folder, self.subTest(
                    proc=proc):
                namespace.reference = parser_factory(ref)
                namespace.prediction = parser_factory(ref)
                namespace.processes = proc
                namespace.log = None
                namespace.out = os.path.join(folder, "compare_internal")
                self.assertTrue(os.path.exists(folder))
                self.assertIsInstance(folder, str)
                with self.assertLogs("main_compare", level="INFO") as cm:
                    compare(namespace)
                refmap = "{}.refmap".format(namespace.out)
                tmap = "{}.tmap".format(namespace.out)
                stats = "{}.stats".format(namespace.out)
                # Stats should not exist
                log = [str(_).rstrip() for _ in cm.records]
                self.assertFalse(os.path.exists(stats))
                self.assertTrue(os.path.exists(refmap))
                self.assertGreater(os.stat(refmap).st_size, 0, "\n".join(log))
                self.assertTrue(os.path.exists(tmap))
                self.assertGreater(os.stat(tmap).st_size, 0, "\n".join(log))

                with open(refmap) as _:
                    reader = csv.DictReader(_, delimiter="\t")
                    for counter, line in enumerate(reader, start=1):
                        self.assertNotIn(line["ccode"], ("=", "_"))

                self.assertEqual(counter, len(transcripts))

                with open(tmap) as _:
                    reader = csv.DictReader(_, delimiter="\t")
                    for counter, line in enumerate(reader, start=1):
                        self.assertNotIn(line["ccode"], ("=", "_"))

                self.assertEqual(counter, len(transcripts))

class CompareFusionCheck(unittest.TestCase):

    """WARNING!!! This test seems to cause a deadlock if put in the same class as other compare tests,
    e.g. CompareCheck above. KEEP SEPARATE"""

    # @unittest.skip
    def test_false_fusion(self):

        """
        System test to verify that the false fusion is not called.
        WARNING: this test is quite brittle as I am creating the namespace myself
        instead of loading it from the subprograms.compare module.
        :return:
        """

        ref_file = pkg_resources.resource_filename("Mikado.tests", os.path.join("fusion_test",
                                                                                "fusion_test_ref.gff3"))
        pred_file = pkg_resources.resource_filename("Mikado.tests", os.path.join("fusion_test",
                                                                                 "fusion_test_pred.gtf"))
        with tempfile.TemporaryDirectory() as out, \
                parser_factory(ref_file) as reference, parser_factory(pred_file) as prediction:
            args = Namespace()
            args.no_save_index = True
            args.reference = reference
            args.prediction = prediction
            args.log = None
            args.out = os.path.join(out, "fusion_test", "fusion_test")
            args.distance = 2000
            args.verbose = True
            args.exclude_utr = False
            args.protein_coding = False
            args.index = False
            args.self = False
            args.extended_refmap = False
            args.gzip = False
            args.processes = 1
            with self.assertLogs("main_compare") as cmo:
                t = threading.Thread(target=compare, args=(args,))
                t.start()
                t.join(timeout=4)

            out_refmap = os.path.join(out, "fusion_test", "fusion_test.refmap")
            self.assertTrue(os.path.exists(out_refmap), cmo.output)
            self.assertGreater(os.stat(out_refmap).st_size, 0, cmo.output)
            with open(out_refmap) as refmap:
                for line in csv.DictReader(refmap, delimiter="\t"):
                    if line["ref_id"] not in ("AT1G78880.1", "AT1G78882.1"):
                        continue
                    self.assertEqual(line["ccode"], "=", (line, cmo.output))

        return


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
        namespace.reference = self.fai.filename.decode()
        namespace.external = None
        namespace.mode = ["permissive"]
        namespace.threads = 1
        namespace.blast_targets = []
        namespace.junctions = []
        namespace.new_scoring = None
        namespace.seed = 0
        namespace.min_clustering_cds_overlap = 0.2
        namespace.min_clustering_cdna_overlap = 0.2
        namespace.no_files = True
        namespace.scheduler = ""
        namespace.mode = ["stringent"]
        namespace.blast_chunks = 10
        namespace.exclude_redundant = False
        folder = tempfile.TemporaryDirectory()

        for output_format in ("json", "yaml", "toml", "config"):
            with self.subTest(output_format=output_format):
                out = os.path.join(folder.name, "configuration.{}".format(output_format))
                with open(out, "w") as out_handle:
                    namespace.out = out_handle
                    sub_configure.create_config(namespace)
                self.assertGreater(os.stat(out).st_size, 0)
                if output_format in ("toml", "config"):
                    raw = toml.load(open(out))
                elif output_format == "yaml":
                    raw = yaml.load(open(out), Loader=yaml.SafeLoader)
                else:
                    raw = json.load(open(out))
                self.assertIn("threads", raw)
                loaded_raw = MikadoConfiguration.Schema().load(raw)
                conf = configuration.configurator.load_and_validate_config(out)
                conf = configuration.configurator.check_and_load_scoring(conf)
                conf = configuration.configurator.check_and_load_scoring(conf)
        folder.cleanup()

    def test_external_file(self):
        ext_file = pkg_resources.resource_filename("Mikado.tests", "external.yaml")
        for external in (True, False):
            with self.subTest(external=external):
                if external == False:
                    with self.assertRaises(InvalidConfiguration):
                        load_and_validate_config(ext_file, external=external)
                else:
                    config = load_and_validate_config(ext_file, external=external)
                    self.assertIsInstance(config, (MikadoConfiguration, DaijinConfiguration))
                    self.assertEqual(config.pick.alternative_splicing.max_isoforms, 15)

    def test_seed(self):
        namespace = Namespace(default=False)
        namespace.scoring = None
        namespace.intron_range = None
        namespace.reference = self.fai.filename.decode()
        namespace.external = None
        namespace.threads = 1
        namespace.multiprocessing_method = "spawn"
        namespace.blast_targets = []
        namespace.junctions = []
        namespace.new_scoring = None
        namespace.min_clustering_cds_overlap = 0.2
        namespace.min_clustering_cdna_overlap = 0.2
        with tempfile.TemporaryDirectory() as folder:
            out = os.path.join(folder, "configuration.yaml")
            for trial in (1066, 175108):  # (None, 1066, 175108):
                with self.subTest(trial=trial):
                    namespace.mode = ["permissive"]
                    namespace.seed = trial
                    assert namespace.seed is not False
                    namespace.multiprocessing_method = "spawn"
                    with open(out, "w") as out_handle:
                        namespace.out = out_handle
                        sub_configure.create_config(namespace)
                    self.assertGreater(os.stat(out).st_size, 0)
                    conf = configuration.configurator.load_and_validate_config(out)
                    conf = configuration.configurator.check_and_load_scoring(conf)
                    conf = configuration.configurator.check_and_load_scoring(conf)
                    if trial is not None:
                        self.assertEqual(conf.seed, trial)
                    else:
                        self.assertNotEqual(conf.seed, trial)
                        self.assertIsInstance(conf.seed, int)

            with self.subTest(mistake=False):
                with self.assertRaises(OSError):
                    namespace.seed = False
                    namespace.daijin = False
                    namespace.mode = ["permissive"]
                    self.assertFalse(namespace.seed)
                    with open(out, "w") as out_handle:
                        namespace.out = out_handle
                        sub_configure.create_config(namespace)

            with self.subTest(mistake="hello"):
                with self.assertRaises((OSError, InvalidConfiguration)):
                    namespace.seed = "hello"
                    namespace.daijin = False
                    namespace.mode = ["permissive"]
                    with open(out, "w") as out_handle:
                        namespace.out = out_handle
                        sub_configure.create_config(namespace)

            with self.subTest(mistake=b"890"):
                with self.assertRaises((OSError, InvalidConfiguration)):
                    namespace.seed = b"890"
                    namespace.daijin = False
                    namespace.mode = ["permissive"]
                    with open(out, "w") as out_handle:
                        namespace.out = out_handle
                        sub_configure.create_config(namespace)

            with self.subTest(mistake=10.5):
                with self.assertRaises((OSError, InvalidConfiguration)):
                    namespace.seed = 10.5
                    namespace.daijin = False
                    namespace.mode = ["permissive"]
                    with open(out, "w") as out_handle:
                        namespace.out = out_handle
                        sub_configure.create_config(namespace)

    def test_mikado_config_full(self):
        namespace = Namespace(default=False)
        namespace.scoring = None
        namespace.intron_range = None
        namespace.reference = self.fai.filename.decode()
        namespace.external = None
        namespace.mode = ["permissive"]
        namespace.threads = 1
        namespace.blast_targets = []
        namespace.junctions = []
        namespace.new_scoring = None
        namespace.full = True
        namespace.daijin = False
        namespace.seed = 0
        namespace.min_clustering_cds_overlap = 0.2
        namespace.min_clustering_cdna_overlap = 0.2
        dir = tempfile.TemporaryDirectory()
        out = os.path.join(dir.name, "configuration.yaml")
        with open(out, "w") as out_handle:
            namespace.out = out_handle
            sub_configure.create_config(namespace)
        self.assertGreater(os.stat(out).st_size, 0)
        conf = configuration.configurator.load_and_validate_config(out)
        conf = configuration.configurator.check_and_load_scoring(conf)
        conf = configuration.configurator.check_and_load_scoring(conf)
        dir.cleanup()

    def test_mikado_config_daijin(self):
        namespace = Namespace(default=False)
        namespace.scoring = None
        namespace.seed = 0
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
        namespace.min_clustering_cds_overlap = 0.2
        namespace.min_clustering_cdna_overlap = 0.2
        namespace.blast_chunks = 10
        namespace.scheduler = ""
        dir = tempfile.TemporaryDirectory()
        out = os.path.join(dir.name, "configuration.yaml")
        with open(out, "w") as out_handle:
            namespace.out = out_handle
            sub_configure.create_config(namespace)
        self.assertGreater(os.stat(out).st_size, 0)
        conf = configuration.configurator.load_and_validate_config(out)
        conf = configuration.configurator.check_and_load_scoring(conf)
        conf = configuration.configurator.check_and_load_scoring(conf)
        dir.cleanup()

    def test_mikado_config_daijin_set_from_mode(self):
        namespace = Namespace(default=False)
        namespace.scoring = None
        namespace.blast_chunks = 10
        namespace.scheduler = ""
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
        namespace.seed = 0
        namespace.min_clustering_cds_overlap = 0.2
        namespace.min_clustering_cdna_overlap = 0.2
        dir = tempfile.TemporaryDirectory()
        out = os.path.join(dir.name, "configuration.yaml")
        with open(out, "w") as out_handle:
            namespace.out = out_handle
            sub_configure.create_config(namespace)
        self.assertGreater(os.stat(out).st_size, 0)
        conf = configuration.configurator.load_and_validate_config(out)
        conf = configuration.configurator.check_and_load_scoring(conf)
        conf = configuration.configurator.check_and_load_scoring(conf)
        dir.cleanup()

    @mark.slow
    def test_daijin_config(self):

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

        scorers = []
        score__folder = pkg_resources.resource_filename("Mikado.configuration", "scoring_files")
        for root, _, files in os.walk(score__folder):
            for fname in files:
                scorers.append(os.path.join(root, fname))

        for iteration in range(20):
            with self.subTest(iteration=iteration), tempfile.TemporaryDirectory() as folder:
                namespace.out_dir = folder
                namespace.scoring = scorers[np.random.choice(len(scorers))]
                out = os.path.join(folder, "configuration.yaml")
                config = DaijinConfiguration()
                with open(out, "wt") as out_handle:
                    namespace.out = out_handle
                    daijin_configurator.create_daijin_config(namespace, config, level="ERROR")
                self.assertGreater(os.stat(out).st_size, 0)

                with open(out) as out_handle:
                    config = yaml.load(out_handle, Loader=yLoader)

                config = load_and_validate_config(config)
                self.assertIsInstance(config, DaijinConfiguration)


class DaijinTest(unittest.TestCase):
    """This test case will check that daijin can load configurations and snakemake correctly.
    All the tests will be done in 'dry-run' mode."""

    @classmethod
    def setUpClass(cls):
        cls.fai = pysam.FastaFile(pkg_resources.resource_filename("Mikado.tests", "chr5.fas.gz"))

    def test_mikado_dry_run(self):
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
        namespace.dryrun = True
        folder = tempfile.TemporaryDirectory()
        namespace.out_dir = folder.name
        score__folder = pkg_resources.resource_filename("Mikado.configuration", "scoring_files")
        scorers = []
        for root, _, files in os.walk(score__folder):
            for fname in files:
                scorers.append(os.path.join(root, fname))
        namespace.scoring = scorers[np.random.choice(len(scorers))]
        out = os.path.join(folder.name, "configuration.yaml")
        config = DaijinConfiguration()
        with open(out, "wt") as out_handle:
            namespace.out = out_handle
            daijin_configurator.create_daijin_config(namespace, config, level="ERROR")
        # Now the dry run
        namespace.config = out_handle.name
        mikado_pipeline(namespace)

    def test_assemble_dry_run(self):
        # "daijin configure -as stringtie scallop -lal gmap --scheduler local -al \
        # hisat gsnap --sample-sheet samples.txt -o {output} -g {input.genome} -od daijin_test \
        #  --prot-db {input.prots} --scoring plant.yaml"
        namespace = Namespace(default=False)
        folder = tempfile.TemporaryDirectory()
        file_list = os.path.join(folder.name, "samples.txt")
        # ERR588038.R1.fq.gz      ERR588038.R2.fq.gz      ERR588038       fr-firststrand  False
        # SRR5037293.pacbio.fastq.gz              SRR5037293      fr-firststrand  True
        r1 = pkg_resources.resource_filename("Mikado.tests", "ERR588038.R1.fq.gz")
        r2 = pkg_resources.resource_filename("Mikado.tests", "ERR588038.R2.fq.gz")
        pb = pkg_resources.resource_filename("Mikado.tests", "SRR5037293.pacbio.fastq.gz")
        with open(file_list, "wt") as flist:
            print(r1, r2, "ERR588038", "fr-firststrand", False, sep="\t", file=flist)
            print(pb, "", "SRR5037293", "fr-firststrand", True, sep="\t", file=flist)

        namespace.r1 = [r1]
        namespace.r2 = [r2]
        namespace.samples = ["ERR588038"]
        namespace.strandedness = ["fr-firststrand"]
        namespace.asm_methods = ["stringtie", "scallop"]
        namespace.aligners = ["hisat", "gsnap"]
        namespace.modes = ["nosplit"]
        namespace.cluster_config = None
        namespace.scheduler = "local"
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
        namespace.dryrun = True
        namespace.out_dir = folder.name
        out = os.path.join(folder.name, "configuration.yaml")
        config = DaijinConfiguration()
        with open(out, "wt") as out_handle:
            namespace.out = out_handle
            daijin_configurator.create_daijin_config(namespace, config, level="ERROR")
        # Now the dry run
        namespace.config = out_handle.name
        assemble_transcripts_pipeline(namespace)
        namespace.r1, namespace.r2 = [], []
        namespace.sample_sheet = file_list
        assemble_transcripts_pipeline(namespace)


class PickUtilsTest(unittest.TestCase):
    def test_check_regions_parsing(self):
        self.assertEqual(_parse_regions(None), None)
        it_d = _parse_regions("Chr1:1000-1500")
        self.assertIsInstance(it_d["Chr1"], IntervalTree)
        self.assertEqual((it_d["Chr1"].start, it_d["Chr1"].end), (1000, 1500))
        it_d = _parse_regions("Chr1:1000..1500")
        self.assertIsInstance(it_d["Chr1"], IntervalTree)
        self.assertEqual((it_d["Chr1"].start, it_d["Chr1"].end), (1000, 1500))
        with self.assertRaises(ValueError):
            _ = _parse_regions("Chr1:1000.1500")

        with tempfile.NamedTemporaryFile(suffix=".txt", mode="wt") as region_temp:
            print("Chr1:1000-1500", file=region_temp)
            print("Chr1:10000..15000", file=region_temp)
            print("Chr2:30000..35000", file=region_temp)
            print("Chr2:40000-45000", file=region_temp)
            region_temp.flush()
            regions = _parse_regions(region_temp.name)
            self.assertEqual(sorted(regions.keys()), ["Chr1", "Chr2"])
            self.assertEqual((regions["Chr1"].start, regions["Chr1"].end), (1000, 15000))
            self.assertEqual((regions["Chr2"].start, regions["Chr2"].end), (30000, 45000))

        with tempfile.NamedTemporaryFile(suffix=".txt", mode="wt") as region_err:
            print("Chr1:1000-1500", file=region_err)
            print("Chr1:10000..15000", file=region_err)
            print("Chr2:30000..35000", file=region_err)
            print("Chr2:40000-45000", file=region_err)
            print("Chr2:50000.55000", file=region_err)
            region_err.flush()
            with self.assertRaises(ValueError) as exc:
                _ = _parse_regions(region_err.name)
            self.assertEqual(str(exc.exception).rstrip(), "Invalid region line, no. 5: Chr2:50000.55000")

    def test_set_mode(self):
        conf = MikadoConfiguration()
        for invalid in [10, None, "hello", "lenient2", "str1ngent"]:
            with self.assertRaises(InvalidConfiguration):
                _set_pick_mode(conf, "wrong")

        for modebase in ["lenient", "stringent", "permissive", "nosplit", "split"]:
            for mode in [modebase.lower(), modebase.upper(), modebase.capitalize()]:
                with self.subTest(mode=mode):
                    conf = _set_pick_mode(conf, mode)
                    if mode.lower() == "nosplit":
                        self.assertFalse(conf.pick.chimera_split.execute)
                    else:
                        self.assertTrue(conf.pick.chimera_split.execute)
                    if mode.lower() in ["nosplit", "split"]:
                        self.assertFalse(conf.pick.chimera_split.blast_check)
                    else:
                        self.assertTrue(conf.pick.chimera_split.blast_check)
                        self.assertEqual(conf.pick.chimera_split.blast_params.leniency, mode.upper())

    def test_create_log(self):
        # Test for _utils.check_log_settings_and_create_logger
        conf = MikadoConfiguration()
        for invalid in [10, 50.0, dict(), None, b"nonexistent", b"pick2", "nonexistent", "pick2"]:
            with self.subTest(invalid=invalid):
                with self.assertRaises(InvalidConfiguration):
                    check_log_settings_and_create_logger(invalid, None, "INFO", section=None)
                if invalid is not None:
                    with self.assertRaises((AttributeError, TypeError)):
                        check_log_settings_and_create_logger(conf, None, "INFO", section=invalid)
                if invalid is not None and not isinstance(invalid, (bytes, str)):
                    with self.assertRaises(TypeError):
                        check_log_settings_and_create_logger(conf, invalid, "INFO", section=None)
                if invalid is not None:
                    with self.assertRaises((TypeError, AttributeError, InvalidConfiguration)):
                        check_log_settings_and_create_logger(conf, "mikado.log", invalid, section=None)

        # Let's check that we are setting fields correctly now.
        curr_dir = os.getcwd()
        os.chdir(tempfile.gettempdir())
        levels = []
        [levels.extend([level.upper(), level.lower(), level.capitalize()]) for level in
         ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]]

        conf = MikadoConfiguration()
        try:
            for section, handle, level in itertools.product(
                    [None, "prepare", "serialise", "pick"],
                    [None, "stderr", "mikado.log", b"mikado.log", open("mikado.log", "wt"),
                     "folder/mikado.log"],
                    levels):
                with self.subTest(section=section, level=level, handle=handle):
                    conf.pick.files.output_dir = "pick"
                    conf.pick.files.log = "pick.log"
                    conf.prepare.files.output_dir = "prepare"
                    conf.prepare.files.log = "prepare.log"
                    conf.serialise.files.output_dir = "serialise"
                    conf.serialise.files.log = "serialise.log"
                    conf, logger = check_log_settings_and_create_logger(conf, handle, level, section)
                    self.assertEqual(len(logger.handlers), 1)
                    handler = logger.handlers[0]
                    self.assertEqual(logger.level, logging.getLevelName(level.upper()), logger.level)
                    if (section is None and handle is None) or handle == "stderr":
                        self.assertIsInstance(handler, logging.StreamHandler)
                        if section == "pick":
                            self.assertIsNone(conf.pick.files.log)
                        elif section == "serialise":
                            self.assertIsNone(conf.serialise.files.log)
                        elif section == "prepare":
                            self.assertIsNone(conf.prepare.files.log)
                    else:
                        self.assertIsInstance(handler, logging.FileHandler)
                        fname = handler.stream.name
                        if handle is None:
                            if section == "prepare":
                                self.assertEqual(fname, os.path.join(os.getcwd(), "prepare", "prepare.log"))
                            elif section == "serialise":
                                self.assertEqual(fname, os.path.join(os.getcwd(), "serialise", "serialise.log"))
                            elif section == "pick":
                                self.assertEqual(fname, os.path.join(os.getcwd(), "pick", "pick.log"))
                        else:
                            if handle == b"mikado.log":
                                check = "mikado.log"
                            elif isinstance(handle, (io.BufferedWriter, io.TextIOWrapper)):
                                check = "mikado.log"
                            else:
                                check = handle
                            if section == "prepare":
                                if not os.path.dirname(check):
                                    self.assertEqual(fname, os.path.join(os.getcwd(), conf.prepare.files.output_dir,
                                                                         check))
                                else:
                                    self.assertEqual(fname, os.path.join(os.getcwd(), handle))
                                other = conf.prepare.files.log
                            elif section == "serialise":
                                if not os.path.dirname(check):
                                    self.assertEqual(fname, os.path.join(os.getcwd(),
                                                                         conf.serialise.files.output_dir,
                                                                         check))
                                else:
                                    self.assertEqual(fname, os.path.join(os.getcwd(), handle))
                                other = conf.serialise.files.log
                            elif section == "pick":
                                if not os.path.dirname(check):
                                    self.assertEqual(fname, os.path.join(os.getcwd(),
                                                                         conf.pick.files.output_dir,
                                                                         check))
                                else:
                                    self.assertEqual(fname, os.path.join(os.getcwd(), handle))
                                other = conf.pick.files.log
                            else:
                                other = None
                            if other is not None:
                                self.assertEqual(conf.log_settings.log, other)

                        self.assertTrue(os.path.exists(fname))
                        if os.path.dirname(fname) != tempfile.gettempdir():
                            shutil.rmtree(os.path.dirname(fname))
                        else:
                            os.remove(fname)
        except AssertionError:
            raise
        finally:
            os.chdir(curr_dir)


@mark.slow
class PickTest(unittest.TestCase):

    """This unit test will check that pick functions correctly."""

    def setUp(self):
        
        self.configuration = configurator.load_and_validate_config(None)
        self.configuration.reference.genome = self.fai.filename.decode()

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
    def test_single_and_multi_proc(self):

        self.configuration.threads = 1
        self.configuration.db_settings.db = pkg_resources.resource_filename("Mikado.tests", "mikado.db")

        self.configuration.pick.files.input = pkg_resources.resource_filename("Mikado.tests",
                                                                                   "mikado_prepared.gtf")
        with tempfile.TemporaryDirectory() as folder:
            self.configuration.pick.files.output_dir = folder
            self.configuration.pick.files.loci_out = "mikado.monoproc.loci.gff3"
            self.configuration.pick.files.subloci_out = "mikado.monoproc.subloci.gff3"
            self.configuration.pick.files.monoloci_out = "mikado.monoproc.monoloci.gff3"
            self.configuration.pick.files.log = "mikado.monoproc.log"
            self.configuration.pick.alternative_splicing.pad = False
            self.configuration.log_settings.log_level = "WARNING"

            pick_caller = picker.Picker(configuration=self.configuration)
            with self.assertRaises(SystemExit), self.assertLogs("main_logger", "INFO") as cm:
                pick_caller()
            self.assertTrue(os.path.exists(os.path.join(folder, "mikado.monoproc.loci.gff3")))
            single_scores = pd.read_csv(os.path.join(folder, "mikado.monoproc.loci.scores.tsv"),
                                        delimiter="\t").sort_values("alias")
            single_metrics = pd.read_csv(os.path.join(folder, "mikado.monoproc.loci.metrics.tsv"),
                                         delimiter="\t").sort_values("alias")
            with parser_factory(os.path.join(folder, "mikado.monoproc.loci.gff3")) as inp_gff:
                lines = [_ for _ in inp_gff if not _.header is True]
                self.assertGreater(len(lines), 0)
                self.assertGreater(len([_ for _ in lines if _.is_transcript is True]), 0)
                self.assertGreater(len([_ for _ in lines if _.feature == "mRNA"]), 0,
                                   [_ for _ in cm.output if "WARNING" in _])
                self.assertGreater(len([_ for _ in lines if _.feature == "CDS"]), 0)
            self.assertTrue(np.array_equal(single_scores["score"].values, single_metrics["score"].values))
            self.assertFalse(single_metrics["source_score"].isna().any())
            self.assertFalse(single_scores["source_score"].isna().any())

        self.configuration.threads = 2
        self.configuration.pick.files.input = pkg_resources.resource_filename("Mikado.tests",
                                                                              "mikado_prepared.gtf")
        with tempfile.TemporaryDirectory() as folder:
            self.configuration.pick.files.output_dir = folder
            self.configuration.pick.files.loci_out = "mikado.multiproc.loci.gff3"
            self.configuration.pick.files.subloci_out = "mikado.multiproc.subloci.gff3"
            self.configuration.pick.files.monoloci_out = "mikado.multiproc.monoloci.gff3"
            self.configuration.pick.files.log = "mikado.multiproc.log"
            self.configuration.db_settings.db = pkg_resources.resource_filename("Mikado.tests", "mikado.db")
            self.configuration.pick.alternative_splicing.pad = False
            self.configuration.log_settings.log_level = "INFO"
            pick_caller = picker.Picker(configuration=self.configuration)
            with self.assertRaises(SystemExit), self.assertLogs("main_logger", "INFO"):
                pick_caller()
            self.assertTrue(os.path.exists(os.path.join(folder, "mikado.multiproc.loci.gff3")))
            with parser_factory(os.path.join(folder, "mikado.multiproc.loci.gff3")) as inp_gff:
                lines = [_ for _ in inp_gff if not _.header is True]
                self.assertGreater(len(lines), 0)
                self.assertGreater(len([_ for _ in lines if _.is_transcript is True]), 0)
                self.assertGreater(len([_ for _ in lines if _.feature == "mRNA"]), 0)
                self.assertGreater(len([_ for _ in lines if _.feature == "CDS"]), 0)
            log = [line for line in open(os.path.join(folder, self.configuration.pick.files.log))]
            self.assertTrue(os.path.exists(os.path.join(folder, "mikado.multiproc.loci.gff3")))
            self.assertGreater(os.stat(os.path.join(folder, "mikado.multiproc.loci.gff3")).st_size, 0,
                               log)
            self.assertGreater(os.stat(os.path.join(folder, "mikado.multiproc.loci.scores.tsv")).st_size, 0,
                               [(name, os.stat(os.path.join(folder, name)).st_size) for name in os.listdir(folder)])
            multi_scores = pd.read_csv(os.path.join(folder, "mikado.multiproc.loci.scores.tsv"),
                                       delimiter="\t").sort_values("alias")
            multi_metrics = pd.read_csv(os.path.join(folder, "mikado.multiproc.loci.metrics.tsv"),
                                        delimiter="\t").sort_values("alias")
            self.assertTrue(np.array_equal(multi_scores["score"].values, multi_metrics["score"].values))
            self.assertFalse(multi_metrics["source_score"].isna().any())
            self.assertFalse(multi_scores["source_score"].isna().any())
        self.assertTrue(np.array_equal(multi_scores, single_scores))
        self.assertTrue(np.array_equal(multi_metrics, single_metrics))

    @mark.slow
    def test_subprocess(self):
                
        self.configuration.pick.files.input = pkg_resources.resource_filename("Mikado.tests",
                                                                                   "mikado_prepared.gtf")
        with tempfile.TemporaryDirectory() as folder:
            self.configuration.pick.files.output_dir = folder
            self.configuration.pick.files.loci_out = "mikado.subproc.loci.gff3"
            self.configuration.pick.files.subloci_out = "mikado.subproc.subloci.gff3"
            self.configuration.pick.files.monoloci_out = "mikado.subproc.monoloci.gff3"
            self.configuration.pick.alternative_splicing.pad = False
            self.configuration.pick.files.log = "mikado.subproc.log"
            self.configuration.db_settings.db = str(pkg_resources.resource_filename("Mikado.tests", "mikado.db"))
            self.configuration.log_settings.log_level = "WARNING"
    
            for num in (1, 2):
                with self.subTest(num=num):
                    self.configuration.pick.run_options.single_thread = (num == 1)
                    json_file = os.path.join(folder, "mikado.yaml")

                    with open(json_file, "wt") as json_handle:
                        print_config(self.configuration, json_handle, output_format="yaml")
    
                    sys.argv = ["mikado", "pick", "--json-conf", json_file, "--seed", "1078"]
                    with self.assertRaises(SystemExit):
                        pkg_resources.load_entry_point("Mikado", "console_scripts", "mikado")()
    
                    self.assertTrue(os.path.exists(os.path.join(folder, "mikado.subproc.loci.gff3")))
                    with parser_factory(os.path.join(folder, "mikado.subproc.loci.gff3")) as inp_gff:
                        lines = [_ for _ in inp_gff if not _.header is True]
                        self.assertGreater(len(lines), 0)
                        self.assertGreater(len([_ for _ in lines if _.is_transcript is True]), 0)
                        self.assertGreater(len([_ for _ in lines if _.feature == "mRNA"]), 0)
                        self.assertGreater(len([_ for _ in lines if _.feature == "CDS"]), 0)
                    [os.remove(_) for _ in glob.glob(os.path.join(folder, "mikado.subproc.") + "*")]

    @mark.slow
    @unittest.skipUnless(os.path.exists("/dev/shm") and os.access("/dev/shm", os.W_OK),
                         "/dev/shm not present or not writeable")
    def test_subprocess_shm(self):
        self.configuration.pick.files.input = pkg_resources.resource_filename("Mikado.tests",
                                                                                   "mikado_prepared.gtf")
        self.configuration.pick.files.loci_out = "mikado.subproc.loci.gff3"
        self.configuration.pick.files.subloci_out = "mikado.subproc.subloci.gff3"
        self.configuration.pick.files.monoloci_out = "mikado.subproc.monoloci.gff3"
        self.configuration.pick.alternative_splicing.pad = False
        self.configuration.pick.files.log = "mikado.subproc.log"
        self.configuration.db_settings.db = str(pkg_resources.resource_filename("Mikado.tests", "mikado.db"))
        self.configuration.log_settings.log_level = "WARNING"

        for num, shm in itertools.product((1, 2), (True,)):
            with self.subTest(num=num, shm=shm), tempfile.TemporaryDirectory() as folder:
                self.configuration.pick.files.output_dir = folder
                self.configuration.pick.run_options.single_thread = (num == 1)
                json_file = os.path.join(folder, "mikado.yaml")

                with open(json_file, "wt") as json_handle:
                    print_config(self.configuration, json_handle, output_format="yaml")

                log = "pick.log"
                if os.path.exists(os.path.join(folder, log)):
                    os.remove(os.path.join(folder, log))
                sys.argv = ["mikado", "pick", "--configuration", json_file, "--seed", "1078", "--log", log]
                if shm is True:
                    sys.argv.append("--shm")
                with self.assertRaises(SystemExit):
                    pkg_resources.load_entry_point("Mikado", "console_scripts", "mikado")()

                self.assertTrue(os.path.exists(os.path.join(folder, "mikado.subproc.loci.gff3")))
                with parser_factory(os.path.join(folder, "mikado.subproc.loci.gff3")) as inp_gff:
                    lines = [_ for _ in inp_gff if not _.header is True]
                    self.assertGreater(len(lines), 0)
                    self.assertGreater(len([_ for _ in lines if _.is_transcript is True]), 0)
                    self.assertGreater(len([_ for _ in lines if _.feature == "mRNA"]), 0)
                    self.assertGreater(len([_ for _ in lines if _.feature == "CDS"]), 0)
                with open(os.path.join(folder, log)) as hlog:
                    log_lines = [_.rstrip() for _ in hlog]
                if shm is True:
                    self.assertTrue(any("Copying Mikado database into a SHM db" in _ for _ in log_lines))

    @mark.slow
    def test_different_scoring(self):

        with tempfile.TemporaryDirectory(prefix="test_different_scoring") as folder:
            self.configuration.pick.files.output_dir = os.path.abspath(folder)
            self.configuration.pick.files.input = pkg_resources.resource_filename("Mikado.tests",
                                                                                       "mikado_prepared.gtf")

            self.configuration.pick.files.loci_out = "mikado.test_diff.loci.gff3"
            self.configuration.pick.files.subloci_out = "mikado.test_diff.subloci.gff3"
            self.configuration.pick.files.monoloci_out = "mikado.test_diff.monoloci.gff3"
            self.configuration.pick.files.log = "mikado.test_diff.log"
            self.configuration.pick.alternative_splicing.pad = False
            self.configuration.log_settings.log_level = "DEBUG"

            self.assertEqual(os.path.basename(self.configuration.pick.scoring_file), "plant.yaml")
            shutil.copy(pkg_resources.resource_filename("Mikado.tests", "mikado.db"),
                        os.path.join(self.configuration.pick.files.output_dir, "mikado.db"))
            self.configuration.db_settings.db = os.path.join(self.configuration.pick.files.output_dir, "mikado.db")
            json_file = os.path.join(self.configuration.pick.files.output_dir, "mikado.yaml")
            with open(json_file, "wt") as json_handle:
                print_config(self.configuration, json_handle, output_format="yaml")
            sys.argv = ["mikado", "pick", "--json-conf", json_file, "--single", "--seed", "1078"]
            with self.assertRaises(SystemExit):
                pkg_resources.load_entry_point("Mikado", "console_scripts", "mikado")()
            with open(os.path.join(self.configuration.pick.files.output_dir, "mikado.test_diff.loci.scores.tsv")) as tsv:
                reader = csv.DictReader(tsv, delimiter="\t")
                score_names = [_ for _ in self.configuration.scoring.scoring]
                score_header = [_ for _ in reader.fieldnames if _ not in
                                ("tid", "alias", "parent", "score", "source_score")]
                self.assertEqual(sorted(score_names), sorted(score_header))

    @mark.slow
    def test_different_scoring_2(self):

        self.configuration.pick.files.input = pkg_resources.resource_filename("Mikado.tests", "mikado_prepared.gtf")

        self.configuration.pick.files.loci_out = "mikado.test_diff.loci.gff3"
        self.configuration.pick.files.subloci_out = "mikado.test_diff.subloci.gff3"
        self.configuration.pick.files.monoloci_out = "mikado.test_diff.monoloci.gff3"
        self.configuration.pick.files.log = "mikado.test_diff.log"
        self.configuration.pick.alternative_splicing.pad = False
        self.configuration.log_settings.log_level = "DEBUG"

        self.assertEqual(os.path.basename(self.configuration.pick.scoring_file), "plant.yaml")
        with tempfile.TemporaryDirectory(prefix="test_different_scoring_2") as outdir:
            shutil.copy(pkg_resources.resource_filename("Mikado.tests", "mikado.db"),
                        os.path.join(outdir, "mikado.db"))
            self.configuration.db_settings.db = os.path.join(outdir, "mikado.db")
            self.configuration.pick.files.output_dir = os.path.join(outdir)
            json_file = os.path.join(outdir, "mikado.yaml")
            with open(json_file, "wt") as json_handle:
                print_config(self.configuration, json_handle, output_format="yaml")
            self.configuration.pick.files.output_dir = os.path.join(outdir)
            scoring_file = pkg_resources.resource_filename("Mikado.tests", "scoring_only_cds.yaml")
            sys.argv = ["mikado", "pick", "--json-conf", json_file, "--single",
                        "--scoring-file", scoring_file, "--seed", "1078"]

            with self.assertRaises(SystemExit):
                pkg_resources.load_entry_point("Mikado", "console_scripts", "mikado")()

            import csv
            with open(os.path.join(self.configuration.pick.files.output_dir, "mikado.test_diff.loci.scores.tsv")) as tsv:
                reader = csv.DictReader(tsv, delimiter="\t")
                score_header = [_ for _ in reader.fieldnames if _ not in
                                ("tid", "alias", "parent", "score", "source_score")]
                self.assertEqual(score_header, ["selected_cds_length"])

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

        self.configuration.pick.files.input = temp_gtf.name
        self.configuration.db_settings.db = os.path.join(dir.name, "mikado.db")
        self.configuration.pick.files.output_dir = dir.name
        self.configuration.log_settings.log_level = "WARNING"
        self.configuration.pick.alternative_splicing.pad = False  # Necessary!
        # del self.configuration["scoring"]
        # del self.configuration["requirements"]
        # del self.configuration["as_requirements"]
        # del self.configuration["not_fragmentary"]
        scoring = dict()

        scoring["requirements"] = dict()
        scoring["requirements"]["expression"] = ["exon_num"]
        scoring["requirements"]["parameters"] = dict()
        scoring["requirements"]["parameters"]["exon_num"] = dict()
        scoring["requirements"]["parameters"]["exon_num"]["name"] = "exon_num"
        scoring["requirements"]["parameters"]["exon_num"]["operator"] = "gt"
        scoring["requirements"]["parameters"]["exon_num"]["value"] = 1
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
        self.configuration.pick.scoring_file = scoring_file.name

        for purging in (False, True):
            with self.subTest(purging=purging):
                self.configuration.pick.files.loci_out = "mikado.purging_{}.loci.gff3".format(purging)
                self.configuration.pick.files.log = os.path.join(dir.name, "mikado.purging_{}.log".format(purging))
                self.configuration.pick.clustering.purge = purging
                self.configuration.pick.scoring_file = scoring_file.name
                self.configuration = configurator.check_and_load_scoring(self.configuration)
                self.assertEqual(len(self.configuration.scoring.scoring.keys()), 1,
                                 self.configuration.scoring.scoring.keys())

                pick_caller = picker.Picker(configuration=self.configuration)
                with self.assertRaises(SystemExit), self.assertLogs("main_logger", "INFO"):
                    pick_caller()

                with parser_factory(os.path.join(dir.name,
                                                 self.configuration.pick.files.loci_out)) as gff:
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
        self.configuration.pick.scoring_file = scoring_file.name

        for purging in (False, True):
            with self.subTest(purging=purging):
                self.configuration.pick.files.loci_out = "mikado.purging_{}.loci.gff3".format(purging)
                self.configuration.pick.files.subloci_out = "mikado.purging_{}.subloci.gff3".format(purging)
                self.configuration.pick.files.log = os.path.join(
                    folder.name,
                    "mikado.purging_{}.log".format(purging))
                self.configuration.pick.clustering.purge = purging
                self.configuration.pick.scoring_file = scoring_file.name
                self.configuration = configurator.check_and_load_scoring(self.configuration)
                self.assertEqual(len(self.configuration.scoring.scoring.keys()), 2,
                                 self.configuration.scoring.scoring.keys())

                pick_caller = picker.Picker(configuration=self.configuration)
                with self.assertRaises(SystemExit), self.assertLogs("main_logger", "INFO"):
                    pick_caller()

                with parser_factory(os.path.join(folder.name,
                                                 self.configuration.pick.files.loci_out)) as gff:
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
        self.configuration.pick.scoring_file = scoring_file.name

        temp_gtf = tempfile.NamedTemporaryFile(mode="wt", suffix=".gtf", delete=True, dir=dir.name)

        gtf = "\n".join([_ for _ in gtf.split("\n") if "foo1.1" not in _])

        temp_gtf.write(gtf)
        temp_gtf.flush()
        self.configuration.pick.files.input = temp_gtf.name

        for purging in (False, True):
            with self.subTest(purging=purging):
                self.configuration.pick.files.loci_out = "mikado.purging_{}.loci.gff3".format(purging)
                self.configuration.pick.files.subloci_out = "mikado.purging_{}.subloci.gff3".format(purging)
                self.configuration.pick.files.log = "mikado.purging_{}.log".format(purging)
                self.configuration.pick.clustering.purge = purging
                self.configuration.pick.scoring_file = scoring_file.name
                self.configuration = configurator.check_and_load_scoring(self.configuration)
                self.assertEqual(len(self.configuration.scoring.scoring.keys()), 2,
                                 self.configuration.scoring.scoring.keys())

                pick_caller = picker.Picker(configuration=self.configuration)
                with self.assertRaises(SystemExit), self.assertLogs("main_logger", "INFO"):
                    pick_caller()

                with parser_factory(os.path.join(dir.name,
                                                 self.configuration.pick.files.loci_out)) as gff:
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
        self.configuration = configurator.load_and_validate_config(None)
        self.configuration.reference.genome = self.fai.filename.decode()

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

        with tempfile.TemporaryDirectory(prefix="test_subprocess_single") as folder:
            json_file = os.path.join(folder, "mikado.yaml")
            db = os.path.join(folder, "mikado.db")
            log = os.path.join(folder, "serialise.log")
            uni_out = os.path.join(folder, "uniprot_sprot_plants.fasta")
            with gzip.open(uniprot, "rb") as uni, open(uni_out, "wb") as uni_out_handle:
                uni_out_handle.write(uni.read())

            with open(json_file, "wt") as json_handle:
                print_config(self.configuration, json_handle, output_format="yaml")
            # Set up the command arguments
            for procs in (1, ):
                with self.subTest(proc=procs):
                    sys.argv = [str(_) for _ in ["mikado", "serialise", "--configuration", json_file,
                                "--transcripts", transcripts, "--blast_targets", uni_out,
                                "--orfs", orfs, "--junctions", junctions, "--xml", xml, "-od", folder,
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

    def test_subprocess_multi(self):

        xml = pkg_resources.resource_filename("Mikado.tests", "chunk-001-proteins.xml.gz")
        transcripts = pkg_resources.resource_filename("Mikado.tests", "mikado_prepared.fasta")
        junctions = pkg_resources.resource_filename("Mikado.tests", "junctions.bed")
        orfs = pkg_resources.resource_filename("Mikado.tests", "transcripts.fasta.prodigal.gff3")
        uniprot = pkg_resources.resource_filename("Mikado.tests", "uniprot_sprot_plants.fasta.gz")
        mobjects = 300  # Let's test properly the serialisation for BLAST

        # Set up the command arguments
        with tempfile.TemporaryDirectory(suffix="test_subprocess_multi_{}".format(1)) as folder_one, \
                tempfile.TemporaryDirectory(suffix="test_subprocess_multi_{}".format(1)) as folder_three:
            for procs, folder in [(1, folder_one), (3, folder_three)]:
                with self.subTest(proc=procs):
                    json_file = os.path.join(folder, "mikado.yaml")
                    db = os.path.join(folder, "mikado.db")
                    log = os.path.join(folder, "serialise.log")
                    uni_out = os.path.join(folder, "uniprot_sprot_plants.fasta")
                    with gzip.open(uniprot, "rb") as uni, open(uni_out, "wb") as uni_out_handle:
                        uni_out_handle.write(uni.read())
                    with open(json_file, "wt") as json_handle:
                        print_config(self.configuration, json_handle, output_format="yaml")
                    sys.argv = [str(_) for _ in ["mikado", "serialise", "--json-conf", json_file,
                                "--transcripts", transcripts, "--blast_targets", uni_out,
                                "--orfs", orfs, "--junctions", junctions, "--xml", xml, "-od", folder,
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

    @mark.slow
    def test_xml_vs_tsv(self):
        xml = pkg_resources.resource_filename("Mikado.tests", os.path.join("blast_data", "diamond.0.9.30.xml.gz"))
        tsv = pkg_resources.resource_filename("Mikado.tests", os.path.join("blast_data", "diamond.0.9.30.tsv.gz"))
        queries = pkg_resources.resource_filename("Mikado.tests", os.path.join("blast_data", "transcripts.fasta"))
        prots = pkg_resources.resource_filename("Mikado.tests", "uniprot_sprot_plants.fasta.gz")
        logs = collections.defaultdict(dict)
        dbs = collections.defaultdict(dict)

        with tempfile.TemporaryDirectory(prefix="test_xml_vs_tsv") as test_xml_vs_tsv_folder:
            assert "test_xml_vs_tsv" in test_xml_vs_tsv_folder
            args = Namespace(default=None)
            args.seed = 1078
            args.max_objects = 1000
            args.output_dir = test_xml_vs_tsv_folder
            args.transcripts = queries
            args.blast_targets = [prots]
            for name, blast in zip(["xml", "tsv"], [xml, tsv]):
                for proc in (1, 3):
                    with self.subTest(name=name, blast=blast, proc=proc):
                        args.db = "{}_{}.db".format(name, proc)
                        args.log = "{}_{}.log".format(name, proc)
                        args.xml = blast
                        args.procs = proc
                        args.start_adjustment = True
                        serialise(args)
                        dbs[name][proc] = os.path.join(test_xml_vs_tsv_folder, args.db)
                        logged = [_.rstrip() for _ in open(os.path.join(test_xml_vs_tsv_folder, args.log))]
                        logs[name][proc] = logged

            def prep_dbs(name):
                import sqlalchemy.exc
                try:
                    _ = pd.read_sql(sql="hsp", con=name)
                except sqlalchemy.exc.OperationalError:
                    try:
                        with sqlite3.connect(name) as conn:
                            tables = conn.execute("select sql from sqlite_master where type = 'table'").fetchall()
                        raise sqlite3.OperationalError(tables)
                    except sqlite3.OperationalError:
                        raise sqlite3.OperationalError(logs)

                hsp, hit, query, target = [pd.read_sql(table, name) for table in ["hsp", "hit", "query", "target"]]
                hit = hit.join(target.set_index("target_id"), on=["target_id"], how="inner").join(
                    query.set_index("query_id"), on=["query_id"], how="inner")
                hsp = hsp.join(target.set_index("target_id"), on=["target_id"], how="inner").join(
                    query.set_index("query_id"), on=["query_id"], how="inner")
                hsp.set_index(["query_name", "target_name", "counter"], inplace=True)
                hit.set_index(["query_name", "target_name"], inplace=True)
                return hit, hsp

            try:
                xml_hit, xml_hsp = prep_dbs("sqlite:///" + dbs["xml"][1])
            except KeyError:
                raise KeyError(dbs)
            xml_hit_multi, xml_hsp_multi = prep_dbs("sqlite:///" + dbs["xml"][3])
            tsv_hit, tsv_hsp = prep_dbs("sqlite:///" + dbs["tsv"][1])
            tsv_hit_multi, tsv_hsp_multi = prep_dbs("sqlite:///" + dbs["tsv"][3])

            hit = pd.merge(xml_hit, tsv_hit, left_index=True, right_index=True, suffixes=("_xml", "_tsv"))
            hit_multi = pd.merge(xml_hit_multi, tsv_hit_multi, left_index=True, right_index=True, suffixes=["_xml_m",
                                                                                                            "_tsv_m"])
            hit = pd.merge(hit, hit_multi, left_index=True, right_index=True)
            hsp = pd.merge(xml_hsp, tsv_hsp, left_index=True, right_index=True, suffixes=("_xml", "_tsv"))
            hsp_multi = pd.merge(xml_hsp_multi, tsv_hsp_multi, left_index=True, right_index=True, suffixes=("_xml_m",
                                                                                                            "_tsv_m"))
            hsp = pd.merge(hsp, hsp_multi, left_index=True, right_index=True)
            self.assertTrue(xml_hit_multi.shape[0] == xml_hit.shape[0] == tsv_hit.shape[0] > 0)
            self.assertTrue(tsv_hit_multi.shape[0] == tsv_hit.shape[0] > 0,
                            (tsv_hit_multi.shape[0], tsv_hit.shape[0], logs["tsv"][3]))

            self.assertTrue(hit.shape[0] == xml_hit.shape[0] == tsv_hit.shape[0] > 0)
            self.assertTrue(hsp.shape[0] == xml_hsp.shape[0] == tsv_hsp.shape[0] > 0)
            self.assertTrue(hsp.shape[0] == xml_hsp_multi.shape[0] == tsv_hsp_multi.shape[0] > 0)
            # Get the columns
            hitcols, hspcols = collections.defaultdict(list), collections.defaultdict(list)
            pat = re.compile(r"_(tsv|xml)($|_m$)")
            for d, df in zip([hitcols, hspcols], [hit, hsp]):
                for col in df.columns:
                    name = re.sub(pat, '', col)
                    d[name].append(col)
                failed = []
                for col in d:
                    if col in ("query_id", "target_id"):
                        continue
                    if len(d[col]) == 1:
                        raise ValueError(col)
                    with self.subTest(col=col):
                        catch = df[d[col]].apply(lambda row: row[0] == row[1] or
                                                             np.isclose(row[0], row[1], atol=.01, rtol=.01), axis=1)
                    if not (catch).all():
                        failed_rows = df.loc[~catch, d[col]]
                        failed.append((col, failed_rows.head()))
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
        with tempfile.TemporaryDirectory(prefix="has_to_fail") as folder_one, \
                tempfile.TemporaryDirectory(prefix="has_to_fail") as folder_two:
            for procs, folder in [(3, folder_one), (1, folder_two)]:
                with self.subTest(procs=procs):
                    json_file = os.path.join(folder, "mikado.yaml")
                    db = os.path.join(folder, "mikado.db")
                    log = "failed_serialise.log"
                    uni_out = os.path.join(folder, "uniprot_sprot_plants.fasta")
                    self.configuration.serialise.files.log = os.path.basename(log)
                    self.configuration.multiprocessing_method = "fork"
                    with gzip.open(uniprot, "rb") as uni, open(uni_out, "wb") as uni_out_handle:
                        uni_out_handle.write(uni.read())
                    with open(json_file, "wt") as json_handle:
                        print_config(self.configuration, json_handle, output_format="yaml")
                    with self.subTest(proc=procs):
                        sys.argv = [str(_) for _ in ["mikado", "serialise", "--json-conf", json_file,
                                                     "--transcripts", transcripts, "--blast_targets", uni_out,
                                                     "--log", log,
                                                     "-od", folder,
                                                     "--orfs", tmp_orf.name, "--junctions", junctions, "--xml", xml,
                                                     "-p", procs, "-mo", mobjects, db,
                                                     "--seed", "1078"]]
                        log = os.path.join(folder, log)
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

    def test_serialise_external(self):

        base = pkg_resources.resource_filename("Mikado.tests", "test_external_aphid")
        external_conf = os.path.join(base, "mikado.configuration.testds.yaml")
        external_scores = os.path.join(base, "annotation_run1.metrics.testds.txt")
        fasta = os.path.join(base, "mikado_prepared.testds.fasta")

        with tempfile.TemporaryDirectory(prefix="test_serialise_external_1") as folder_one,\
                tempfile.TemporaryDirectory(prefix="test_serialise_external_3") as folder_two:
            for procs, folder in [(1, folder_one), (3, folder_two)]:
                with self.subTest(procs=procs):
                    log = "serialise.log"
                    sys.argv = [str(_) for _ in ["mikado", "serialise", "--configuration", external_conf,
                                                 "--transcripts", fasta, "-od", folder,
                                                 "-l", log,
                                                 "--external-scores", external_scores,
                                                 "--seed", 10, "--procs", procs, "mikado.db"]]
                    log = os.path.join(folder, log)
                    pkg_resources.load_entry_point("Mikado", "console_scripts", "mikado")()
                    conn = sqlite3.connect(os.path.join(folder, "mikado.db"))
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
                namespace.gff = parser_factory(filename)
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
                namespace.gff = parser_factory(filename)
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

        with tempfile.NamedTemporaryFile("wt", suffix=".txt") as id_temp_file:
            to_write = [ids[_] for _ in np.random.choice(len(ids), 10, replace=False)]
            [print(*idline, sep="\t", file=id_temp_file) for idline in to_write]
            id_temp_file.flush()

            for fname in files:
                form = os.path.splitext(fname)[1]
                with self.subTest(fname=fname), tempfile.NamedTemporaryFile("wt", suffix=form) as outfile:
                    sys.argv = ["mikado", "util", "grep", id_temp_file.name, fname, outfile.name]
                    pkg_resources.load_entry_point("Mikado", "console_scripts", "mikado")()
                    self.assertTrue(os.path.exists(outfile.name))
                    found = set()
                    with parser_factory(outfile.name, input_format=form[1:]) as stream:
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

        with tempfile.NamedTemporaryFile("wt", suffix=".txt") as id_temp_file:
            to_write = [ids[_] for _ in np.random.choice(len(ids), 10, replace=False)]
            others = [_ for _ in ids if _ not in to_write]
            [print(*idline, sep="\t", file=id_temp_file) for idline in to_write]
            id_temp_file.flush()

            for fname in files:
                form = os.path.splitext(fname)[1]
                with self.subTest(fname=fname), tempfile.NamedTemporaryFile("wt", suffix=form) as outfile:
                    sys.argv = ["mikado", "util", "grep", "-v", id_temp_file.name, fname, outfile.name]
                    pkg_resources.load_entry_point("Mikado", "console_scripts", "mikado")()
                    self.assertTrue(os.path.exists(outfile.name))
                    found = set()
                    with parser_factory(outfile.name, input_format=form[1:]) as stream:
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
            form = os.path.splitext(fname)[1]
            with self.subTest(flag=flag), tempfile.NamedTemporaryFile("wt", suffix=form) as outfile:
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

                with parser_factory(outfile.name, input_format=form[1:]) as stream:
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
