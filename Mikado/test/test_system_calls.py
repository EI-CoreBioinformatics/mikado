import unittest
from Mikado.preparation import prepare
from Mikado.configuration import configurator, daijin_configurator
import pkg_resources
import tempfile
from Mikado.loci.transcript import Namespace
from Mikado.utilities.log_utils import create_null_logger
from Mikado.scales.compare import compare, load_index
import logging
import gzip
import pyfaidx
import os
import itertools
import csv
from Mikado.parsers import to_gff
from Mikado.subprograms.util.stats import Calculator
from Mikado.picking import picker
import Mikado.subprograms.configure
import Mikado.daijin
import yaml
import random
import sys


class PrepareChek(unittest.TestCase):

    __genomefile__ = None

    @classmethod
    def setUpClass(cls):
        cls.__genomefile__ = tempfile.NamedTemporaryFile(mode="wb", delete=False, suffix=".fa", prefix="prepare")

        with pkg_resources.resource_stream("Mikado.test", "chr5.fas.gz") as _:
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
                self.conf["prepare"]["files"]["gff"] = [pkg_resources.resource_filename("Mikado.test",
                                                                                        test_file)]

                prepare.prepare(args, self.logger)

                # Now that the program has run, let's check the output
                fa = pyfaidx.Fasta(os.path.join(self.conf["prepare"]["files"]["output_dir"],
                                                "mikado_prepared.fasta"))
                res = dict((_, len(fa[_])) for _ in fa.keys())
                fa.close()
                self.assertEqual(res, self.trinity_res)

    def test_prepare_trinity_and_cufflinks(self):

        self.conf["prepare"]["files"]["labels"] = ["cl", "tr"]
        self.conf["prepare"]["files"]["gff"].append(pkg_resources.resource_filename("Mikado.test",
                                                                                    "cufflinks.gtf"))
        self.conf["prepare"]["files"]["gff"].append("")
        self.conf["prepare"]["files"]["output_dir"] = tempfile.gettempdir()
        self.conf["prepare"]["files"]["out_fasta"] = "mikado_prepared.fasta"
        self.conf["prepare"]["files"]["out"] = "mikado_prepared.gtf"
        args = Namespace()
        args.json_conf = self.conf

        for test_file in ("trinity.gff3",
                          "trinity.match_matchpart.gff3",
                          "trinity.cDNA_match.gff3",
                          "trinity.gtf"):
            with self.subTest(test_file=test_file):
                self.conf["prepare"]["files"]["gff"][1] = pkg_resources.resource_filename("Mikado.test",
                                                                                          test_file)
                self.conf["prepare"]["files"]["out_fasta"] = "mikado_prepared.fasta"
                self.conf["prepare"]["files"]["out"] = "mikado_prepared.gtf"

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


class CompareCheck(unittest.TestCase):

    """Test to check that compare interacts correctly with match, match_part, cDNA_match"""

    def tearDown(self):
        logging.shutdown()

    def test_index(self):

        # Create the list of files
        files = ["trinity.gtf",
                 "trinity.gff3",
                 "trinity.cDNA_match.gff3",
                 "trinity.match_matchpart.gff3"]
        # files = [pkg_resources.resource_filename("Mikado.test", filename) for filename in files]

        namespace = Namespace(default=False)
        namespace.distance = 2000
        namespace.index = True
        namespace.prediction = None
        namespace.log = os.path.join(tempfile.gettempdir(), "index.log")
        logger = create_null_logger("null")

        for ref in files:
            with self.subTest(ref=ref):
                temp_ref = os.path.join(tempfile.gettempdir(), ref)
                with pkg_resources.resource_stream("Mikado.test", ref) as ref_handle,\
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
                os.remove("{}.midx".format(namespace.reference.name))

    def test_compare_trinity(self):

        # Create the list of files
        files = ["trinity.gtf",
                 "trinity.gff3",
                 "trinity.cDNA_match.gff3",
                 "trinity.match_matchpart.gff3"]
        files = [pkg_resources.resource_filename("Mikado.test", filename) for filename in files]

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

        for permutation in itertools.permutations(range(len(files)), 2):
            for suff in ["log", "refmap", "tmap", "stats"]:
                fname = os.path.join(tempfile.gettempdir(),
                                   "compare_{}_{}.{}".format(permutation[0], permutation[1], suff))
                if os.path.exists(fname):
                    os.remove(fname)


class StatCheck(unittest.TestCase):

    """This unit test takes care of verifying that statistics are generated correctly when
    considering four different inputs. Output will be checked against a standard file."""

    def test_stat(self):

        files = ["trinity.gtf",
                 "trinity.gff3",
                 "trinity.cDNA_match.gff3",
                 "trinity.match_matchpart.gff3"]
        files = [pkg_resources.resource_filename("Mikado.test", filename) for filename in files]

        std_lines = []
        with pkg_resources.resource_stream("Mikado.test", "trinity_stats.txt") as t_stats:
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


class ConfigureCheck(unittest.TestCase):

    """Test for creating configuration files"""

    __genomefile__ = None

    @classmethod
    def setUpClass(cls):
        cls.__genomefile__ = tempfile.NamedTemporaryFile(mode="wb", delete=False, suffix=".fa", prefix="configure")

        with pkg_resources.resource_stream("Mikado.test", "chr5.fas.gz") as _:
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
        namespace.reference = None
        namespace.external = None
        namespace.mode = "permissive"
        namespace.blast_targets = []
        namespace.junctions = []
        out = os.path.join(tempfile.gettempdir(), "configuration.yaml")
        with open(out, "w") as out_handle:
            namespace.out = out_handle
            Mikado.subprograms.configure.create_config(namespace)
        self.assertGreater(os.stat(out).st_size, 0)
        conf = Mikado.configuration.configurator.to_json(out)
        conf = Mikado.configuration.configurator.check_json(conf)
        conf = Mikado.configuration.configurator.check_json(conf)
        os.remove(out)

    @unittest.skipUnless((sys.version_info.minor>4),
                         "Due to a bug in JSONSCHEMA, Daijin configure fails with Python versions lower than 3.5.")
    def test_daijin_config(self):

        # Check the basic function actually functions
        _ = daijin_configurator.create_daijin_base_config()

        namespace = Namespace(default=False)
        namespace.r1 = []
        namespace.r2 = []
        namespace.samples = []
        namespace.strandedness = []
        namespace.asm_methods = ["cufflinks"]
        namespace.aligners = ["hisat"]
        namespace.modes = ["nosplit"]
        namespace.cluster_config = None
        namespace.scheduler = ""
        namespace.flank = None
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
        json_conf["db_settings"]["db"] = pkg_resources.resource_filename("Mikado.test", "mikado.db")

        json_conf["pick"]["files"]["input"] = pkg_resources.resource_filename("Mikado.test",
                                                                              "mikado_prepared.gtf")
        json_conf["pick"]["files"]["output_dir"] = tempfile.gettempdir()
        json_conf["pick"]["files"]["subloci_out"] = "mikado.subloci.gff3"
        json_conf["pick"]["files"]["monoloci_out"] = "mikado.monoloci.gff3"
        pick_caller = picker.Picker(json_conf=json_conf)
        with self.assertRaises(SystemExit):
            pick_caller()
        self.assertTrue(os.path.exists(os.path.join(tempfile.gettempdir(), "mikado.loci.gff3")))
        os.remove(os.path.join(tempfile.gettempdir(), "mikado.loci.gff3"))

    def test_multi_proc(self):
        json_conf = configurator.to_json(None)
        json_conf["pick"]["run_options"]["procs"] = 2
        json_conf["pick"]["files"]["input"] = pkg_resources.resource_filename("Mikado.test",
                                                                              "mikado_prepared.gtf")
        json_conf["pick"]["files"]["output_dir"] = tempfile.gettempdir()
        json_conf["pick"]["files"]["subloci_out"] = "mikado.subloci.gff3"
        json_conf["pick"]["files"]["monoloci_out"] = "mikado.monoloci.gff3"
        json_conf["db_settings"]["db"] = pkg_resources.resource_filename("Mikado.test", "mikado.db")
        pick_caller = picker.Picker(json_conf=json_conf)
        with self.assertRaises(SystemExit):
            pick_caller()
        self.assertTrue(os.path.exists(os.path.join(tempfile.gettempdir(), "mikado.loci.gff3")))
        os.remove(os.path.join(tempfile.gettempdir(), "mikado.loci.gff3"))


if __name__ == "__main__":
    unittest.main()
