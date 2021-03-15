#!/usr/bin/env python3
import itertools
import re

import pkg_resources

from Mikado.configuration import MikadoConfiguration, parse_list_file
from Mikado.utilities.log_utils import LoggingConfiguration, create_logger_from_conf, create_default_logger
from Mikado import utilities
import unittest
import os
import tempfile
import logging
import queue


class UtilTester(unittest.TestCase):

    def test_path(self):

        nonabs = os.path.basename(__file__)
        absolute = os.path.abspath(__file__)
        cur_dir = os.path.basename(os.path.dirname(__file__))

        self.assertEqual(os.path.join(cur_dir, nonabs),
                         utilities.path_join(cur_dir, nonabs))
        self.assertEqual(absolute, utilities.path_join(".", absolute))

    def test_memoizer(self):

        @utilities.memoize
        def tester(num):

            return num ** 2

        self.assertEqual(tester.cache, {})
        tester(10)
        self.assertIn("(10,){}", tester.cache)
        self.assertNotIn("(20,){}", tester.cache)
        tester(20)
        self.assertIn("(20,){}", tester.cache)
        tester.cache.clear()
        self.assertNotIn("(20,){}", tester.cache)
        self.assertNotIn("(10,){}", tester.cache)

    def test_check_to_regions(self):
        for invalid in [0, b"Chr1:100.1000", "Chr1:100.1000", "Chr1:1500-1000", None, "Chr1::1400-2000"]:
            with self.assertRaises(ValueError), self.subTest(invalid=invalid):
                utilities.to_region(invalid)

        for reg, result in [("Chr1:1000-2000", ("Chr1", 1000, 2000)),
                            ("ath_Chr1:1500-10000", ("ath_Chr1", 1500, 10000)),
                            ("hello_chrom:5000..50000", ("hello_chrom", 5000, 50000))]:
            chrom, start, end = utilities.to_region(reg)
            self.assertEqual(result, (chrom, start, end))
            chrom, start, end = utilities.to_region(reg.encode())
            self.assertEqual(result, (chrom, start, end))

        chroms = ("Chr1", "Chr4", "Chr5")
        starts = (100, 500, 1000, 5000)
        ends = (500, 800, 1500, 10000)
        ini_seps = (":", "-")
        seps = ("-", "..", "foo")

        for chrom, start, end, ini_sep, sep in itertools.product(chroms, starts, ends, ini_seps, seps):
            to_fail = (end < start) or (sep not in ("-", "..")) or (ini_sep != ":")
            string = "{chrom}{ini_sep}{start}{sep}{end}".format(**locals())
            with self.subTest(chrom=chrom, start=start, end=end, sep=sep):
                if to_fail:
                    with self.assertRaises(ValueError):
                        utilities.to_region(string)
                else:
                    uchrom, ustart, uend = utilities.to_region(string)
                    self.assertEqual(chrom, uchrom)
                    self.assertEqual(start, ustart)
                    self.assertEqual(end, uend)

    def test_percentage(self):
        for incorrect in [None, "a", b"a"]:
            with self.assertRaises((TypeError, ValueError)):
                utilities.percentage(incorrect)
        for incorrect in [-1, -0.1, 101, 10**4]:
            with self.assertRaises(ValueError):
                utilities.percentage(incorrect)
        for num in [0.4, 0.3, 30, 36, 80.4, 100, 1.1, 1]:
            unum = utilities.percentage(num)
            if num > 1:
                self.assertEqual(unum, num / 100)
            else:
                self.assertEqual(unum, num)
            self.assertIsInstance(unum, float)

    def test_grouper(self):

        objects = list(range(9))

        grouped = list(utilities.grouper(objects, 3))
        self.assertEqual(grouped, [[0, 1, 2], [3, 4, 5], [6, 7, 8]])
        grouped = list(utilities.grouper(objects, 4))
        self.assertEqual(grouped, [[0, 1, 2, 3], [4, 5, 6, 7], [8]])

    def test_merger(self):

        with tempfile.NamedTemporaryFile(suffix=".tmp", mode="wt", delete=False) as first, \
                tempfile.NamedTemporaryFile(suffix=".tmp", mode="wt", delete=False) as second, \
                tempfile.NamedTemporaryFile(suffix=".tmp", mode="wt", delete=False) as third,\
                tempfile.NamedTemporaryFile(suffix=".out", mode="wt", delete=False) as out:
            print("1/first case", file=first)
            print("3/third case", file=first)
            print("4/fourth case", file=first)
            print("2/second case", file=second)
            print("5/fifth case", file=third)
            print("6/sixth case", file=second)
            first.flush()
            second.flush()
            third.flush()
            logger = create_default_logger("test_merger", level="DEBUG")
            utilities.merge_partial([first.name, second.name, third.name], out, logger=logger)
            with open(out.name) as _:
                lines = [l for l in _]
            self.assertEqual(len(lines), 6, (out.name, lines))
            self.assertEqual(lines[0], "first case\n", lines)
            self.assertEqual(lines[1], "second case\n")
            self.assertEqual(lines[2], "third case\n")
            self.assertEqual(lines[3], "fourth case\n")
            self.assertEqual(lines[4], "fifth case\n")
            self.assertEqual(lines[5], "sixth case\n")
            # Verify the temporary files have been deleted correctly
            existing = []
            for handle in (first, second, third):
                if os.path.exists(handle.name):
                    existing.append(handle.name)
            self.assertEqual(len(existing), 0, [os.remove(name) for name in existing])
            self.assertTrue(os.path.exists(out.name))
            os.remove(out.name)


class ParseListFileTest(unittest.TestCase):

    def setUp(self) -> None:
        self.config = MikadoConfiguration()
        self.ref = pkg_resources.resource_filename("Mikado.tests", "reference.gff3")
        self.trinity = pkg_resources.resource_filename("Mikado.tests", "trinity.bed12")
        self.pacbio = pkg_resources.resource_filename("Mikado.tests", "pacbio.bam")

        self.files = {
            "ref": {"filename": self.ref,
                    "label": "ref", "stranded": True, "score": 10, "exclude_redundant": False,
                    "strip_cds": False, "skip_split": True, "reference": True},
            "pacbio": {"filename": self.pacbio,
                       "label": "pacbio", "stranded": True, "score": 5, "exclude_redundant": True,
                        "strip_cds": False, "skip_split": True, "reference": False},
            "trinity": {"filename": self.trinity,
                        "label": "trinity", "stranded": False, "score": -1, "exclude_redundant": True,
                        "strip_cds": True, "skip_split": True, "reference": False},
        }

    def test_parse_failing_basic_list_file(self):
        # 0 Filename
        # 1 Label
        # 2 Strandedness (True/False)
        # 3 Score (Float)
        # 4 Reference
        # 5 Exclude redundant
        # 6 strip cds
        # 7 skip split

        with tempfile.NamedTemporaryFile(mode="wt") as fail:
            print(*[self.files[_]["filename"] for _ in self.files], sep="\n", file=fail)
            fail.flush()
            with self.assertRaises(IndexError):
                parse_list_file(self.config, fail.name)

        with tempfile.NamedTemporaryFile("wt") as fail:
            for key, item in self.files.items():
                print(item["filename"], item["label"], file=fail, sep="\t")
            fail.flush()
            with self.assertRaises(IndexError):
                parse_list_file(self.config, fail.name)

        with tempfile.NamedTemporaryFile("wt") as fail:
            for key, item in self.files.items():
                print(item["filename"], item["label"], "foo", file=fail, sep="\t")
            fail.flush()
            with self.assertRaises(ValueError):
                parse_list_file(self.config, fail.name)

        with tempfile.NamedTemporaryFile("wt") as fail:
            for key, item in self.files.items():
                # Double labels
                print(item["filename"], self.files["pacbio"]["label"], item["stranded"], file=fail, sep="\t")
            fail.flush()
            with self.assertRaises(ValueError):
                parse_list_file(self.config, fail.name)

        with tempfile.NamedTemporaryFile("wt") as fail, tempfile.TemporaryDirectory() as tmp:
            for key, item in self.files.items():
                fname = os.path.join(tmp, os.path.basename(item["filename"]))
                self.assertFalse(os.path.exists(fname))
                print(fname, item["label"], item["stranded"], file=fail, sep="\t")
            fail.flush()
            with self.assertRaises(ValueError):
                parse_list_file(self.config, fail.name)

        with tempfile.NamedTemporaryFile("wt") as fail:
            for key, item in self.files.items():
                print(item["filename"], item["label"], item["stranded"], "wrong", file=fail, sep="\t")
            fail.flush()
            with self.assertRaises(ValueError) as exc:
                parse_list_file(self.config, fail.name)
            self.assertIsNotNone(re.search(r"Invalid score specified for", str(exc.exception)))

    def test_valid(self):
        # 0 Filename
        # 1 Label
        # 2 Strandedness (True/False)
        # 3 Score (Float)
        # 4 Reference
        # 5 Exclude redundant
        # 6 strip cds
        # 7 skip split

        with tempfile.NamedTemporaryFile(mode="wt") as out:
            for _, item in self.files.items():
                print(item["filename"], item["label"], item["stranded"], item["score"],
                      item["reference"], item["exclude_redundant"], item["strip_cds"], item["skip_split"],
                      sep="\t", file=out)
            out.flush()
            configuration = parse_list_file(self.config, out.name)

        for arr in (configuration.prepare.files.source_score,
                    configuration.prepare.files.reference,
                    configuration.prepare.files.exclude_redundant,
                    configuration.prepare.files.reference,
                    configuration.prepare.files.strip_cds,
                    configuration.pick.chimera_split.skip):
            self.assertEqual(len(arr), 3)

        for _, item in self.files.items():
            label = item["label"]
            self.assertIn(item["filename"], configuration.prepare.files.gff)
            self.assertIn(item["label"], configuration.prepare.files.labels)
            pos = configuration.prepare.files.gff.index(item["filename"])
            self.assertEqual(pos,configuration.prepare.files.labels.index(item["label"]))
            if item["stranded"] is True:
                self.assertIn(item["filename"], configuration.prepare.files.strand_specific_assemblies)
            else:
                self.assertNotIn(item["filename"], configuration.prepare.files.strand_specific_assemblies)
            self.assertEqual(item["score"], configuration.prepare.files.source_score[label])
            self.assertEqual(item["exclude_redundant"], configuration.prepare.files.exclude_redundant[pos])
            self.assertEqual(item["reference"], configuration.prepare.files.reference[pos])
            self.assertEqual(item["strip_cds"], configuration.prepare.files.strip_cds[pos])
            self.assertEqual(item["skip_split"], configuration.pick.chimera_split.skip[pos])

    def test_default(self):
        with tempfile.NamedTemporaryFile(mode="wt") as out:
            # Test that header lines are skipped
            print("# Filename label stranded", file=out)
            print("", file=out)
            for _, item in self.files.items():
                print(item["filename"], item["label"], item["stranded"], sep="\t", file=out)
            out.flush()
            configuration = parse_list_file(self.config, out.name)

        for name, arr in zip(("source_score", "reference", "exclude_redundant", "strip_cds", "skip_split"),
                             (configuration.prepare.files.source_score,
                              configuration.prepare.files.reference,
                              configuration.prepare.files.exclude_redundant,
                              configuration.prepare.files.strip_cds,
                              configuration.pick.chimera_split.skip)):
            self.assertEqual(len(arr), 3, name)

        for _, item in self.files.items():
            label = item["label"]
            self.assertIn(item["filename"], configuration.prepare.files.gff)
            self.assertIn(item["label"], configuration.prepare.files.labels)
            pos = configuration.prepare.files.gff.index(item["filename"])
            self.assertEqual(pos,configuration.prepare.files.labels.index(item["label"]))
            if item["stranded"] is True:
                self.assertIn(item["filename"], configuration.prepare.files.strand_specific_assemblies)
            else:
                self.assertNotIn(item["filename"], configuration.prepare.files.strand_specific_assemblies)

        # for arr, pos, default in [("reference", 4, False), ("exclude_redundant", 5, False),
        #                                   ("strip_cds", 6, False), ("skip_split", 7, False)]:

        self.assertEqual(configuration.prepare.files.source_score,
                         dict((label, 0) for label in configuration.prepare.files.labels))
        self.assertEqual(configuration.prepare.files.reference, [False] * 3)
        self.assertEqual(configuration.prepare.files.exclude_redundant, [False] * 3)
        self.assertEqual(configuration.prepare.files.strip_cds, [False] * 3)
        self.assertEqual(configuration.pick.chimera_split.skip, [False] * 3)


class LogUtilsTester(unittest.TestCase):

    def test_null_log(self):
        """Tester for the default null logger. The function returns
        a single null logger created at run time."""

        null_logger = utilities.log_utils.create_null_logger()
        self.assertIsInstance(null_logger, logging.Logger)
        self.assertEqual(null_logger.name, "null")
        null2 = utilities.log_utils.create_null_logger("test", "will_nilly",
                                                              useless_val=None)
        self.assertIsNot(null2, null_logger)
        self.assertIsInstance(null2.handlers[0], logging.NullHandler)
        self.assertIsInstance(null_logger.handlers[0], logging.NullHandler)

    def test_default_logger(self):

        """
        Tester for the default non-null logger. The function returns
        a logger created on-demand with the given name.
        :return:
        """

        first_logger = utilities.log_utils.create_default_logger("first_case")
        self.assertEqual(first_logger.name, "first_case")
        self.assertIsInstance(first_logger, logging.Logger)
        with self.assertLogs("first_case", level="WARNING") as test_log:
            first_logger.warning("Test log")
        self.assertIn("WARNING:first_case:Test log", test_log.output)

        second_logger = utilities.log_utils.create_default_logger("second_case")
        self.assertIsNot(first_logger, second_logger)
        self.assertEqual(second_logger.name, "second_case")
        self.assertEqual(second_logger.level, 30)
        self.assertEqual(second_logger.level, first_logger.level)
        with self.assertLogs("second_case", level="DEBUG") as test_log:
            second_logger.debug("Test log")
        self.assertIn("DEBUG:second_case:Test log", test_log.output)

    def test_queue_logger_error(self):

        class inc_obj:
            def __init__(self):
                pass

        with self.assertRaises(AttributeError):
            inc_instance = inc_obj()
            utilities.log_utils.create_queue_logger(inc_instance)

    def test_queue_logger(self):

        class obj:
            def __init__(self):
                self.logging_queue = queue.Queue()
                self.log_level = "DEBUG"

        instance = obj()

        utilities.log_utils.create_queue_logger(instance, prefix="prefix_test")

        self.assertTrue(hasattr(instance, "logger"))
        self.assertTrue(hasattr(instance, "_log_handler"))
        self.assertEqual(instance.logger.level, 10)
        self.assertEqual(instance.logger.name, "prefix_test.default")

        instance.name = "test"
        utilities.log_utils.create_queue_logger(instance)
        self.assertTrue(hasattr(instance, "logger"))
        self.assertTrue(hasattr(instance, "_log_handler"))
        self.assertEqual(instance.logger.level, 10)
        self.assertEqual(instance.logger.name, "test")

    def test_queue_logger_no_default_level(self):

        class obj:
            def __init__(self):
                self.logging_queue = queue.Queue()

        instance = obj()
        utilities.log_utils.create_queue_logger(instance, prefix="prefix_test")
        self.assertEqual(instance.logger.level, 30)

    def test_check_logger(self):
        logger = logging.getLogger("test_validity")
        self.assertEqual(logger, utilities.log_utils.check_logger(logger))
        logger = "foo"
        with self.assertRaises(ValueError):
            utilities.log_utils.check_logger(logger)

    def test_merge_dictionaries(self):

        a = {10: 20, 30: 50}
        b = {20: 40, 30: 10}

        merged = utilities.merge_dictionaries(a, b)
        self.assertEqual(merged,
                         {10: 20, 20:40,
                          30: 10},  # 30 is updated to the value of dictionary b
                         merged)

    def test_logger_from_conf(self):
        with tempfile.TemporaryDirectory() as folder:
            for logger_name, (mode, level, initialiser, streaming) in enumerate(
                    itertools.product(["a", "w"], ["WARNING", "INFO", "DEBUG"],
                                      (LoggingConfiguration, MikadoConfiguration),
                                      (False, True))):
                with tempfile.NamedTemporaryFile(mode="wt", dir=folder, suffix=".log") as log:
                    if streaming is True:
                        name = None
                    else:
                        name = log.name
                    if initialiser == LoggingConfiguration:
                        conf = LoggingConfiguration()
                        conf.log = name
                        conf.log_level = level
                    else:
                        conf = MikadoConfiguration()
                        conf.log_settings.log = name
                        conf.log_settings.log_level = level

                logger = create_logger_from_conf(conf, name="test_logger_from_conf" + str(logger_name), mode=mode)
                self.assertTrue(os.path.exists(log.name) or streaming)
                self.assertTrue(os.path.exists(log.name) != streaming)  # Either we have created a streaming log or not
                self.assertEqual(logging.getLevelName(logger.level), level, conf)
                if streaming:
                    self.assertIsInstance(logger.handlers[0], logging.StreamHandler)
                else:
                    self.assertIsInstance(logger.handlers[0], logging.FileHandler)
                    self.assertEqual(logger.handlers[0].mode, mode)
                if streaming is False and os.path.exists(log.name):
                    os.remove(log.name)


if __name__ == "__main__":
    unittest.main()
