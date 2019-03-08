#!/usr/bin/env python3

from .. import utilities
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

    def test_grouper(self):

        objects = list(range(9))

        grouped = list(utilities.grouper(objects, 3))
        self.assertEqual(grouped, [[0, 1, 2], [3, 4, 5], [6, 7, 8]])
        grouped = list(utilities.grouper(objects, 4))
        self.assertEqual(grouped, [[0, 1, 2, 3], [4, 5, 6, 7], [8]])

    def test_merger(self):

        first_name = tempfile.mktemp(suffix=".tmp", dir=tempfile.tempdir)
        second_name = tempfile.mktemp(suffix=".tmp", dir=tempfile.tempdir)
        third_name = tempfile.mktemp(suffix=".tmp", dir=tempfile.tempdir)

        with open(first_name, "wt") as first, open(second_name, "wt") as second,\
                open(third_name, "wt") as third:

            print("1/first case", file=first)
            print("3/third case", file=first)
            print("4/fourth case", file=first)
            print("2/second case", file=second)
            print("5/fifth case", file=third)
            print("6/sixth case", file=second)

        out_name = tempfile.mktemp(suffix=".out", dir=tempfile.tempdir)
        with open(out_name, "wt") as out:
            utilities.merge_partial([first.name,
                                            second.name,
                                            third.name], out)
        with open(out_name) as out:
            lines = [l for l in out]
        self.assertEqual(lines[0], "first case\n")
        self.assertEqual(lines[1], "second case\n")
        self.assertEqual(lines[2], "third case\n")
        self.assertEqual(lines[3], "fourth case\n")
        self.assertEqual(lines[4], "fifth case\n")
        self.assertEqual(lines[5], "sixth case\n")
        # Verify the temporary files have been deleted correctly
        self.assertFalse(os.path.exists(first_name))
        self.assertFalse(os.path.exists(second_name))
        self.assertFalse(os.path.exists(third_name))
        self.assertTrue(os.path.exists(out_name))
        os.remove(out_name)
        self.assertFalse(os.path.exists(out_name))


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


if __name__ == "__main__":
    unittest.main()
