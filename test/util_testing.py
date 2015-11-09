import unittest
import mikado_lib.utilities
import logging

__author__ = 'Luca Venturini'

class TestUtils(unittest.TestCase):

    def test_check_logger(self):

        logger = logging.getLogger("test_validity")
        self.assertEqual(logger, mikado_lib.utilities.log_utils.check_logger(logger))
        logger = "foo"
        with self.assertRaises(ValueError):
            mikado_lib.utilities.log_utils.check_logger(logger)

if __name__ == "__main__":
    unittest.main()
