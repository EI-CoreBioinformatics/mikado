__author__ = 'Luca Venturini'


import unittest
from Mikado.loci import Transcript
import intervaltree
from Mikado.utilities.log_utils import create_default_logger


class TestMetricsEndDistances(unittest.TestCase):

    logger = create_default_logger("End")
    logger.setLevel("ERROR")

    def setUp(self):

        self.tr = Transcript()
        self.tr.logger = self.logger
        self.tr.start = 101
        self.tr.end = 10000
        self.tr.add_exons([(101, 300),
                           (501, 800),
                           (1001, 1200),
                           (1301, 2000),
                           (3501, 5000),
                           (5501, 6000),
                           (6201, 7000),
                           (7301, 7700),
                           (8201, 9000),
                           (9101, 9300),
                           (9501, 9700),
                           (9801, 10000)])
        self.tr.id = "test1"
        self.tr.parent = "test1.gene"

    def test_end_positive(self):

        self.tr.strand = "+"

        cds = [(1161, 1200),  # 40 % 3 == 1
               (1301, 2000),  # 700 % 3 == 1
               (3501, 5000),  # 1500 % 3 == 0
               (5501, 6000),  # 500 % 3 == 2
               (6201, 7000),  # 800 % 3 == 2
               (7301, 7700),  # 400 % 3 == 1
               (8201, 9000),  # 800 % 3 == 2
               (9101, 9130)]

        self.tr.add_exons(cds, features="CDS")
        self.tr.finalize()
        self.assertEqual(self.tr.selected_cds_end,
                         self.tr.combined_cds_end)
        self.assertEqual(self.tr.selected_cds_end,
                         9130)
        self.assertEqual(self.tr.end_distance_from_junction,
                         (9300 - 9131 + 1) + (9700 - 9501 + 1)
                         )
        self.assertEqual(self.tr.end_distance_from_tes,
                         (9300 - 9131 + 1) + (9700 - 9501 + 1) + (10000 - 9801 + 1)
                         )

        self.tr.strip_cds()
        self.assertEqual(len(self.tr.internal_orfs), 0, self.tr.internal_orfs)
        self.tr.finalized = False
        cds = [(1161, 1200),  # 40 % 3 == 1
               (1301, 2000),  # 700 % 3 == 1
               (3501, 5000),  # 1500 % 3 == 0
               (5501, 6000),  # 500 % 3 == 2
               (6201, 7000),  # 800 % 3 == 2
               (7301, 7700),  # 400 % 3 == 1
               (8201, 9000),  # 800 % 3 == 2
               (9101, 9300),  # 200 % 3 == 2
               (9501, 9690)  # 190 % 3 == 1
               ]
        self.tr.add_exons(cds, features="CDS")

        self.tr.finalize()
        self.assertEqual(self.tr.combined_cds_end,
                         9690)
        self.assertEqual(self.tr.selected_cds_end,
                         self.tr.combined_cds_end)

        self.assertEqual(self.tr.end_distance_from_junction,
                         (9700 - 9691 + 1)
                         )
        self.assertEqual(self.tr.end_distance_from_tes,
                         (9700 - 9691 + 1) + (10000 - 9801 + 1)
                         )

        self.tr.strip_cds()
        self.assertEqual(self.tr.combined_cds_end,
                         self.tr.selected_cds_end,
                         self.tr.combined_cds)
        self.assertEqual(self.tr.combined_cds_end,
                         None,
                         self.tr.combined_cds_end)

        self.tr.finalized = False
        cds = [(1161, 1200),  # 40 % 3 == 1
               (1301, 2000),  # 700 % 3 == 1
               (3501, 5000),  # 1500 % 3 == 0
               (5501, 6000),  # 500 % 3 == 2
               (6201, 7000),  # 800 % 3 == 2
               (7301, 7700),  # 400 % 3 == 1
               (8201, 9000),  # 800 % 3 == 2
               (9101, 9300),  # 200 % 3 == 2
               (9501, 9700),  # 200 % 3 == 2
               (9801, 9820),  # 20 % 2 == 2
               ]
        self.tr.add_exons(cds, features="CDS")

        self.tr.finalize()
        self.assertEqual(self.tr.combined_cds_end,
                         9820)
        self.assertEqual(self.tr.selected_cds_end,
                         self.tr.combined_cds_end)
        self.assertEqual(self.tr.end_distance_from_tes,
                         180)
        self.assertEqual(self.tr.end_distance_from_junction,
                         0)

        pass

if __name__ == "__main__":
    unittest.main()