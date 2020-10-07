#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for the scales library
"""

import unittest
from .. import loci, utilities
from .. import scales
import argparse
import os
from .. import parsers
import csv
import tempfile
from ..scales.assignment.assigner import Assigner


class AssignerTest(unittest.TestCase):
    """
    This unit test has the purpose of testing the scales module of Mikado.py.
    """

    def test_get_f1(self):

        # __slots__ = ["ref_id", "ref_gene", "ccode",
        #              "tid", "gid",
        #              "tid_num_exons", "ref_num_exons",
        #              "n_prec", "n_recall", "n_f1",
        #              "j_prec", "j_recall", "j_f1",
        #              "e_prec", "e_recall", "e_f1",
        #              "distance"]

        result_perfect = scales.resultstorer.ResultStorer(
            "t1", "g1", "=",
            "p1", "pg1", "2", "2",
            100, 100, 100,
            100, 100, 100,
            100, 100, 100,
            0, "chr1:100..10000")

        result_perfect_j = scales.resultstorer.ResultStorer(
            "t1", "g1", "=",
            "p1", "pg1", "2", "2",
            80, 80, 80,
            100, 100, 100,
            0, 0, 0,
            0, "chr1:100..10000")

        # This does not make any sense, but it's only for the tests
        result_perfect_n = scales.resultstorer.ResultStorer(
            "t1", "g1", "j",
            "p1", "pg1", "2", "2",
            100, 100, 100,
            80, 80, 80,
            0, 0, 0,
            0, "chr1:100..10000")

        result_imperfect = scales.resultstorer.ResultStorer(
            "t1", "g1", "j",
            "p1", "pg1", "2", "2",
            80, 80, 80,
            80, 80, 80,
            0, 0, 0,
            0, "chr1:100..10000")

        result_near = scales.resultstorer.ResultStorer(
            "t1", "g1", "p",
            "p1", "pg1", "2", "2",
            0, 0, 0,
            0, 0, 0,
            0, 0, 0,
            10, "chr1:100..10000")

        result_middle = scales.resultstorer.ResultStorer(
            "t1", "g1", "p",
            "p1", "pg1", "2", "2",
            0, 0, 0,
            0, 0, 0,
            0, 0, 0,
            500, "chr1:100..10000")

        result_far = scales.resultstorer.ResultStorer(
            "t1", "g1", "p",
            "p1", "pg1", "2", "2",
            0, 0, 0,
            0, 0, 0,
            0, 0, 0,
            1000, "chr1:100..10000")

        result_x = scales.resultstorer.ResultStorer(
            "t1", "g1", "x",
            "p1", "pg1", "2", "2",
            100, 100, 100,
            100, 100, 100,
            100, 100, 100,
            1000, "chr1:100..10000")

        self.assertEqual(sorted(
            [result_perfect, result_imperfect, result_perfect_j, result_perfect_n],
            key=Assigner.get_f1, reverse=True),
            [result_perfect, result_perfect_j, result_perfect_n, result_imperfect]
        )

        self.assertEqual(sorted([result_perfect, result_far], key=Assigner.get_f1, reverse=True),
                         [result_perfect, result_far])

        self.assertEqual(sorted(
            [result_far, result_near, result_middle, result_imperfect, result_perfect],
            key=Assigner.get_f1, reverse=True),
            [result_perfect, result_imperfect, result_near, result_middle,  result_far]
        )

        self.assertEqual(sorted(
            [result_perfect, result_x],
            key=Assigner.get_f1, reverse=True),
            [result_perfect, result_x]
        )

    def test_self(self):
        reference = loci.Transcript()
        reference.start = 100
        reference.end = 2000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.add_exons([(100, 300), (500, 1000), (1500, 2000)])
        reference.finalize()

        result, _ = Assigner.compare(reference, reference)
        self.assertEqual(result.ccode, ("=",))
        self.assertEqual(result.n_f1, (100,))
        self.assertEqual(result.n_prec, (100,))
        self.assertEqual(result.n_recall, (100,))
        self.assertEqual(result.j_f1, (100,))
        self.assertEqual(result.j_prec, (100,))
        self.assertEqual(result.j_recall, (100,))

    def test_mono_self(self):
        reference = loci.Transcript()
        reference.start = 100
        reference.end = 1000
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.strand = "+"
        reference.exons = [(100, 1000)]
        reference.finalize()

        result, _ = Assigner.compare(reference, reference)
        self.assertEqual(result.ccode, ("_",))
        self.assertEqual(result.n_f1, (100,))
        self.assertEqual(result.n_prec, (100,))
        self.assertEqual(result.n_recall, (100,))
        self.assertEqual(result.j_f1, (100,))
        self.assertEqual(result.j_prec, (100,))
        self.assertEqual(result.j_recall, (100,))

    def test_equal(self):
        reference = loci.Transcript()
        reference.start = 100
        reference.end = 2000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(100, 300), (500, 1000), (1500, 2000)]
        reference.finalize()

        prediction = loci.Transcript()
        prediction.start = 200
        prediction.end = 1800
        prediction.strand = "+"
        prediction.chrom = "Chr1"
        prediction.id = "P1.1"
        prediction.parent = "P1"
        prediction.exons = [(200, 300), (500, 1000), (1500, 1800)]
        prediction.finalize()

        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("=",))
        self.assertEqual(result.j_f1, (100,))
        self.assertEqual(result.j_prec, (100,))
        self.assertEqual(result.j_recall, (100,))
        self.assertAlmostEqual(result.n_prec[0], 100, delta=0.1)
        self.assertAlmostEqual(result.n_recall[0],
                               100 * (300 - 200 + 1 + 1000 - 500 + 1 + 1800 - 1500 + 1) / reference.cdna_length,
                               delta=0.1)

        prediction.strand = "-"
        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("x",))

    def test_mono_equal(self):
        reference = loci.Transcript()
        reference.start = 100
        reference.end = 1000
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.strand = "+"
        reference.exons = [(100, 1000)]
        reference.finalize()

        prediction = loci.Transcript()
        prediction.start = 105
        prediction.end = 995
        prediction.chrom = "Chr1"
        prediction.id = "G1.1"
        prediction.parent = "G1"
        prediction.strand = "+"
        prediction.exons = [(105, 995)]
        prediction.finalize()

        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("_",))
        self.assertEqual(result.j_f1, (100,))
        self.assertEqual(result.j_prec, (100,))
        self.assertEqual(result.j_recall, (100,))
        self.assertAlmostEqual(result.n_prec[0], 100, delta=0.1)
        self.assertAlmostEqual(result.n_recall[0], 100 * (995 - 105 + 1) / reference.cdna_length, delta=0.1)

        prediction.strand = "-"
        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("x",))

    def test_mono_semiequal(self):
        reference = loci.Transcript()
        reference.start = 100
        reference.end = 1000
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.strand = "+"
        reference.exons = [(100, 1000)]
        reference.finalize()

        prediction = loci.Transcript()
        prediction.start = 200
        prediction.end = 500
        prediction.chrom = "Chr1"
        prediction.id = "G1.1"
        prediction.parent = "G1"
        prediction.strand = "+"
        prediction.exons = [(200, 500)]
        prediction.finalize()

        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("c",))
        self.assertEqual(result.j_f1, (100,))
        self.assertEqual(result.j_prec, (100,))
        self.assertEqual(result.j_recall, (100,))
        self.assertAlmostEqual(result.n_prec[0], 100, delta=0.1)
        self.assertAlmostEqual(result.n_recall[0], 100 * (500 - 200 + 1) / reference.cdna_length,
                               delta=0.1)

        prediction.strand = "-"
        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("x",))

    def test_mono_overlap(self):

        """Test that two monoexonic overlapping genes get a m"""

        reference = loci.Transcript()
        reference.start = 100
        reference.end = 1000
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.strand = "+"
        reference.exons = [(100, 1000)]
        reference.finalize()

        prediction = loci.Transcript()
        prediction.start = 200
        prediction.end = 3000
        prediction.chrom = "Chr1"
        prediction.id = "G1.1"
        prediction.parent = "G1"
        prediction.strand = "+"
        prediction.exons = [(200, 3000)]
        prediction.finalize()

        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("m",))
        self.assertEqual(result.j_f1, (100,))
        self.assertEqual(result.j_prec, (100,))
        self.assertEqual(result.j_recall, (100,))
        self.assertAlmostEqual(result.n_prec[0], 100 * 801 / prediction.cdna_length, delta=0.1)
        self.assertAlmostEqual(result.n_recall[0], 100 * 801 / reference.cdna_length, delta=0.1)

    def test_contained(self):

        reference = loci.Transcript()
        reference.start = 100
        reference.end = 2000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(100, 300), (500, 1000), (1500, 2000)]
        reference.finalize()

        prediction = loci.Transcript()
        prediction.start = 600
        prediction.end = 1800
        prediction.strand = "+"
        prediction.chrom = "Chr1"
        prediction.id = "P1.1"
        prediction.parent = "P1"
        prediction.exons = [(600, 1000), (1500, 1800)]
        prediction.finalize()

        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("c",))
        self.assertAlmostEqual(result.j_f1[0], (100 * (2 / 3),)[0], delta=0.1)
        self.assertAlmostEqual(result.j_prec, (100,))
        self.assertAlmostEqual(result.j_recall, (100 * 1 / 2,), delta=0.1)
        self.assertAlmostEqual(result.n_prec[0], 100, delta=0.1)

        prediction.strand = "-"
        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("x",))

    def test_exon_skipping(self):

        """Test to verify that we call correctly an exon skipping as j"""

        reference = loci.Transcript()
        reference.chrom = 12
        reference.strand = "+"
        reference.add_exons([(7014003, 7014067),
                             (7014749, 7014923),
                             (7015573, 7015826),
                             (7016479, 7016609),
                             (7019054, 7019190),
                             (7021894, 7022191),
                             (7023055, 7023555),
                             (7015009, 7015118)  # This is our additional internal exon
                             ])
        reference.start = min(_[0] for _ in reference.exons)
        reference.end = max(_[1] for _ in reference.exons)
        reference.id = "t5"

        reference.finalize()

        prediction = loci.Transcript()
        prediction.chrom, prediction.strand = reference.chrom, reference.strand
        prediction.start, prediction.end = reference.start, reference.end
        prediction.add_exons([(7014003, 7014067),
                             (7014749, 7014923),
                             (7015573, 7015826),
                             (7016479, 7016609),
                             (7019054, 7019190),
                             (7021894, 7022191),
                             (7023055, 7023555)
                              ])
        prediction.finalize()
        prediction.id = "t3"

        self.assertEqual(len(reference.introns), 7)
        self.assertEqual(len(prediction.introns), 6)

        ref_vs_pred, _ = Assigner.compare(prediction, reference)
        self.assertEqual(ref_vs_pred.ccode, ("j", ), prediction.introns - reference.introns)
        pred_vs_ref, _ = Assigner.compare(reference, prediction)
        self.assertEqual(pred_vs_ref.ccode, ("j", ), prediction.introns - reference.introns)

    def test_alternative(self):

        reference = loci.Transcript()
        reference.start = 100
        reference.end = 2000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(100, 300), (500, 1000), (1500, 2000)]
        reference.finalize()

        prediction = loci.Transcript()
        prediction.start = 100
        prediction.end = 2000
        prediction.strand = "+"
        prediction.chrom = "Chr1"
        prediction.id = "P1.1"
        prediction.parent = "P1"
        prediction.exons = [(100, 300), (600, 1000), (1500, 2000)]
        prediction.finalize()

        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("j",))
        self.assertAlmostEqual(result.j_f1[0], (100 * ((3 / 4) / 1),)[0], delta=0.1)
        self.assertAlmostEqual(result.j_prec, (100 * 3 / 4,), delta=0.1)
        self.assertAlmostEqual(result.j_recall, (100 * 3 / 4,), delta=0.1)

        prediction.strand = "-"
        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("X",))

    def test_mono_intronic(self):

        """
        R   |xxxxx|-------------------|xxxxxxxxx|
        P            |xxxxxxx|
        
        Expected class code: i
        """

        reference = loci.Transcript()
        reference.start = 100
        reference.end = 2000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(100, 300), (500, 1000), (1500, 2000)]
        reference.finalize()

        prediction = loci.Transcript()
        prediction.start = 1100
        prediction.end = 1400
        prediction.strand = "+"
        prediction.chrom = "Chr1"
        prediction.id = "P1.1"
        prediction.parent = "P1"
        prediction.exons = [(1100, 1400)]
        prediction.finalize()

        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("i",))
        self.assertEqual(result.j_f1, (0,))
        self.assertEqual(result.j_prec, (0,))
        self.assertEqual(result.j_recall, (0,))
        if result.j_f1[0] > 0:
            self.assertAlmostEqual(result.j_f1[0], 2 * (result.j_recall[0] * result.j_prec[0]) /
                                   (result.j_prec[0] + result.j_recall[0]), delta=0.1)
        self.assertEqual(result.n_prec, (0,))
        self.assertEqual(result.n_recall, (0,))

        prediction.strand = "-"
        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("i",))

    def test_multi_intronic(self):

        """
        R   |xxxxx|-------------------|xxxxxxxxx|
        P            |xxx|----|xxxx|
        
        OR
        
        R   |xxxxx|----------|xxxxxxxxx|---------|xxxxxxxxx|
        P            |xxx|----------------|xxxx|
        
        Expected class code: I
        """

        reference = loci.Transcript()
        reference.start = 100
        reference.end = 2000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(100, 300), (500, 1000), (1500, 2000)]
        reference.finalize()

        prediction = loci.Transcript()
        prediction.start = 1100
        prediction.end = 1400
        prediction.strand = "+"
        prediction.chrom = "Chr1"
        prediction.id = "P1.1"
        prediction.parent = "P1"
        prediction.exons = [(1100, 1200), (1300, 1400)]
        prediction.finalize()

        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("I",))
        self.assertEqual(result.j_f1, (0,))
        self.assertEqual(result.j_prec, (0,))
        self.assertEqual(result.j_recall, (0,))
        if result.j_f1[0] > 0:
            self.assertAlmostEqual(result.j_f1, 2 * (result.j_recall[0] * result.j_prec[0]) / (
                result.j_prec[0] + result.j_recall[0]), delta=0.1)
        self.assertEqual(result.n_prec, (0,))
        self.assertEqual(result.n_recall, (0,))

        prediction.strand = "-"
        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("I",))

        # Now the reference spans two introns
        prediction = loci.Transcript()
        prediction.start = 350
        prediction.end = 1400
        prediction.strand = "+"
        prediction.chrom = "Chr1"
        prediction.id = "P1.1"
        prediction.parent = "P1"
        prediction.exons = [(350, 450), (1300, 1400)]
        prediction.finalize()

        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("I",))
        self.assertEqual(result.j_f1, (0,))
        self.assertEqual(result.j_prec, (0,))
        self.assertEqual(result.j_recall, (0,))
        if result.j_f1[0] > 0:
            self.assertAlmostEqual(result.j_f1, 2 * (result.j_recall[0] * result.j_prec[0]) / (
                result.j_prec[0] + result.j_recall[0]), delta=0.1)
        self.assertEqual(result.n_prec, (0,))
        self.assertEqual(result.n_recall, (0,))

        prediction.strand = "-"
        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("I",))

    def test_overlap(self):

        """
        R   |xxxxx|-------|xxxxx|
        P      |xxxxx|---|xxxx|
        
        No junction in common
        
        """

        reference = loci.Transcript()
        reference.start = 100
        reference.end = 2000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(100, 300), (500, 1000), (1500, 2000)]
        reference.finalize()

        prediction = loci.Transcript()
        prediction.start = 200
        prediction.end = 1300
        prediction.strand = "+"
        prediction.chrom = "Chr1"
        prediction.id = "P1.1"
        prediction.parent = "P1"
        prediction.exons = [(200, 700), (900, 1300)]
        prediction.finalize()

        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("o",))
        self.assertEqual(result.j_f1, (0,))
        self.assertEqual(result.j_prec, (0,))
        self.assertEqual(result.j_recall, (0,))
        n_overlap = (300 - 200 + 1 + 700 - 500 + 1 + 1000 - 900 + 1)
        self.assertAlmostEqual(result.n_prec[0], 100 * n_overlap / prediction.cdna_length, delta=0.1)
        self.assertAlmostEqual(result.n_recall[0], 100 * n_overlap / reference.cdna_length, delta=0.1)
        n_f1 = 2 * (n_overlap ** 2 / (reference.cdna_length * prediction.cdna_length)) / (
            n_overlap / reference.cdna_length + n_overlap / prediction.cdna_length)
        self.assertAlmostEqual(result.n_f1[0], 100 * n_f1, delta=0.1)

        prediction.strand = "-"
        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("X",))

    def test_ccode_e(self):

        """Case:
        
        R ---|xxxxx|-------|xxxxxx|---
        P      |xxxxxx|
        
        Exonic and intronic overlap
        
        """

        reference = loci.Transcript()
        reference.start = 100
        reference.end = 2000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(100, 300), (500, 1000), (1500, 2000)]
        reference.finalize()

        prediction = loci.Transcript()
        prediction.start = 200
        prediction.end = 400
        prediction.strand = "+"
        prediction.chrom = "Chr1"
        prediction.id = "P1.1"
        prediction.parent = "P1"
        prediction.exons = [(200, 400)]
        prediction.finalize()

        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("e",), result)
        self.assertEqual(result.j_f1, (0,))
        self.assertEqual(result.j_prec, (0,))
        self.assertEqual(result.j_recall, (0,))
        n_overlap = (300 - 200 + 1)
        self.assertAlmostEqual(result.n_prec[0], 100 * n_overlap / prediction.cdna_length, delta=0.1)
        self.assertAlmostEqual(result.n_recall[0], 100 * n_overlap / reference.cdna_length, delta=0.1)
        n_f1 = 2 * (n_overlap ** 2 / (reference.cdna_length * prediction.cdna_length)) / (
            n_overlap / reference.cdna_length + n_overlap / prediction.cdna_length)
        self.assertAlmostEqual(result.n_f1[0], 100 * n_f1, delta=0.1)

        prediction.strand = "-"
        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("x",))

    def test_not_ccode_e(self):

        """Case:
        
        R ---|xxxxx|-------|xxxxxx|---
        P |xxxxxx|
        
        Exonic overlap only
        
        """
        reference = loci.Transcript()
        reference.start = 100
        reference.end = 2000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(100, 300), (500, 1000), (1500, 2000)]
        reference.finalize()

        prediction = loci.Transcript()
        prediction.start = 50
        prediction.end = 310
        prediction.strand = "+"
        prediction.chrom = "Chr1"
        prediction.id = "P1.1"
        prediction.parent = "P1"
        prediction.exons = [(50, 310)]
        prediction.finalize()

        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("g",))
        self.assertEqual(result.j_f1, (0,))
        self.assertEqual(result.j_prec, (0,))
        self.assertEqual(result.j_recall, (0,))
        n_overlap = (300 - 100 + 1)
        self.assertAlmostEqual(result.n_prec[0], 100 * n_overlap / prediction.cdna_length, delta=0.1)
        self.assertAlmostEqual(result.n_recall[0], 100 * n_overlap / reference.cdna_length, delta=0.1)
        n_f1 = 2 * (n_overlap ** 2 / (reference.cdna_length * prediction.cdna_length)) / (
            n_overlap / reference.cdna_length + n_overlap / prediction.cdna_length)
        self.assertAlmostEqual(result.n_f1[0], 100 * n_f1, delta=0.1)

        prediction.strand = "-"
        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("x",))

    def test_left_extension(self):

        """
        R   |xxxxxxx|-------|xxxxx|
        P     |xxxxx|-------|xxxx|-------|xxxx|
        
        Expected ccode: j
        
        """

        reference = loci.Transcript()
        reference.start = 100
        reference.end = 2000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(100, 300), (500, 1000), (1500, 2000)]
        reference.finalize()

        prediction = loci.Transcript()
        prediction.start = 50
        prediction.end = 1800
        prediction.strand = "+"
        prediction.chrom = "Chr1"
        prediction.id = "P1.1"
        prediction.parent = "P1"
        prediction.exons = [(50, 150), (200, 300), (500, 1000), (1500, 1800)]
        prediction.finalize()

        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("J",))
        self.assertEqual(result.j_f1, (4 / 5 * 100,))
        self.assertAlmostEqual(result.j_prec[0], 2 / 3 * 100, delta=0.1)
        self.assertEqual(result.j_recall, (100,))

    def test_left_extension_n(self):

        """
        R   |xxxxxxx|-------|xxxxx|
        P     |xxxxx|-------|xxxxx|-------|xxxx|

        Expected ccode: j

        """

        reference = loci.Transcript()
        reference.start = 100
        reference.end = 2000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(100, 300), (500, 1000), (1500, 2000)]
        reference.finalize()

        prediction = loci.Transcript()
        prediction.start = 200
        prediction.end = 3000
        prediction.strand = "+"
        prediction.chrom = "Chr1"
        prediction.id = "P1.1"
        prediction.parent = "P1"
        prediction.exons = [(200, 300), (500, 1000), (1500, 2050), (2200, 3000)]
        prediction.finalize()

        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("n",))
        self.assertEqual(result.j_f1, (4 / 5 * 100,))
        self.assertAlmostEqual(result.j_prec[0], 2 / 3 * 100, delta=0.1)
        self.assertEqual(result.j_recall, (100,))

    def test_right_extension(self):
        """
        R                |xxxxxxx|-------|xxxxx|
        P     |xxxxx|-------|xxxx|-------|xxxx|
        
        Expected ccode: j
        """

        reference = loci.Transcript()
        reference.start = 100
        reference.end = 2000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(100, 300), (500, 1000), (1500, 2000)]
        reference.finalize()

        prediction = loci.Transcript()
        prediction.start = 200
        prediction.end = 3000
        prediction.strand = "+"
        prediction.chrom = "Chr1"
        prediction.id = "P1.1"
        prediction.parent = "P1"
        prediction.exons = [(200, 300), (500, 1000), (1500, 1800), (2500, 3000)]
        prediction.finalize()

        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("J",))
        self.assertEqual(result.j_f1, (4 / 5 * 100,))
        self.assertAlmostEqual(result.j_prec[0], 2 / 3 * 100, delta=0.1)
        self.assertEqual(result.j_recall, (100,))

    def test_left_right_extension(self):
        """
        R                |xxxxxxx|-------|xxxxx|
        P     |xxxxx|-------|xxxx|-------|xxxx|------|xxxxx|
        
        Expected ccode: j
        """

        reference = loci.Transcript()
        reference.start = 100
        reference.end = 2000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(100, 300), (500, 1000), (1500, 2000)]
        reference.finalize()

        prediction = loci.Transcript()
        prediction.start = 50
        prediction.end = 3000
        prediction.strand = "+"
        prediction.chrom = "Chr1"
        prediction.id = "P1.1"
        prediction.parent = "P1"
        prediction.exons = [(50, 150), (200, 300), (500, 1000), (1500, 1800), (2500, 3000)]
        prediction.finalize()

        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("J",))
        self.assertAlmostEqual(result.j_f1[0], 2 / 3 * 100, delta=0.1)
        self.assertAlmostEqual(result.j_prec[0], 1 / 2 * 100, delta=0.1)
        self.assertEqual(result.j_recall, (100,))

    def test_left_right_extension_novel(self):
        """
        R                   |xxxx|-------|xxxx|
        P     |xxxxx|----|xxxxxxx|-------|xxxxx|------|xxxxx|

        Expected ccode: n
        """

        reference = loci.Transcript()
        reference.start = 1000
        reference.end = 3000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(1000, 1300), (1500, 2000), (2500, 3000)]
        reference.finalize()

        prediction = loci.Transcript()
        prediction.start = 500
        prediction.end = 3500
        prediction.strand = "+"
        prediction.chrom = "Chr1"
        prediction.id = "P1.1"
        prediction.parent = "P1"
        prediction.exons = [(500, 700), (900, 1300), (1500, 2000), (2500, 3100), (3200, 3500)]
        prediction.finalize()

        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("n",))
        self.assertAlmostEqual(result.j_f1[0], 2 / 3 * 100, delta=0.1)
        self.assertAlmostEqual(result.j_prec[0], 1 / 2 * 100, delta=0.1)
        self.assertEqual(result.j_recall, (100,))

    def test_internal_extension(self):
        """
        
        R    |xxxxxxxx|-----------------|xxxxxxxx|
        P       |xxxxx|----|xxxxxx|-----|xxxxx|
        
        Expected ccode: j, junction recall: 100%
        """

        reference = loci.Transcript()
        reference.start = 100
        reference.end = 2000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(100, 300), (500, 1000), (1500, 2000)]
        reference.finalize()

        prediction = loci.Transcript()
        prediction.start = 200
        prediction.end = 2100
        prediction.strand = "+"
        prediction.chrom = "Chr1"
        prediction.id = "P1.1"
        prediction.parent = "P1"
        prediction.exons = [(200, 300), (500, 1000), (1200, 1300), (1500, 2100)]
        prediction.finalize()

        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("j",))
        self.assertAlmostEqual(result.j_f1[0], 4 / 5 * 100, delta=0.1)
        self.assertAlmostEqual(result.j_prec[0], 2 / 3 * 100, delta=0.1)
        self.assertEqual(result.j_recall, (100,))

    def test_internal_external_extension(self):
        """
        
        R      |xxxxxx|-----------------|xxxxxxxx|
        P       |xxxxx|----|xxxxxx|-----|xxxxx|------|xxxxxxx|
        
        Expected ccode: j, junction recall: 100%
        """

        reference = loci.Transcript()
        reference.start = 100
        reference.end = 2000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(100, 300), (500, 1000), (1500, 2000)]
        reference.finalize()

        prediction = loci.Transcript()
        prediction.start = 200
        prediction.end = 3000
        prediction.strand = "+"
        prediction.chrom = "Chr1"
        prediction.id = "P1.1"
        prediction.parent = "P1"
        prediction.exons = [(200, 300), (500, 1000), (1200, 1300), (1500, 1800), (2500, 3000)]
        prediction.finalize()

        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("j",))
        self.assertAlmostEqual(result.j_f1[0], 2 / 3 * 100, delta=0.1)
        self.assertAlmostEqual(result.j_prec[0], 1 / 2 * 100, delta=0.1)
        self.assertEqual(result.j_recall, (100,))

    def test_left_right_internal_extension_novel(self):
        """
        R                   |xxxx|-------|xxxx|
        P     |xxxxx|----|xxxxxxx|--|x|--|xxxxx|------|xxxxx|

        Expected ccode: j
        """

        reference = loci.Transcript()
        reference.start = 1000
        reference.end = 3000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(1000, 1300), (1500, 2000), (2500, 3000)]
        reference.finalize()

        prediction = loci.Transcript()
        prediction.start = 500
        prediction.end = 3500
        prediction.strand = "+"
        prediction.chrom = "Chr1"
        prediction.id = "P1.1"
        prediction.parent = "P1"
        prediction.exons = [(500, 700), (900, 1300), (1500, 2000),
                            (2200,2300),
                            (2500, 3100), (3200, 3500)]
        prediction.finalize()

        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("j",))
        self.assertAlmostEqual(result.j_f1[0], 2 * 0.4 / 1.4 * 100, delta=0.1)
        self.assertAlmostEqual(result.j_prec[0], 2 / 5 * 100, delta=0.1)
        self.assertEqual(result.j_recall, (100,))

    def test_contained_bleeding(self):

        """
        R     |xxxxxxx|--------|xxxxxx|----|xxxx|--------|xxxxxxx|
        P                  |xxxxxxxxxx|----|xxxxxx|

        Expected class code: C
        :return:
        """
        reference = loci.Transcript()
        reference.start = 1000
        reference.end = 3000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(1000, 1300), (1500, 2000), (2200, 2400), (2500, 3000)]
        reference.finalize()

        prediction = loci.Transcript()
        prediction.start = 1400
        prediction.end = 2450
        prediction.strand = "+"
        prediction.chrom = "Chr1"
        prediction.id = "P1.1"
        prediction.parent = "P1"
        prediction.exons = [(1400, 2000), (2200, 2450)]
        prediction.finalize()

        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("C",), result)

    def test_contained_alternative(self):

        """
        R     |xxxxxxx|--------|xxxxxx|----|xxxx|--------|xxxxxxx|
        P           |xxxxxxxxxxxxxxxxx|----|xxxxxx|

        Expected class code: C
        :return:
        """
        reference = loci.Transcript()
        reference.start = 1000
        reference.end = 3000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(1000, 1300), (1500, 2000), (2200, 2400), (2500, 3000)]
        reference.finalize()

        prediction = loci.Transcript()
        prediction.start = 1200
        prediction.end = 2450
        prediction.strand = "+"
        prediction.chrom = "Chr1"
        prediction.id = "P1.1"
        prediction.parent = "P1"
        prediction.exons = [(1200, 2000), (2200, 2450)]
        prediction.finalize()

        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("j",))

    def test_mono_overlap_nostrand(self):

        reference = loci.Transcript()
        reference.start = 100
        reference.end = 2000
        reference.strand = "-"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(100, 2000)]
        reference.finalize()

        prediction = loci.Transcript()
        prediction.start = 300
        prediction.end = 3000
        prediction.strand = "."
        prediction.chrom = "Chr1"
        prediction.id = "G1.1"
        prediction.parent = "G1"
        prediction.exons = [(300, 3000)]
        prediction.finalize()

        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("m",))

        result, _ = Assigner.compare(reference, prediction)
        self.assertEqual(result.ccode, ("m",))

    def test_mono_multi_overlap_nostrand(self):

        """Test a monoexonic transcript against a multiexonic transcript
         to verify the assignment of correct ccodes (o and O)
        """

        reference = loci.Transcript()
        reference.start = 100
        reference.end = 2000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(100, 700), (1000, 2000)]
        reference.finalize()

        prediction = loci.Transcript()
        prediction.start = 50
        prediction.end = 600
        prediction.strand = "."
        prediction.chrom = "Chr1"
        prediction.id = "G1.1"
        prediction.parent = "G1"
        prediction.exons = [(50, 600)]
        prediction.finalize()

        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("g",))

        result, _ = Assigner.compare(reference, prediction)
        self.assertEqual(result.ccode, ("G",))

    def test_neighbors(self):

        """
        Test for the assignment of transcripts inside the index. Still a stub
        """
        keys = [(10, 200), (350, 500)]
        keys = utilities.intervaltree.IntervalTree.from_tuples(keys)
        self.assertEqual(Assigner.find_neighbours(
                         keys, (350,500), distance=2000),
                         [((350, 500), 0), ((10, 200), 150)],
                         )
        self.assertEqual(Assigner.find_neighbours(
                         keys, (5350,5500), distance=1000),
                         []
                         )

        self.assertEqual(Assigner.find_neighbours(
                         keys, (5350,5500), distance=10000),
                         [((350, 500), 4850), ((10, 200), 5150)]
                         )

    def test_false_fusion(self):

        """
        System test to verify that the false fusion is not called.
        WARNING: this test is quite brittle as I am creating the namespace myself
        instead of loading it from the subprograms.compare module.
        :return:
        """

        master = os.path.dirname(os.path.abspath(__file__))
        out = tempfile.mkdtemp()
        os.makedirs(out, exist_ok=True)

        args = argparse.Namespace()
        args.no_save_index = True
        args.reference = parsers.GFF.GFF3(
            os.path.join(master, "fusion_test", "fusion_test_ref.gff3"))
        args.prediction = parsers.GTF.GTF(
            os.path.join(master, "fusion_test", "fusion_test_pred.gtf"))
        args.log = os.path.join(out, "fusion_test", "fusion_test.log")
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

        scales.compare.compare(args)

        out_refmap = os.path.join(out, "fusion_test", "fusion_test.refmap")
        self.assertTrue(os.path.exists(out_refmap))
        self.assertGreater(os.stat(out_refmap).st_size, 0)
        with open(out_refmap) as refmap:
            for line in csv.DictReader(refmap, delimiter="\t"):
                if line["ref_id"] not in ("AT1G78880.1", "AT1G78882.1"):
                    continue
                self.assertEqual(line["ccode"], "=", line)
        args.reference.close()
        args.prediction.close()
        import shutil
        try:
            shutil.rmtree(out)
        except PermissionError:
            raise PermissionError(out)

    def test_monoexonic_contained(self):
        t1 = loci.Transcript()
        t1.chrom, t1.start, t1.end, t1.strand, t1.id = "Chr1", 9192977, 9193471, "+", "t1"
        t1.add_exons([(t1.start, t1.end)])
        t1.add_exons([(t1.start, t1.end)], features=["CDS"])
        t1.finalize()
        t2 = loci.Transcript()
        t2.chrom, t2.start, t2.end, t2.strand, t2.id = "Chr1", 9182621, 9194526, "+", "t2"
        t2.add_exons([(9192977, 9193480)], features=["CDS"])
        t2.add_exons([(9186643, 9186780), (9182621, 9182770), (9192959, 9194526)])
        t2.finalize()
        mono_ref = Assigner.compare(reference=t1, prediction=t2)[0]
        mono_pred = Assigner.compare(reference=t2, prediction=t1)[0]
        self.assertEqual(mono_ref.n_recall, mono_pred.n_prec)
        self.assertEqual(mono_pred.ccode[0], "c")
        self.assertEqual(mono_ref.ccode[0], "n")

    def test_h_case(self):

        """
        |xxxxxxxxxxxxxxx|-------------------|xxxxxxxxxxxx|
        |xxxxxxxxxxxxxxxxxxx|-----------------|xxxxxxxxx|
        :return:
        """

        reference = loci.Transcript()
        reference.start = 1000
        reference.end = 3000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(1000, 1300), (2500, 3000)]
        reference.finalize()

        prediction = loci.Transcript()
        prediction.start = 1200
        prediction.end = 2450
        prediction.strand = "+"
        prediction.chrom = "Chr1"
        prediction.id = "P1.1"
        prediction.parent = "P1"
        prediction.exons = [(1200, 2000), (2200, 2450)]
        prediction.finalize()

        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("h",))

    def test_double_h_case(self):

        """
        |xxxxxxxxxxxxxxx|---------|xxxxxxx|-------|xxxxxxxxxxxx|
        |xxxxxxxxxxxxxxxxxxx|-----------------|xxxxxxxxx|
        :return:
        """

        reference = loci.Transcript()
        reference.start = 1000
        reference.end = 3000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(1000, 1300), (1700, 2000), (2500, 3000)]
        reference.finalize()

        prediction = loci.Transcript()
        prediction.start = 1200
        prediction.end = 2450
        prediction.strand = "+"
        prediction.chrom = "Chr1"
        prediction.id = "P1.1"
        prediction.parent = "P1"
        prediction.exons = [(1200, 1500), (2200, 2450)]
        prediction.finalize()

    def test_non_h_case(self):

        """
        |xxxxxxxxxxxxxxx|-------------------|xxxxxxxxxxxx|
        |xxxxx|-----|xxxxxxxxx|
        :return:
        """

        reference = loci.Transcript()
        reference.start = 1000
        reference.end = 3000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(1000, 1800), (2500, 3000)]
        reference.finalize()

        prediction = loci.Transcript()
        prediction.start = 1200
        prediction.end = 2050
        prediction.strand = "+"
        prediction.chrom = "Chr1"
        prediction.id = "P1.1"
        prediction.parent = "P1"
        prediction.exons = [(1200, 1600), (1700, 2050)]
        prediction.finalize()

        result, _ = Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("o",))

    def test_J_and_C_case(self):

        """

        1  xxxxxxxxx------------xxxxxxx-----------xxxxxxxx-------------xxxxxxxxxxxxx
        2                 xxxxxxxxxxxxx-----------xxxxxxxxxxxx

        We do expect the comparison to be:

        1 reference, 2 prediction: C
        1 prediction, 2 reference: J

        :return:
        """

        t1 = loci.Transcript()
        t1.chrom = "Chr1"
        t1.strand = "+"
        t1.score = 10
        t1.source = "mikado"
        t1.id = "t1.1"
        t1.parent = "t1"
        t1.start = 101
        t1.end = 5000
        t1.add_exons([(101, 1000), (1501, 2000), (2301, 4000), (4501, 5000)])
        t1.finalize()

        t2 = loci.Transcript()
        for attr in ["chrom", "source", "strand", "score"]:
            setattr(t2, attr, getattr(t1, attr))
        t2.id = "t2.1"
        t2.parent = "t2"
        t2.start = 1300
        t2.end = 4300
        t2.add_exons([(1300,2000), (2301, 4300)])
        t2.finalize()

        result, _ = Assigner.compare(t2, t1)
        self.assertEqual(result.ccode, ("C",))

        result, _ = Assigner.compare(t1, t2)
        self.assertEqual(result.ccode, ("J",))

    def test_J_and_C_case_in_exon(self):

        """

        1  xxxxxxxxx------------xxxxxxx-----------xxxxxxxx-------------xxxxxxxxxxxxx
        2         xxxxxxxxxxxxxxxxxxxxx-----------xxxxxxxxxxxx

        We do expect the comparison to be:

        1 reference, 2 prediction: j
        1 prediction, 2 reference: J

        Notice the class code switch from C to j due to the first exon being internal to the terminal exon,
        rather than internal to the intron.

        :return:
        """

        t1 = loci.Transcript()
        t1.chrom = "Chr1"
        t1.strand = "+"
        t1.score = 10
        t1.source = "mikado"
        t1.id = "t1.1"
        t1.parent = "t1"
        t1.start = 101
        t1.end = 5000
        t1.add_exons([(101, 1000), (1501, 2000), (2301, 4000), (4501, 5000)])
        t1.finalize()

        t2 = loci.Transcript()
        for attr in ["chrom", "source", "strand", "score"]:
            setattr(t2, attr, getattr(t1, attr))
        t2.id = "t2.1"
        t2.parent = "t2"
        t2.start = 900
        t2.end = 4300
        t2.add_exons([(900,2000), (2301, 4300)])
        t2.finalize()

        result, _ = Assigner.compare(t2, t1)
        self.assertEqual(result.ccode, ("j",))

        result, _ = Assigner.compare(t1, t2)
        self.assertEqual(result.ccode, ("J",))

    def test_fuzzy_match(self):

        t1 = loci.Transcript()
        t1.chrom = "Chr1"
        t1.strand = "+"
        t1.score = 10
        t1.source = "mikado"
        t1.id = "t1.1"
        t1.parent = "t1"
        t1.start = 101
        t1.end = 5000
        t1.add_exons([(101, 1000), (1501, 2000), (2301, 4000), (4501, 5000)])
        t1.finalize()

        t2 = loci.Transcript()
        t2.chrom = "Chr1"
        t2.strand = "+"
        t2.score = 10
        t2.source = "mikado"
        t2.id = "t2.1"
        t2.parent = "t2"
        t2.start = 401
        t2.end = 5130
        t2.add_exons([(401, 1010), (1491, 2005), (2305, 4002), (4509, 5130)])
        t2.finalize()

        for fuzzymatch in (0, 1, 5, 10, 20, 30):
            with self.subTest(fuzzymatch=fuzzymatch):
                result, _ = Assigner.compare(t2, t1, fuzzy_match=fuzzymatch)
                if fuzzymatch < 5:
                    self.assertEqual(result.ccode, ("h",))
                elif fuzzymatch >= 10:
                    self.assertEqual(result.ccode, ("=",))
                else:
                    self.assertEqual(result.ccode, ("j",), fuzzymatch)


if __name__ == '__main__':
    unittest.main()
