#!/usr/bin/env python3
# coding: utf-8

"""
Unit tests for the scales library
"""

import unittest
import mikado_lib.loci_objects
import mikado_lib.scales
import intervaltree

class AssignerTester(unittest.TestCase):
    """
    This unit test has the purpose of testing the scales module of mikado_lib.
    """

    def test_self(self):
        reference = mikado_lib.loci_objects.transcript.Transcript()
        reference.start = 100
        reference.end = 2000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(100, 300), (500, 1000), (1500, 2000)]
        reference.exons = [intervaltree.Interval(*_) for _ in reference.exons]
        reference.finalize()

        result, _ = mikado_lib.scales.assigner.Assigner.compare(reference, reference)
        self.assertEqual(result.ccode, ("=",))
        self.assertEqual(result.n_f1, (100,))
        self.assertEqual(result.n_prec, (100,))
        self.assertEqual(result.n_recall, (100,))
        self.assertEqual(result.j_f1, (100,))
        self.assertEqual(result.j_prec, (100,))
        self.assertEqual(result.j_recall, (100,))

    def test_mono_self(self):
        reference = mikado_lib.loci_objects.transcript.Transcript()
        reference.start = 100
        reference.end = 1000
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.strand = "+"
        reference.exons = [(100, 1000)]
        reference.finalize()

        result, _ = mikado_lib.scales.assigner.Assigner.compare(reference, reference)
        self.assertEqual(result.ccode, ("_",))
        self.assertEqual(result.n_f1, (100,))
        self.assertEqual(result.n_prec, (100,))
        self.assertEqual(result.n_recall, (100,))
        self.assertEqual(result.j_f1, (100,))
        self.assertEqual(result.j_prec, (100,))
        self.assertEqual(result.j_recall, (100,))

    def test_equal(self):
        reference = mikado_lib.loci_objects.transcript.Transcript()
        reference.start = 100
        reference.end = 2000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(100, 300), (500, 1000), (1500, 2000)]
        reference.finalize()

        prediction = mikado_lib.loci_objects.transcript.Transcript()
        prediction.start = 200
        prediction.end = 1800
        prediction.strand = "+"
        prediction.chrom = "Chr1"
        prediction.id = "P1.1"
        prediction.parent = "P1"
        prediction.exons = [(200, 300), (500, 1000), (1500, 1800)]
        prediction.finalize()

        result, _ = mikado_lib.scales.assigner.Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("=",))
        self.assertEqual(result.j_f1, (100,))
        self.assertEqual(result.j_prec, (100,))
        self.assertEqual(result.j_recall, (100,))
        self.assertAlmostEqual(result.n_prec[0], 100, delta=0.1)
        self.assertAlmostEqual(result.n_recall[0],
                               100 * (300 - 200 + 1 + 1000 - 500 + 1 + 1800 - 1500 + 1) / reference.cdna_length,
                               delta=0.1)

        prediction.strand = "-"
        result, _ = mikado_lib.scales.assigner.Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("x",))

    def test_mono_equal(self):
        reference = mikado_lib.loci_objects.transcript.Transcript()
        reference.start = 100
        reference.end = 1000
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.strand = "+"
        reference.exons = [(100, 1000)]
        reference.finalize()

        prediction = mikado_lib.loci_objects.transcript.Transcript()
        prediction.start = 105
        prediction.end = 995
        prediction.chrom = "Chr1"
        prediction.id = "G1.1"
        prediction.parent = "G1"
        prediction.strand = "+"
        prediction.exons = [(105, 995)]
        prediction.finalize()

        result, _ = mikado_lib.scales.assigner.Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("_",))
        self.assertEqual(result.j_f1, (100,))
        self.assertEqual(result.j_prec, (100,))
        self.assertEqual(result.j_recall, (100,))
        self.assertAlmostEqual(result.n_prec[0], 100, delta=0.1)
        self.assertAlmostEqual(result.n_recall[0], 100 * (995 - 105 + 1) / reference.cdna_length, delta=0.1)

        prediction.strand = "-"
        result, _ = mikado_lib.scales.assigner.Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("x",))

    def test_mono_semiequal(self):
        reference = mikado_lib.loci_objects.transcript.Transcript()
        reference.start = 100
        reference.end = 1000
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.strand = "+"
        reference.exons = [(100, 1000)]
        reference.finalize()

        prediction = mikado_lib.loci_objects.transcript.Transcript()
        prediction.start = 200
        prediction.end = 900
        prediction.chrom = "Chr1"
        prediction.id = "G1.1"
        prediction.parent = "G1"
        prediction.strand = "+"
        prediction.exons = [(200, 900)]
        prediction.finalize()

        result, _ = mikado_lib.scales.assigner.Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("c",))
        self.assertEqual(result.j_f1, (100,))
        self.assertEqual(result.j_prec, (100,))
        self.assertEqual(result.j_recall, (100,))
        self.assertAlmostEqual(result.n_prec[0], 100, delta=0.1)
        self.assertAlmostEqual(result.n_recall[0], 100 * (900 - 200 + 1) / reference.cdna_length,
                               delta=0.1)

        prediction.strand = "-"
        result, _ = mikado_lib.scales.assigner.Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("x",))

    def test_mono_overlap(self):

        """Test that two monoexonic overlapping genes get a m"""

        reference = mikado_lib.loci_objects.transcript.Transcript()
        reference.start = 100
        reference.end = 1000
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.strand = "+"
        reference.exons = [(100, 1000)]
        reference.finalize()

        prediction = mikado_lib.loci_objects.transcript.Transcript()
        prediction.start = 200
        prediction.end = 1200
        prediction.chrom = "Chr1"
        prediction.id = "G1.1"
        prediction.parent = "G1"
        prediction.strand = "+"
        prediction.exons = [(200, 1200)]
        prediction.finalize()

        result, _ = mikado_lib.scales.assigner.Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("m",))
        self.assertEqual(result.j_f1, (100,))
        self.assertEqual(result.j_prec, (100,))
        self.assertEqual(result.j_recall, (100,))
        self.assertAlmostEqual(result.n_prec[0], 100 * 801 / prediction.cdna_length, delta=0.1)
        self.assertAlmostEqual(result.n_recall[0], 100 * 801 / reference.cdna_length, delta=0.1)

    def test_contained(self):

        reference = mikado_lib.loci_objects.transcript.Transcript()
        reference.start = 100
        reference.end = 2000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(100, 300), (500, 1000), (1500, 2000)]
        reference.finalize()

        prediction = mikado_lib.loci_objects.transcript.Transcript()
        prediction.start = 600
        prediction.end = 1800
        prediction.strand = "+"
        prediction.chrom = "Chr1"
        prediction.id = "P1.1"
        prediction.parent = "P1"
        prediction.exons = [(600, 1000), (1500, 1800)]
        prediction.finalize()

        result, _ = mikado_lib.scales.assigner.Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("c",))
        self.assertAlmostEqual(result.j_f1[0], (100 * (2 / 3),)[0], delta=0.1)
        self.assertAlmostEqual(result.j_prec, (100,))
        self.assertAlmostEqual(result.j_recall, (100 * 1 / 2,), delta=0.1)
        self.assertAlmostEqual(result.n_prec[0], 100, delta=0.1)

        prediction.strand = "-"
        result, _ = mikado_lib.scales.assigner.Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("x",))

    def test_alternative(self):

        reference = mikado_lib.loci_objects.transcript.Transcript()
        reference.start = 100
        reference.end = 2000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(100, 300), (500, 1000), (1500, 2000)]
        reference.finalize()

        prediction = mikado_lib.loci_objects.transcript.Transcript()
        prediction.start = 100
        prediction.end = 2000
        prediction.strand = "+"
        prediction.chrom = "Chr1"
        prediction.id = "P1.1"
        prediction.parent = "P1"
        prediction.exons = [(100, 300), (600, 1000), (1500, 2000)]
        prediction.finalize()

        result, _ = mikado_lib.scales.assigner.Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("j",))
        self.assertAlmostEqual(result.j_f1[0], (100 * ((3 / 4) / 1),)[0], delta=0.1)
        self.assertAlmostEqual(result.j_prec, (100 * 3 / 4,), delta=0.1)
        self.assertAlmostEqual(result.j_recall, (100 * 3 / 4,), delta=0.1)

        prediction.strand = "-"
        result, _ = mikado_lib.scales.assigner.Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("j",))

    def test_mono_intronic(self):

        """
        R   |=====|-------------------|=========|
        P            |=======|
        
        Expected class code: i
        """

        reference = mikado_lib.loci_objects.transcript.Transcript()
        reference.start = 100
        reference.end = 2000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(100, 300), (500, 1000), (1500, 2000)]
        reference.finalize()

        prediction = mikado_lib.loci_objects.transcript.Transcript()
        prediction.start = 1100
        prediction.end = 1400
        prediction.strand = "+"
        prediction.chrom = "Chr1"
        prediction.id = "P1.1"
        prediction.parent = "P1"
        prediction.exons = [(1100, 1400)]
        prediction.finalize()

        result, _ = mikado_lib.scales.assigner.Assigner.compare(prediction, reference)
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
        result, _ = mikado_lib.scales.assigner.Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("i",))

    def test_multi_intronic(self):

        """
        R   |=====|-------------------|=========|
        P            |===|----|====|
        
        OR
        
        R   |=====|----------|=========|---------|=========|
        P            |===|----------------|====|
        
        Expected class code: I
        """

        reference = mikado_lib.loci_objects.transcript.Transcript()
        reference.start = 100
        reference.end = 2000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(100, 300), (500, 1000), (1500, 2000)]
        reference.finalize()

        prediction = mikado_lib.loci_objects.transcript.Transcript()
        prediction.start = 1100
        prediction.end = 1400
        prediction.strand = "+"
        prediction.chrom = "Chr1"
        prediction.id = "P1.1"
        prediction.parent = "P1"
        prediction.exons = [(1100, 1200), (1300, 1400)]
        prediction.finalize()

        result, _ = mikado_lib.scales.assigner.Assigner.compare(prediction, reference)
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
        result, _ = mikado_lib.scales.assigner.Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("I",))

        # Now the reference spans two introns
        prediction = mikado_lib.loci_objects.transcript.Transcript()
        prediction.start = 350
        prediction.end = 1400
        prediction.strand = "+"
        prediction.chrom = "Chr1"
        prediction.id = "P1.1"
        prediction.parent = "P1"
        prediction.exons = [(350, 450), (1300, 1400)]
        prediction.finalize()

        result, _ = mikado_lib.scales.assigner.Assigner.compare(prediction, reference)
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
        result, _ = mikado_lib.scales.assigner.Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("I",))

    def test_overlap(self):

        """
        R   |=====|-------|=====|
        P      |=====|---|====|
        
        No junction in common
        
        """

        reference = mikado_lib.loci_objects.transcript.Transcript()
        reference.start = 100
        reference.end = 2000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(100, 300), (500, 1000), (1500, 2000)]
        reference.finalize()

        prediction = mikado_lib.loci_objects.transcript.Transcript()
        prediction.start = 200
        prediction.end = 1300
        prediction.strand = "+"
        prediction.chrom = "Chr1"
        prediction.id = "P1.1"
        prediction.parent = "P1"
        prediction.exons = [(200, 700), (900, 1300)]
        prediction.finalize()

        result, _ = mikado_lib.scales.assigner.Assigner.compare(prediction, reference)
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
        result, _ = mikado_lib.scales.assigner.Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("x",))

    def test_ccode_e(self):

        """Case:
        
        R ---|=====|-------|======|---
        P      |======|
        
        Exonic and intronic overlap
        
        """

        reference = mikado_lib.loci_objects.transcript.Transcript()
        reference.start = 100
        reference.end = 2000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(100, 300), (500, 1000), (1500, 2000)]
        reference.finalize()

        prediction = mikado_lib.loci_objects.transcript.Transcript()
        prediction.start = 200
        prediction.end = 400
        prediction.strand = "+"
        prediction.chrom = "Chr1"
        prediction.id = "P1.1"
        prediction.parent = "P1"
        prediction.exons = [(200, 400)]
        prediction.finalize()

        result, _ = mikado_lib.scales.assigner.Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("e",))
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
        result, _ = mikado_lib.scales.assigner.Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("x",))

    def test_not_ccode_e(self):

        """Case:
        
        R ---|=====|-------|======|---
        P |======|
        
        Exonic overlap only
        
        """
        reference = mikado_lib.loci_objects.transcript.Transcript()
        reference.start = 100
        reference.end = 2000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(100, 300), (500, 1000), (1500, 2000)]
        reference.finalize()

        prediction = mikado_lib.loci_objects.transcript.Transcript()
        prediction.start = 50
        prediction.end = 310
        prediction.strand = "+"
        prediction.chrom = "Chr1"
        prediction.id = "P1.1"
        prediction.parent = "P1"
        prediction.exons = [(50, 310)]
        prediction.finalize()

        result, _ = mikado_lib.scales.assigner.Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("o",))
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
        result, _ = mikado_lib.scales.assigner.Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("x",))

    def test_left_extension(self):

        """
        R   |=======|-------|=====|
        P     |=====|-------|====|-------|====|
        
        Expected ccode: n
        
        """

        reference = mikado_lib.loci_objects.transcript.Transcript()
        reference.start = 100
        reference.end = 2000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(100, 300), (500, 1000), (1500, 2000)]
        reference.finalize()

        prediction = mikado_lib.loci_objects.transcript.Transcript()
        prediction.start = 50
        prediction.end = 1800
        prediction.strand = "+"
        prediction.chrom = "Chr1"
        prediction.id = "P1.1"
        prediction.parent = "P1"
        prediction.exons = [(50, 150), (200, 300), (500, 1000), (1500, 1800)]
        prediction.finalize()

        result, _ = mikado_lib.scales.assigner.Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("n",))
        self.assertEqual(result.j_f1, (4 / 5 * 100,))
        self.assertAlmostEqual(result.j_prec[0], 2 / 3 * 100, delta=0.1)
        self.assertEqual(result.j_recall, (100,))

    def test_right_extension(self):
        """
        R                |=======|-------|=====|
        P     |=====|-------|====|-------|====|
        
        Expected ccode: n
        """

        reference = mikado_lib.loci_objects.transcript.Transcript()
        reference.start = 100
        reference.end = 2000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(100, 300), (500, 1000), (1500, 2000)]
        reference.finalize()

        prediction = mikado_lib.loci_objects.transcript.Transcript()
        prediction.start = 200
        prediction.end = 3000
        prediction.strand = "+"
        prediction.chrom = "Chr1"
        prediction.id = "P1.1"
        prediction.parent = "P1"
        prediction.exons = [(200, 300), (500, 1000), (1500, 1800), (2500, 3000)]
        prediction.finalize()

        result, _ = mikado_lib.scales.assigner.Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("n",))
        self.assertEqual(result.j_f1, (4 / 5 * 100,))
        self.assertAlmostEqual(result.j_prec[0], 2 / 3 * 100, delta=0.1)
        self.assertEqual(result.j_recall, (100,))

    def test_left_right_extension(self):
        """
        R                |=======|-------|=====|
        P     |=====|-------|====|-------|====|------|=====|
        
        Expected ccode: n
        """

        reference = mikado_lib.loci_objects.transcript.Transcript()
        reference.start = 100
        reference.end = 2000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(100, 300), (500, 1000), (1500, 2000)]
        reference.finalize()

        prediction = mikado_lib.loci_objects.transcript.Transcript()
        prediction.start = 50
        prediction.end = 3000
        prediction.strand = "+"
        prediction.chrom = "Chr1"
        prediction.id = "P1.1"
        prediction.parent = "P1"
        prediction.exons = [(50, 150), (200, 300), (500, 1000), (1500, 1800), (2500, 3000)]
        prediction.finalize()

        result, _ = mikado_lib.scales.assigner.Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("n",))
        self.assertAlmostEqual(result.j_f1[0], 2 / 3 * 100, delta=0.1)
        self.assertAlmostEqual(result.j_prec[0], 1 / 2 * 100, delta=0.1)
        self.assertEqual(result.j_recall, (100,))

    def test_internal_extension(self):
        """
        
        R    |========|-----------------|========|
        P       |=====|----|======|-----|=====|
        
        Expected ccode: j, junction recall: 100%
        """

        reference = mikado_lib.loci_objects.transcript.Transcript()
        reference.start = 100
        reference.end = 2000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(100, 300), (500, 1000), (1500, 2000)]
        reference.finalize()

        prediction = mikado_lib.loci_objects.transcript.Transcript()
        prediction.start = 200
        prediction.end = 2100
        prediction.strand = "+"
        prediction.chrom = "Chr1"
        prediction.id = "P1.1"
        prediction.parent = "P1"
        prediction.exons = [(200, 300), (500, 1000), (1200, 1300), (1500, 2100)]
        prediction.finalize()

        result, _ = mikado_lib.scales.assigner.Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("j",))
        self.assertAlmostEqual(result.j_f1[0], 4 / 5 * 100, delta=0.1)
        self.assertAlmostEqual(result.j_prec[0], 2 / 3 * 100, delta=0.1)
        self.assertEqual(result.j_recall, (100,))

    def test_internal_external_extension(self):
        """
        
        R      |======|-----------------|========|
        P       |=====|----|======|-----|=====|------|=======|
        
        Expected ccode: j, junction recall: 100%
        """

        reference = mikado_lib.loci_objects.transcript.Transcript()
        reference.start = 100
        reference.end = 2000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(100, 300), (500, 1000), (1500, 2000)]
        reference.finalize()

        prediction = mikado_lib.loci_objects.transcript.Transcript()
        prediction.start = 200
        prediction.end = 3000
        prediction.strand = "+"
        prediction.chrom = "Chr1"
        prediction.id = "P1.1"
        prediction.parent = "P1"
        prediction.exons = [(200, 300), (500, 1000), (1200, 1300), (1500, 1800), (2500, 3000)]
        prediction.finalize()

        result, _ = mikado_lib.scales.assigner.Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("j",))
        self.assertAlmostEqual(result.j_f1[0], 2 / 3 * 100, delta=0.1)
        self.assertAlmostEqual(result.j_prec[0], 1 / 2 * 100, delta=0.1)
        self.assertEqual(result.j_recall, (100,))

    def test_mono_overlap_nostrand(self):

        reference = mikado_lib.loci_objects.transcript.Transcript()
        reference.start = 100
        reference.end = 2000
        reference.strand = "-"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(100, 2000)]
        reference.finalize()

        prediction = mikado_lib.loci_objects.transcript.Transcript()
        prediction.start = 300
        prediction.end = 2300
        prediction.strand = "."
        prediction.chrom = "Chr1"
        prediction.id = "G1.1"
        prediction.parent = "G1"
        prediction.exons = [(300, 2300)]
        prediction.finalize()

        result, _ = mikado_lib.scales.assigner.Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("m",))

        result, _ = mikado_lib.scales.assigner.Assigner.compare(reference, prediction)
        self.assertEqual(result.ccode, ("m",))

    def test_mono_multi_overlap_nostrand(self):

        """Test a monoexonic transcript against a multiexonic transcript
         to verify the assignment of correct ccodes (o and O)
        """

        reference = mikado_lib.loci_objects.transcript.Transcript()
        reference.start = 100
        reference.end = 2000
        reference.strand = "+"
        reference.chrom = "Chr1"
        reference.id = "G1.1"
        reference.parent = "G1"
        reference.exons = [(100, 700), (1000, 2000)]
        reference.finalize()

        prediction = mikado_lib.loci_objects.transcript.Transcript()
        prediction.start = 50
        prediction.end = 600
        prediction.strand = "."
        prediction.chrom = "Chr1"
        prediction.id = "G1.1"
        prediction.parent = "G1"
        prediction.exons = [(50, 600)]
        prediction.finalize()

        result, _ = mikado_lib.scales.assigner.Assigner.compare(prediction, reference)
        self.assertEqual(result.ccode, ("o",))

        result, _ = mikado_lib.scales.assigner.Assigner.compare(reference, prediction)
        self.assertEqual(result.ccode, ("O",))

    def test_neighbors(self):

        """
        Test for the assignment of transcripts inside the index. Still a stub
        """
        keys = [(10, 200), (350, 500)]
        keys = intervaltree.IntervalTree.from_tuples(keys)
        self.assertEqual(mikado_lib.scales.assigner.Assigner.find_neighbours(
                         keys, (350,500)),
                         [((350, 500), 0), ((10, 200), 150)]
                         )
        self.assertEqual(mikado_lib.scales.assigner.Assigner.find_neighbours(
                         keys, (5350,5500), distance=1000),
                         []
                         )

        self.assertEqual(mikado_lib.scales.assigner.Assigner.find_neighbours(
                         keys, (5350,5500), distance=10000),
                         [((350, 500), 4850), ((10, 200), 5150)]
                         )

unittest.main()
