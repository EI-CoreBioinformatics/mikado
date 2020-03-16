#!/usr/bin/env python3


import tempfile
from ..parsers import blast_utils
import unittest
import os
import gzip
import subprocess
from ..serializers.blast_serializer import xml_utils as seri_blast_utils
from ..serializers.blast_serializer.tabular_utils import matrices
from ..serializers.blast_serializer.btop_parser import parse_btop
import numpy as np
import time
import itertools


class BtopTester(unittest.TestCase):

    def test_btop_equals(self):
        matrix=matrices["blosum62"]
        for qmult, tmult in itertools.product([1, 3], [1, 3]):
            for qsize, ssize in itertools.product(range(10 * qmult, 20 * qmult, 10),
                                                  range(10 * tmult, 20 * tmult, 10)):
                for qpos, spos in itertools.product(range(0, qsize, 2), range(0, ssize, 2)):
                    mqlength = (qsize - qpos) // qmult
                    mslength = (ssize - spos) // tmult
                    for mlength in range(1, int(min(mqlength, mslength))):
                        # with self.subTest(qsize=qsize, ssize=ssize, qpos=qpos, spos=spos, mlength=mlength):

                        qar = np.zeros([3, qsize], dtype=np.int)
                        sar = np.zeros([3, ssize], dtype=np.int)
                        self.assertEqual(qar.max(), 0)
                        self.assertEqual(sar.max(), 0)
                        sm = str(mlength)
                        #str btop, long qpos, long spos,
                        # np.ndarray[dtype=np.int, ndim=2, cast=True] query_array,
                        # np.ndarray[dtype=np.int, ndim=2, cast=True] target_array,
                        # dict matrix, long qmult=3, long tmult=1):
                        qar, sar, tot = parse_btop(sm, qpos, spos, qar, sar, matrix, qmult=qmult, tmult=tmult)
                        qfound = np.where(qar > 0)
                        sfound = np.where(sar > 0)
                        self.assertEqual(tot, mlength * min(qmult, tmult),
                                         (tot, mlength, qpos, spos, qmult, tmult))
                        self.assertTrue(
                            ((qfound[1][qfound[0] == 0] == qfound[1][qfound[0] == 1]) &
                             (qfound[1][qfound[0] == 0] == qfound[1][qfound[0] == 2])).all()
                        )
                        self.assertTrue(
                            ((sfound[1][sfound[0] == 0] == sfound[1][sfound[0] == 1]) &
                             (sfound[1][sfound[0] == 0] == sfound[1][sfound[0] == 2])).all()
                        )
                        self.assertEqual(qfound[1][qfound[0] == 0].shape[0], mlength * qmult)
                        self.assertEqual(sfound[1][sfound[0] == 0].shape[0], mlength * tmult,
                                         (tot, tmult, sar, sfound[1][sfound[0] == 0]))
                        # self.assertEqual(qfound[0].shape[0], 3)
                        self.assertEqual(np.where(qar[0] > 0)[0].min(), qpos,
                                         (qfound[0], qpos, mlength))
                        self.assertEqual(np.where(qar[0] > 0)[0].max(), qpos + mlength * qmult - 1,
                                         (qar[0], qpos + tot * qmult, tot))
                        self.assertEqual(np.where(sar[0] > 0)[0].min(), spos)
                        self.assertEqual(np.where(sar[0] > 0)[0].max(), spos + mlength * tmult - 1)


class BlastBasics(unittest.TestCase):

    def test_sniff_correct(self):

        valid_xml = os.path.join(
            os.path.dirname(__file__),
            "mikado.blast.xml"
        )

        valid, header, exc = blast_utils.BlastOpener(valid_xml).sniff()
        self.assertTrue(valid, (valid, exc))
        self.assertIsNone(exc, exc)

        with open(valid_xml, mode="rt") as new_handle:
            valid, header, exc = blast_utils.BlastOpener(new_handle).sniff()
            self.assertTrue(valid, (valid, exc))
            self.assertIsNone(exc, exc)

    def test_sniff_invalid(self):

        invalid_xml = tempfile.NamedTemporaryFile(delete=False)
        invalid_xml.write(b"failing\n")
        invalid_xml.close()
        with self.assertRaises(ValueError):
            valid, header, exc = blast_utils.BlastOpener(invalid_xml.name).sniff()

        os.remove(invalid_xml.name)

    def test_sniff_inexistent(self):

        inexistent_xml = tempfile.mktemp()
        with self.assertRaises(OSError):
            valid, header, exc = blast_utils.BlastOpener(inexistent_xml).sniff()

    def test_sniff_gzip(self):

        new = gzip.open(tempfile.mktemp(suffix=".xml.gz"), mode="wt")
        valid_xml = os.path.join(
            os.path.dirname(__file__),
            "mikado.blast.xml"
        )
        with open(valid_xml) as vx:
            for line in vx:
                new.write(line)
        new.flush()

        valid, header, exc = blast_utils.BlastOpener(new.name).sniff()
        self.assertTrue(valid, (valid, exc))
        self.assertIsNone(exc, exc)

        with gzip.open(new.name, mode="rt") as new_handle:
            valid, header, exc = blast_utils.BlastOpener(new_handle).sniff()
            self.assertTrue(valid, (valid, exc))
            self.assertIsNone(exc, exc)

    def test_fail_closed(self):

        valid_xml = os.path.join(
            os.path.dirname(__file__),
            "mikado.blast.xml"
        )

        opener = blast_utils.BlastOpener(valid_xml)

        opener.close()
        with self.assertRaises(ValueError):
            opener.sniff()

    @unittest.skip
    def test_asn(self):

        master = os.getcwd()
        os.chdir(os.path.dirname(__file__))

        valid_asn = "mikado.blast.asn"

        with gzip.open("{0}.gz".format(valid_asn), "rt") as comp_asn:
            with open(valid_asn, "wt") as asn:
                for line in comp_asn:
                    asn.write(line)

        with open("uniprot_sprot_plants.fasta", "wt") as uni_out:
            with gzip.open("uniprot_sprot_plants.fasta.gz", "rt") as uni:
                for line in uni:
                    uni_out.write(line)
        subprocess.call(
            "makeblastdb -in uniprot_sprot_plants.fasta -dbtype=prot".format(
                os.path.dirname(__file__)
            ),
            shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

        valid, header, exc = blast_utils.BlastOpener(valid_asn).sniff()
        self.assertTrue(valid, (valid, exc))
        self.assertIsNone(exc, exc)

        valid, header, exc = blast_utils.BlastOpener("{0}.gz".format(valid_asn)).sniff()
        self.assertTrue(valid, (valid, exc))
        self.assertIsNone(exc, exc)
        with blast_utils.BlastOpener("{0}.gz".format(valid_asn)) as gzasn:
            self.assertFalse(gzasn.closed)
            self.assertTrue(gzasn.open)
            record = next(gzasn)

        # This should fix the issue of blast_formatter crashing
        time.sleep(1)

        for fname in os.listdir("."):
            if "uniprot_sprot_plants.fasta" in fname and not fname.endswith(".gz"):
                os.remove(fname)

        os.remove(valid_asn)
        os.chdir(master)


class TestMerging(unittest.TestCase):

    """Small class to test basic cases of the merging algorithm."""

    def test_merging_1(self):

        l = [(-10, -5), (-6, 8), (5, 10), (20, 40)]
        tot_length = 51
        corr_merged = [(-10, 10), (20, 40)]
        merged, tot_length = seri_blast_utils.merge(l, query_length=tot_length, offset=1)
        self.assertTrue((merged == corr_merged))
        self.assertEqual(tot_length, 10 - -10 +1 + 40 - 20 + 1)

    def test_merging_2(self):

        l = [(100, 200)]
        for offset in [0, 1, 2]:
            with self.subTest(offset=offset):
                tot_length = l[0][1] - l[0][0] + offset
                if offset == 2:
                    with self.assertRaises(ValueError):
                        _ = seri_blast_utils.merge(l, offset=offset)
                else:
                    merged, length = seri_blast_utils.merge(l, offset=offset)
                    self.assertEqual(length, tot_length)
                    self.assertTrue((merged == l), (merged, l))

    def test_various_merging(self):

        invalid = [
            [('a', 0)],
            [('a', 'b')],
            [(10, 20), ('a', 'b')],
            [(10, 20, 30), (40, 50, 60)]
        ]

        for inv in invalid:
            with self.subTest(inv=inv):
                with self.assertRaises(TypeError, msg=inv):
                    seri_blast_utils.merge(inv)

        valid = {
            0: [[('10', '20')], [(10, 20)]],
            1: [[(10, 30)], [(10, 30)]],
            2: [[(-10.0, 5.5)], [(-10, 5)]],
            3: [[(-4, -10)], [(-10, -4)]],
            4: [[(-5, -10), (-2.2, -7.3)], [(-10, -2)]]
        }

        for val in valid:
            inp, out = valid[val]
            with self.subTest(val=val, msg=valid[val]):
                _ = seri_blast_utils.merge(inp)
                self.assertTrue((out == _[0]), (out, _[0]))

    def test_included(self):

        cases = {
            tuple([(10, 60), (40, 100), (200, 400)]): [(10, 100), (200, 400)],
            tuple([(54, 1194), (110, 790), (950, 1052)]): [(54, 1194)],
            tuple([(54, 1194), (110, 790), (950, 1052), (1200, 1400)]): [(54, 1194), (1200, 1400)]
        }
        for val, out in cases.items():
            with self.subTest(val=val, msg=cases[val]):
                _ = seri_blast_utils.merge(list(val))
                self.assertTrue((out == _[0]), (out, _[0]))

    def test_unordered(self):
        cases = {
            tuple([(10, 60), (40, 100), (200, 400)]): [(10, 100), (200, 400)],
            tuple([(54, 1194), (110, 790), (950, 1052)]): [(54, 1194)],
            tuple([(54, 1194), (110, 790), (950, 1052), (1200, 1400)]): [(54, 1194), (1200, 1400)]
        }

        from random import shuffle
        for num in range(30):
            for val, out in cases.items():
                cval = list(val[:])
                shuffle(cval)
                with self.subTest(val=val, cval=cval, msg=cases[val]):
                    _ = seri_blast_utils.merge(list(cval))
                    self.assertTrue((out == _[0]), (out, _[0]))


if __name__ == '__main__':
    unittest.main()