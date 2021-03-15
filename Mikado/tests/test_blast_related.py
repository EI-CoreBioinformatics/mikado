#!/usr/bin/env python3


import tempfile

import pkg_resources

from Mikado.exceptions import InvalidParsingFormat
from ..parsers import blast_utils
import unittest
import os
import gzip
import subprocess
from ..serializers.blast_serializer import xml_utils as seri_blast_utils
from ..serializers.blast_serializer.tabular_utils import matrices
from ..serializers.blast_serializer import parse_btop
import numpy as np
import time
import itertools
from pytest import mark
from collections import namedtuple
from shutil import which


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
                        with self.subTest(qsize=qsize, ssize=ssize, qpos=qpos, spos=spos, mlength=mlength):
                            qar, sar = np.zeros([3, qsize], dtype=int), np.zeros([3, ssize], dtype=int)
                            sm = str(mlength)
                            qar, sar, tot, match = parse_btop(sm, qpos, spos, qar, sar, matrix, qmult=qmult, tmult=tmult)
                            qfound = np.where(qar > 0)
                            sfound = np.where(sar > 0)
                            self.assertEqual(tot, mlength * min(qmult, tmult),
                                             (tot, mlength, qpos, spos, qmult, tmult))
                            self.assertTrue(((qfound[1][qfound[0] == 0] == qfound[1][qfound[0] == 1]) &
                                 (qfound[1][qfound[0] == 0] == qfound[1][qfound[0] == 2])).all()
                            )
                            self.assertTrue(((sfound[1][sfound[0] == 0] == sfound[1][sfound[0] == 1]) &
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

    def test_only_match(self):
        for qmult, tmult in itertools.product([1, 3], [1, 3]):
            for qsize, ssize in itertools.product(range(10 * qmult, 20 * qmult, 10),
                                                  range(10 * tmult, 20 * tmult, 10)):
                for qpos, spos in itertools.product(range(0, qsize, 2), range(0, ssize, 2)):
                    mqlength = (qsize - qpos) // qmult
                    mslength = (ssize - spos) // tmult
                    for mlength in range(1, int(min(mqlength, mslength))):
                        for score in (-1, 0, 1):
                            with self.subTest():
                                match = "AT" * mlength
                                qar, sar = np.zeros([3, qsize], dtype=int), np.zeros([3, ssize], dtype=int)
                                qar, sar, tot, match = parse_btop(match, qpos, spos, qar, sar, {"AT": score},
                                                           qmult=qmult, tmult=tmult)
                                qfound = np.where(qar > 0)
                                sfound = np.where(sar > 0)
                                self.assertEqual(tot, mlength * min(qmult, tmult),
                                                 (tot, mlength, qpos, spos, qmult, tmult))
                                self.assertTrue(qfound[1][qfound[0] == 0].shape[0] == mlength * qmult)
                                self.assertTrue(sfound[1][sfound[0] == 0].shape[0] == mlength * tmult)
                                self.assertTrue(qfound[1][qfound[0] == 1].shape[0] == 0)
                                self.assertTrue(sfound[1][sfound[0] == 1].shape[0] == 0)
                                if score <= 0:
                                    self.assertTrue(qfound[1][qfound[0] == 2].shape[0] == 0)
                                    self.assertTrue(sfound[1][sfound[0] == 2].shape[0] == 0)
                                else:  # Positives found!
                                    self.assertTrue((qfound[1][qfound[0] == 0] == qfound[1][qfound[0] == 2]).all(),
                                                    (score, qfound))
                                    self.assertTrue((sfound[1][sfound[0] == 0] == sfound[1][sfound[0] == 2]).all())
                                # self.assertEqual(qfound[0].shape[0], 3)
                                self.assertEqual(np.where(qar[0] > 0)[0].min(), qpos, (qfound[0], qpos, mlength))
                                self.assertEqual(np.where(qar[0] > 0)[0].max(), qpos + mlength * qmult - 1)
                                self.assertEqual(np.where(sar[0] > 0)[0].min(), spos)
                                self.assertEqual(np.where(sar[0] > 0)[0].max(), spos + mlength * tmult - 1)

    @mark.slow
    def test_gap(self):

        for qmult, tmult in itertools.product([1, 3], [1, 3]):
            for qsize, ssize in itertools.product(range(10 * qmult, 20 * qmult, 5),
                                                  range(10 * tmult, 20 * tmult, 5)):
                for qpos, spos in itertools.product(range(0, qsize, 5), range(0, ssize, 5)):
                    mqlength = (qsize - qpos) // qmult
                    mslength = (ssize - spos) // tmult
                    for mlength in range(3, int(min(mqlength, mslength))):
                        for gap_pos in range(mlength - 1):
                            for gap in ("A-", "-A"):
                                sm = ""
                                if gap_pos - 1 > -1:
                                    sm += str(gap_pos)
                                sm += gap
                                if mlength - gap_pos - 1:
                                    sm += str(mlength - gap_pos)
                                qar, sar = np.zeros([3, qsize], dtype=int), np.zeros([3, ssize], dtype=int)
                                qar, sar, tot, match = parse_btop(sm, qpos, spos, qar, sar, dict(), qmult=qmult, tmult=tmult)
                                qfound = np.where(qar > 0)
                                sfound = np.where(sar > 0)
                                self.assertEqual(tot, mlength * min(qmult, tmult) + min(qmult, tmult),
                                                 (sm, tot, mlength, gap_pos, qpos, spos, qmult, tmult))
                                if gap[0] == "-":  # gap in query
                                    qadd, sadd = 0, tmult
                                    self.assertIn(spos + gap_pos * tmult, np.where(sar[1] == 0)[0],
                                                  (sm, mlength, spos, gap_pos, tmult, sar))
                                else:
                                    qadd, sadd = qmult, 0
                                    self.assertIn(qpos + gap_pos * qmult, np.where(qar[1] == 0)[0])

                                self.assertEqual(np.where(qar[0] > 0)[0].shape[0], mlength * qmult + qadd)
                                self.assertEqual(np.where(sar[0] > 0)[0].shape[0], mlength * tmult + sadd,
                                                 (sar, sm))
                                self.assertEqual(np.where(sar[1] > 0)[0].shape[0], mlength * tmult)
                                self.assertEqual(np.where(qar[1] > 0)[0].shape[0], mlength * qmult)


class XMLLineTester(unittest.TestCase):

    class HSP(object):
        seq_t = namedtuple("seq", ["seq"])

        def __init__(self, query=None, target=None, mid=None,
                     query_start=0, query_end=0, ident_num=0, pos_num=0,
                     frame=0):
            self.aln_annotation = dict()
            self.query_start, self.query_end = query_start, query_end
            self.ident_num, self.pos_num = ident_num, pos_num
            self.query_frame = frame
            self.set_seqs(query, target, mid)

        def set_seqs(self, query, target, mid):
            self.query = self.seq_t(query)
            self.hit = self.seq_t(target)
            self.aln_annotation["similarity"] = mid

    def test_aln_string_basic(self):
        _qlenght, _slength = 20, 20
        for qpos, spos in itertools.product(range(1, _qlenght), range(1, _qlenght)):
            for qmult, tmult in itertools.product((1, 3), (1, 3)):
                qlength, slength = _qlenght * qmult, _slength * tmult
                max_length = min((qlength - qpos) // qmult,
                                 (slength - spos) // tmult)
                for l in range(1, max_length + 1):
                    qseq, sseq, mid = "A" * l, "A" * l, "A" * l
                    hs = self.HSP(query=qseq, target=sseq, mid=mid,
                                  query_start=qpos, query_end=qpos + l * qmult,
                                  ident_num=l, pos_num=l,
                                  frame=0)
                    match, ident, pos = seri_blast_utils.prepare_aln_strings(hs, qmultiplier=qmult)
                    self.assertEqual(match, "|" * l)
                    self.assertTrue((ident == pos).all())


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

        with tempfile.NamedTemporaryFile(delete=True) as invalid_xml:
            invalid_xml.write(b"failing\n")
            invalid_xml.flush()
            with self.assertRaises((ValueError, InvalidParsingFormat)):
                _ = blast_utils.BlastOpener(invalid_xml.name).sniff()

    def test_sniff_inexistent(self):

        with tempfile.NamedTemporaryFile(delete=True) as inexistent:
            pass

        with self.assertRaises(OSError):
            _ = blast_utils.BlastOpener(inexistent.name).sniff()

    @mark.slow
    def test_sniff_gzip(self):

        with tempfile.NamedTemporaryFile(suffix=".xml.gz") as handle:
            new = gzip.open(handle.name, mode="wb")

            with pkg_resources.resource_stream("Mikado.tests", "mikado.blast.xml") as valid_xml:
                for line in valid_xml:
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

    @unittest.skipUnless(which("blast_formatter") is not None and which("makeblastdb") is not None,
                         "This test can be performed only if NCBI+ is installed")
    @mark.slow
    def test_asn(self):

        folder = tempfile.TemporaryDirectory()

        # valid_asn = tempfile.NamedTemporaryFile(suffix=".asn", mode="wb")
        valid_asn = os.path.join(folder.name, "mikado.blast.asn")
        comp_asn = os.path.join(folder.name, "mikado.blast.asn.gz")
        with pkg_resources.resource_stream("Mikado.tests", "mikado.blast.asn.gz") as compressed_asn, \
                open(valid_asn, "wb") as valid_out, open(comp_asn, "wb") as comp_out:
            stream = compressed_asn.read()
            comp_out.write(stream)
            valid_out.write(gzip.decompress(stream))

        uni_out = os.path.join(folder.name, "uniprot_sprot_plants.fasta")
        with pkg_resources.resource_stream("Mikado.tests", "uniprot_sprot_plants.fasta.gz") as uni, \
                open(uni_out, "wb") as uni_handle:
            uni_handle.write(gzip.decompress(uni.read()))

        master = os.getcwd()
        os.chdir(folder.name)
        subprocess.call(
            "makeblastdb -in {uni_out} -dbtype=prot".format(uni_out=uni_out),
            shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

        valid, header, exc = blast_utils.BlastOpener(valid_asn).sniff()
        self.assertTrue(valid, (valid, exc))
        self.assertIsNone(exc, exc)

        valid, header, exc = blast_utils.BlastOpener("mikado.blast.asn.gz").sniff()
        self.assertTrue(valid, (valid, exc))
        self.assertIsNone(exc, exc)
        with blast_utils.BlastOpener(comp_asn) as gzasn:
            self.assertFalse(gzasn.closed)
            self.assertTrue(gzasn.open)
            _ = next(gzasn)

        # This should fix the issue of blast_formatter crashing
        time.sleep(1)
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