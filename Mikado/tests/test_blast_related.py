#!/usr/bin/env python3


import tempfile
from Mikado.parsers import blast_utils
import unittest
import os
import gzip
import subprocess
from Mikado.serializers.blast_serializer import utils as seri_blast_utils
import time


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
        self.assertEqual(merged, corr_merged)
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
                    self.assertEqual(merged, l)

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
                self.assertEqual(out, _[0])


if __name__ == '__main__':
    unittest.main()