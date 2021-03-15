import re
import shutil
import unittest

import os

from Mikado import create_default_logger
from Mikado.exceptions import CorruptIndex, InvalidParsingFormat
from Mikado.parsers import parser_factory
from Mikado.scales.reference_preparation import indexing
from Mikado.scales.assignment import Assigner
import pkg_resources
import tempfile

from Mikado.transcripts import Transcript
from Mikado.utilities.namespace import Namespace


class TestIndexing(unittest.TestCase):

    # def setUp(self) -> None:

    def test_not_existent(self):
        logger = create_default_logger("test_not_existent")
        with tempfile.NamedTemporaryFile() as foo, self.assertRaises(CorruptIndex):
            indexing.check_index(foo.name, logger)
        with tempfile.NamedTemporaryFile() as foo, self.assertRaises(CorruptIndex):
            namespace = Namespace(default=False)
            namespace.reference = open(foo.name)
            indexing.load_index(namespace, logger)

    def test_index_various(self):
        logger = create_default_logger("test_index_gff3")
        for fname in ("trinity.gff3", "trinity.gtf", "trinity.bed12",
                      "trinity.cDNA_match.gff3", "trinity.match_matchpart.gff3"):
            with self.subTest(fname=fname), tempfile.TemporaryDirectory() as out_folder:
                shutil.copy(pkg_resources.resource_filename("Mikado.tests", fname), out_folder)
                fhandle = os.path.join(out_folder, fname)
                index_name = "{}.midx".format(fhandle)
                indexing.create_index(parser_factory(fhandle), logger, index_name,
                                      ref_gff=(fname.endswith("gff3")))
                self.assertTrue(os.path.exists(index_name))
                self.assertGreater(os.stat(index_name).st_size, 0)
                # Now check the index has been created correctly
                indexing.check_index(index_name, logger)
                namespace = Namespace(default=False)
                namespace.reference = parser_factory(fhandle)
                _ = indexing.load_index(namespace, logger)
                # Now rebuild
                with self.assertLogs(logger, level="INFO") as cmo:
                    indexing.create_index(parser_factory(fhandle), logger, index_name,
                                          ref_gff=(fname.endswith("gff3")))
                self.assertTrue(any([re.search(r"Removing the old index", _.msg) for _ in cmo.records]),
                                cmo.records)

    def test_protein_coding(self):
        fname = pkg_resources.resource_filename("Mikado.tests", "reference.gff3")

    def test_fusion(self):

        t = Transcript()
        t.chrom, t.strand, t.start, t.end, t.id, t.parent = "Chr1", "+", 101, 1000, "foo.1", "foo"
        t.add_exons([(101, 500), (601, 800), (901, 1000)])
        t.finalize()
        t2 = Transcript()
        t2.chrom, t2.strand, t2.start, t2.end, t2.id, t2.parent = "Chr1", "+", 2001, 3000, "bar.1", "bar"
        t2.add_exons([(2001, 2500), (2601, 2800), (2901, 3000)])
        t2.finalize()

        t3 = Transcript()
        t3.chrom, t3.strand, t3.start, t3.end, t3.id, t3.parent = "Chr1", "+", 651, 2703, "faz.1", "faz"
        t3.add_exons([(651, 800), (901, 1300), (2230, 2500), (2601, 2703)])
        t3.finalize()

        logger = create_default_logger("test_fusion")
        with tempfile.TemporaryDirectory() as folder:
            with open(os.path.join(folder, "reference.gtf"), "wt") as reference:
                print(t.format("gtf"), file=reference)
                print(t2.format("gtf"), file=reference)
            self.assertTrue(os.path.exists(reference.name))
            _ = [_ for _ in parser_factory(reference.name)]
            try:
                indexing.create_index(parser_factory(reference.name), logger, "{}.midx".format(reference.name))
            except InvalidParsingFormat:
                self.assertFalse(True, "\n".join([line.rstrip() for line in open(reference.name)]))
            namespace = Namespace(default=False)
            namespace.out = os.path.join(folder, "out")
            for report in (False, True):
                with self.subTest(report=report):
                    namespace.report_fusions = report
                    assigner = Assigner("{}.midx".format(reference.name), args=namespace, printout_tmap=False)
                    result = assigner.get_best(t3)
                    if report:
                        self.assertTrue(len(result), 2)
                        self.assertTrue(result[0].ccode == ("f", "j"), str(result[0]))
                        self.assertTrue(result[1].ccode == ("f", "j"), str(result[1]))
                    else:
                        self.assertTrue(result.ccode == ("j",), str(result))


if __name__ == '__main__':
    unittest.main()
