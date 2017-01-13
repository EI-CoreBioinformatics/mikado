import unittest
from Mikado.preparation import prepare
from Mikado.configuration import configurator
import pkg_resources
import tempfile
from Mikado.loci.transcript import Namespace
from Mikado.utilities.log_utils import create_null_logger
import logging
import gzip
import pyfaidx
import os


class PrepareChek(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.__genomefile__ = tempfile.NamedTemporaryFile(mode="wb", delete=True, suffix=".fa")

        with pkg_resources.resource_stream("Mikado.test", "chr5.fas.gz") as _:
            cls.__genomefile__.write(gzip.decompress(_.read()))
        cls.__genomefile__.flush()

        cls.trinity_res = dict((_[0], _[1]) for _ in [("tr_c73_g1_i1.mrna1.160", 286),
                                                      ("tr_c11_g1_i2.mrna1.111", 844),
                                                      ("tr_c11_g1_i1.mrna1.350", 659),
                                                      ("tr_c3_g1_i1.mrna1.716", 292),
                                                      ("tr_c3_g1_i2.mrna1.181", 213),
                                                      ("tr_c58_g1_i3.mrna1.19", 1718),
                                                      ("tr_c58_g1_i2.mrna1.35", 900),
                                                      ("tr_c58_g1_i7.mrna1.1", 504),
                                                      ("tr_c58_g1_i1.mrna1.190", 417),
                                                      ("tr_c108_g1_i1.mrna1.104", 458),
                                                      ("tr_c58_g1_i11.mrna1", 1744),
                                                      ("tr_c58_g1_i5.mrna1.3", 4707),
                                                      ("tr_c58_g1_i10.mrna1", 1655),
                                                      ("tr_c68_g1_i1.mrna1.173", 275),
                                                      ("tr_c58_g1_i11.mrna2", 404),
                                                      ("tr_c58_g1_i10.mrna2", 383),
                                                      ("tr_c58_g1_i8.mrna2", 383),
                                                      ("tr_c58_g1_i12.mrna1", 1725),
                                                      ("tr_c58_g1_i9.mrna1", 1988),
                                                      ("tr_c58_g1_i4.mrna1.6", 4862),
                                                      ("tr_c58_g1_i6.mrna1.3", 4945),
                                                      ("tr_c58_g1_i8.mrna1", 1819),
                                                      ("tr_c113_g1_i1.mrna1.94", 379),
                                                      ("tr_c77_g1_i1.mrna1.153", 635),
                                                      ("tr_c109_g1_i1.mrna1.102", 310),
                                                      ("tr_c60_g1_i1.mrna1.189", 1916),
                                                      ("tr_c115_g1_i1.mrna1.88", 243),
                                                      ("tr_c152_g1_i1.mrna1.66", 550),
                                                      ("tr_c21_g1_i1.mrna1.302", 475),
                                                      ("tr_c41_g1_i1.mrna1.224", 608),
                                                      ("tr_c6_g1_i1.mrna1.412", 1295),
                                                      ("tr_c74_g1_i1.mrna1.154", 320),
                                                      ("tr_c71_g1_i1.mrna1.167", 203),
                                                      ("tr_c114_g1_i1.mrna1.89", 206),
                                                      ("tr_c137_g1_i1.mrna1.77", 743),
                                                      ("tr_c37_g1_i2.mrna1.66", 417),
                                                      ("tr_c37_g1_i1.mrna1.234", 449),
                                                      ("tr_c120_g1_i1.mrna1.89", 269)])
        # cls.trinity_res = sorted(cls.trinity_res)
        cls.maxDiff = None

    @classmethod
    def tearDownClass(cls):
        """"""
        cls.__genomefile__.close()

    def setUp(self):

        self.conf = configurator.to_json(None)
        self.conf["reference"]["genome"] = self.__genomefile__.name
        self.logger = create_null_logger("prepare")
        self.conf["prepare"]["keep_redundant"] = True

    def tearDown(self):
        logging.shutdown()
        for fname in ["mikado_prepared.fasta", "mikado_prepared.gtf"]:
            fname = os.path.join(tempfile.gettempdir(), fname)
            if os.path.exists(fname):
                os.remove(fname)

    def test_prepare_trinity_gff(self):

        self.conf["prepare"]["files"]["labels"].append("tr")
        self.conf["prepare"]["files"]["output_dir"] = tempfile.gettempdir()
        args = Namespace()
        args.json_conf = self.conf

        for test_file in ("trinity.gff3", "trinity.match_matchpart.gff3", "trinity.cDNA_match.gff3"):
            with self.subTest(test_file=test_file):
                self.conf["prepare"]["files"]["gff"]= [pkg_resources.resource_filename("Mikado.test",
                                                                                       test_file)]

                prepare.prepare(args, self.logger)

                # Now that the program has run, let's check the output
                fa = pyfaidx.Fasta(os.path.join(self.conf["prepare"]["files"]["output_dir"],
                                                "mikado_prepared.fasta"))
                res = dict((_, len(fa[_])) for _ in fa.keys())
                fa.close()
                self.assertEqual(res, self.trinity_res)

if __name__ == "__main__":
    unittest.main()