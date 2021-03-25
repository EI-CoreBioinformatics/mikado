import unittest
from shutil import which

import pkg_resources
import os
import tempfile

import pytest
import sqlalchemy

from Mikado.serializers.blast_serializer import Hit
from Mikado.utilities import dbutils

from Mikado.serializers import blast_serializer
import sqlalchemy.orm
from Mikado import create_default_logger
from Mikado.configuration import MikadoConfiguration
from Mikado.subprograms.serialise import load_blast


class TestBlastParsing(unittest.TestCase):

    def setUp(self) -> None:
        self.master = pkg_resources.resource_filename("Mikado.tests", "blast_input_tests")
        self.queries = os.path.join(self.master, "queries.fasta")
        self.assertTrue(os.path.exists(self.queries))

    def run_loading(self, target, out_db, logger, conf):
        with self.assertLogs(logger, level="DEBUG") as cmo:
            try:
                blast_serializer.BlastSerializer(
                    target,
                    configuration=conf,
                    logger=logger)()
            except (ValueError, KeyError, TypeError) as exc:
                print(*cmo.output, sep="\n")
                raise exc
        self.assertTrue(os.path.exists(out_db), cmo.output)
        self.assertGreater(os.stat(out_db).st_size, 0)
        conf.db_settings.db = out_db
        engine = dbutils.connect(configuration=conf)
        Session = sqlalchemy.orm.sessionmaker()
        session = Session(bind=engine)
        self.assertTrue(session.query(Hit).count() > 0)

    @unittest.skip
    @pytest.mark.triage
    @unittest.skipUnless(which("blast_formatter") is not None, reason="NCBI BLAST+ not installed")
    def test_asn(self):
        # Currently DISABLED because the ASN specifications requires the database to be where indicated by the
        # relative path within the ASN. So for the time being this test is *not active*.

        asns = [os.path.join("blast", "asn", "blast.asn.gz"),
                os.path.join("blast_parse_seqids", "asn", "blast.asn.gz")]

        for folder in ["sanitised", "uniprot"]:
            conf = MikadoConfiguration()
            targets = os.path.join(self.master, folder, "uniprot.fasta")
            self.assertTrue(os.path.exists(targets))
            for asn in asns:
                asn = os.path.join(self.master, folder, asn)
                self.assertTrue(os.path.exists(asn))
                with tempfile.TemporaryDirectory() as out_folder:
                    conf.serialise.files.output_dir = out_folder
                    conf.serialise.files.blast_targets = [targets]
                    conf.serialise.files.transcripts = self.queries
                    conf.serialise.files.xml = [asn]
                    logger = create_default_logger(f"test_asn_{folder}")
                    out_db = os.path.join(out_folder, conf.db_settings.db)
                    conf.db_settings.db = out_db
                    self.run_loading(asn, out_db, logger, conf)

    @unittest.skipUnless(which("diamond") is not None, reason="DIAMOND not installed")
    def test_daa(self):
        daa_base = os.path.join("diamond", "daa", "blast.daa")
        for folder in ["sanitised", "uniprot"]:
            targets = os.path.join(self.master, folder, "uniprot.fasta")
            self.assertTrue(os.path.exists(targets))
            daa = os.path.join(self.master, folder, daa_base)
            self.assertTrue(os.path.exists(daa))
            with tempfile.TemporaryDirectory() as out_folder:
                conf = MikadoConfiguration()
                out_db = os.path.join(out_folder, conf.db_settings.db)
                conf.db_settings.db = out_db
                conf.serialise.files.output_dir = out_folder
                conf.serialise.files.blast_targets = [targets]
                conf.serialise.files.transcripts = self.queries
                conf.serialise.files.xml = [daa]
                logger = create_default_logger(f"test_daa_{folder}")
                self.assertTrue(tempfile.gettempdir() in out_db)
                self.run_loading(daa, out_db, logger, conf)

    def test_xml(self):
        xml_base = os.path.join("xml", "blast.xml.gz")
        for folder in ["sanitised", "uniprot"]:
            for subfolder in ["blast", "blast_parse_seqids", "diamond"]:
                targets = os.path.join(self.master, folder, "uniprot.fasta")
                self.assertTrue(os.path.exists(targets))
                xml = os.path.join(self.master, folder, subfolder, xml_base)
                self.assertTrue(os.path.exists(xml))
                with tempfile.TemporaryDirectory() as out_folder:
                    conf = MikadoConfiguration()
                    conf.serialise.files.output_dir = out_folder
                    out_db = os.path.join(out_folder, conf.db_settings.db)
                    conf.db_settings.db = out_db
                    conf.serialise.files.blast_targets = [targets]
                    conf.serialise.files.transcripts = self.queries
                    conf.serialise.files.xml = [xml]
                    logger = create_default_logger(f"test_xml_{folder}_{subfolder}")
                    self.assertTrue(tempfile.gettempdir() in out_db)
                    self.run_loading(xml, out_db, logger, conf)

    def test_tsv(self):
        tsv_base = os.path.join("tsv", "blast.tsv.gz")
        for folder in ["uniprot", "sanitised"]:
            for subfolder in ["blast", "blast_parse_seqids", "diamond"]:
                targets = os.path.join(self.master, folder, "uniprot.fasta")
                self.assertTrue(os.path.exists(targets))
                tsv = os.path.join(self.master, folder, subfolder, tsv_base)
                self.assertTrue(os.path.exists(tsv))
                with tempfile.TemporaryDirectory() as out_folder:
                    conf = MikadoConfiguration()
                    conf.serialise.files.output_dir = out_folder
                    conf.serialise.files.blast_targets = [targets]
                    conf.serialise.files.transcripts = self.queries
                    conf.serialise.files.xml = [tsv]
                    out_db = os.path.join(out_folder, conf.db_settings.db)
                    conf.db_settings.db = out_db
                    self.assertTrue(tempfile.gettempdir() in out_db)
                    logger = create_default_logger(f"test_tsv_{folder}_{subfolder}", level="DEBUG")
                    self.run_loading(tsv, out_db, logger, conf)


if __name__ == "__main__":
    unittest.main()
