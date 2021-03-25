"""This test get executed before anything else in the GHA."""
import gzip
import subprocess
import sys
import unittest
import os
import pkg_resources
import tempfile
import pysam
from pytest import mark

from Mikado import create_default_logger
from Mikado.configuration import print_config, configurator
from Mikado.exceptions import InvalidSerialization
from Mikado.subprograms.serialise import load_orfs


@mark.slow
class LightTest(unittest.TestCase):
    def test_subprocess_multi_empty_orfs(self):
        print("Started light test")
        self.fai = pysam.FastaFile(pkg_resources.resource_filename("Mikado.tests", "chr5.fas.gz"))
        self.configuration = configurator.load_and_validate_config(None)
        self.configuration.reference.genome = self.fai.filename.decode()

        xml = pkg_resources.resource_filename("Mikado.tests", "chunk-001-proteins.xml.gz")
        transcripts = pkg_resources.resource_filename("Mikado.tests", "mikado_prepared.fasta")
        junctions = pkg_resources.resource_filename("Mikado.tests", "junctions.bed")
        tmp_orf = tempfile.NamedTemporaryFile(suffix=".bed12")
        tmp_orf.write(b"#track\n")
        tmp_orf.write(
            b"cufflinks_star_at.23553.1\t0\t1733\tID=1_1;partial=01;start_type=ATG\t0\t+\t312\t1733\t0,0,0\t1\t1733\t0\n")
        tmp_orf.flush()
        uniprot = pkg_resources.resource_filename("Mikado.tests", "uniprot_sprot_plants.fasta.gz")
        mobjects = 300  # Let's test properly the serialisation for BLAST

        # Set up the command arguments
        with tempfile.TemporaryDirectory(prefix="has_to_fail") as folder_one, \
                tempfile.TemporaryDirectory(prefix="has_to_fail") as folder_two:
            for procs, folder in [(3, folder_one), (1, folder_two)]:
                self.configuration = configurator.load_and_validate_config(None)
                self.configuration.reference.genome = self.fai.filename.decode()
                json_file = os.path.join(folder, "mikado.yaml")
                db = os.path.join(folder, "mikado.db")
                log = "failed_serialise.log"
                uni_out = os.path.join(folder, "uniprot_sprot_plants.fasta")
                self.configuration.serialise.files.log = os.path.basename(log)
                self.configuration.multiprocessing_method = "fork"
                with gzip.open(uniprot, "rb") as uni, open(uni_out, "wb") as uni_out_handle:
                    uni_out_handle.write(uni.read())
                self.configuration.serialise.files.transcripts = transcripts
                self.configuration.serialise.files.blast_targets = uni_out
                self.configuration.serialise.files.log = log
                self.configuration.serialise.files.output_dir = folder
                self.configuration.serialise.force = True
                self.configuration.serialise.files.orfs = [tmp_orf.name]
                self.configuration.serialise.files.junctions = [junctions]
                self.configuration.serialise.files.xml = [xml]
                self.configuration.threads = procs
                self.configuration.serialise.max_objects = mobjects
                self.configuration.db_settings.db = db
                self.configuration.seed = 1078

                self.assertFalse(os.path.exists(db))
                logger = create_default_logger(f"test_light_serialise_{procs}", level="INFO")

                with self.assertRaises(InvalidSerialization), self.assertLogs(logger.name) as cmo:
                    load_orfs(self.configuration, logger)

                self.assertTrue(any(
                    "Mikado serialise failed due to problems with the input data. Please check the logs." in line
                    for line in cmo.output), cmo.output)
                self.assertTrue(any(
                    "The provided ORFs do not match the transcripts provided and "
                    "already present in the database." in line for line in cmo.output),
                    print("\n".join(cmo.output)))
        print("Finished light test")