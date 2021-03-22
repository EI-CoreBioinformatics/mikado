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
from Mikado.configuration import print_config, configurator


@mark.slow
class LightTest(unittest.TestCase):
    def test_subprocess_multi_empty_orfs(self):
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
                json_file = os.path.join(folder, "mikado.yaml")
                db = os.path.join(folder, "mikado.db")
                log = "failed_serialise.log"
                uni_out = os.path.join(folder, "uniprot_sprot_plants.fasta")
                self.configuration.serialise.files.log = os.path.basename(log)
                self.configuration.multiprocessing_method = "fork"
                with gzip.open(uniprot, "rb") as uni, open(uni_out, "wb") as uni_out_handle:
                    uni_out_handle.write(uni.read())
                with open(json_file, "wt") as json_handle:
                    print_config(self.configuration, json_handle, output_format="yaml")

                self.assertFalse(os.path.exists(db))
                sys.argv = [str(_) for _ in ["mikado", "serialise", "--json-conf", json_file,
                                             "--transcripts", transcripts, "--blast_targets", uni_out,
                                             "--log", log,
                                             "-od", folder,
                                             "--force",
                                             "--orfs", tmp_orf.name, "--junctions", junctions, "--xml", xml,
                                             "-p", procs, "-mo", mobjects, db,
                                             "--seed", "1078"]]
                log = os.path.join(folder, log)
                subprocess.call(sys.argv, shell=False)
                # with self.assertRaises(SystemExit):
                #     pkg_resources.load_entry_point("Mikado", "console_scripts", "mikado")()
                self.assertTrue("failed" in log)
                self.assertTrue(os.path.exists(log), log)
                self.assertTrue(os.stat(log).st_size > 0, log)
                logged = [_.rstrip() for _ in open(log)]
                self.assertGreater(len(logged), 0)
                self.assertFalse(os.path.exists(db), logged)
                self.assertTrue(any(
                    "Mikado serialise failed due to problems with the input data. Please check the logs." in line
                    for line in logged), [print(line) for line in logged])
                self.assertTrue(any(
                    "The provided ORFs do not match the transcripts provided and "
                    "already present in the database." in line for line in logged),
                    print("\n".join(logged)))
