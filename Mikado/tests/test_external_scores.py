import re
import unittest

from .. import create_default_logger
from .._transcripts.scoring_configuration import ScoringFile, MinMaxScore
from ..loci import Transcript, Superlocus
import pkg_resources
import os
from ..configuration.configurator import load_and_validate_config
from ..serializers.external import ExternalSerializer
import pandas as pd


class ExternalTester(unittest.TestCase):

    def setUp(self):
        self.conf = load_and_validate_config(None)
        self.transcript = Transcript()
        self.transcript.chrom = "15"
        self.transcript.source = "protein_coding"
        self.transcript.start = 47631264
        self.transcript.end = 48051999

        exons = [(47631264, 47631416),
                 (47704590, 47704669),
                 (47762671, 47762742),
                 (47893062, 47893093),
                 (47895572, 47895655),
                 (48051942, 48051999)]

        self.transcript.strand = "+"
        self.transcript.add_exons(exons)
        self.transcript.id = "ENST00000560636"
        self.transcript.parent = "ENSG00000137872"
        self.transcript2 = self.transcript.copy()
        self.transcript2.id = "ENST00000560637"
        # self.assertIn("scoring", self.conf)

    def test_copying(self):
        self.transcript.external_scores.update({"test": 0, "test1": 1})
        self.assertEqual(self.transcript.external_scores.test, 0)
        self.assertEqual(self.transcript.external_scores.test1, 1)
        transcript = self.transcript.deepcopy()
        self.assertEqual(transcript.external_scores.test, 0)
        self.assertEqual(transcript.external_scores.test1, 1)

    def test_no_handle(self):
        logger = create_default_logger("test_no_handle", level="DEBUG")
        with self.assertLogs(logger.name, level="WARNING") as cmo:
            ExternalSerializer(logger=logger, configuration=load_and_validate_config(None),
                               handle=None)
        self.assertTrue(any([record.msg == "No input file specified. Exiting." for record in cmo.records]))

    def test_no_fasta(self):
        conf = load_and_validate_config(None)
        conf.serialise.files.transcripts = "foo"
        logger = create_default_logger("test_no_fasta", level="DEBUG")
        with self.assertRaises(AssertionError), self.assertLogs(logger.name, level="CRITICAL") as cmo:
            ExternalSerializer(logger=logger, configuration=load_and_validate_config(None),
                               handle="trinity.bed12")  # TODO change
        self.assertTrue(any([record.msg == """I cannot find the mikado prepared FASTA file with the transcripts to analyse.
Please run mikado serialise in the folder with the correct files, and/or modify the configuration or the command line options."""
                             for record in cmo.records]))

    def test_invalid_tsv(self):
        conf = load_and_validate_config(None)
        conf.serialise.files.transcripts = pkg_resources.resource_filename("Mikado.tests",
                                                                           os.path.join("test_external",
                                                                                        "mikado_prepared.fasta"))
        self.assertTrue(os.path.exists(conf.serialise.files.transcripts),
                        conf.serialise.files.transcripts)
        logger = create_default_logger("test_no_fasta", level="DEBUG")
        h5 = pkg_resources.resource_filename("Mikado.tests", os.path.join("test_external", "abundance.h5"))
        with self.assertRaises((pd.errors.ParserError,)):
            ExternalSerializer(logger=logger, configuration=conf, handle=h5)
        wrong_tsv_one = pkg_resources.resource_filename("Mikado.tests", os.path.join("test_external", "wrong_1.tsv"))
        with self.assertRaises((ValueError, pd.errors.ParserError)), \
                self.assertLogs(logger.name, level="CRITICAL") as cmo:
            ExternalSerializer(logger=logger, configuration=conf, handle=wrong_tsv_one)
        self.assertTrue(any([re.search(r"Invalid input file", str(record.msg)) is not None for record in cmo.records]))
        wrong_tsv_two = pkg_resources.resource_filename("Mikado.tests", os.path.join("test_external", "wrong_2.tsv"))
        with self.assertRaises((ValueError, pd.errors.ParserError)), \
                self.assertLogs(logger.name, level="CRITICAL") as cmo:
            ExternalSerializer(logger=logger, configuration=conf, handle=wrong_tsv_two)
        self.assertTrue(any([re.search(r"Invalid input file", str(record.msg)) is not None for record in cmo.records]))

    def test_real(self):
        self.transcript.attributes["tpm"] = 10

        checked_conf = self.conf.copy()
        self.assertTrue(hasattr(self.conf.scoring.requirements, "parameters"))
        self.assertTrue(hasattr(checked_conf.scoring.requirements, "parameters"))
        checked_conf.scoring.scoring["attributes.tpm"] = MinMaxScore.Schema().load(
            {"rescaling": "max", "default": 0, "rtype": "float"})
        # self.assertIn('attributes.tpm', checked_conf.scoring.scoring)

        logger = create_default_logger(name="test_real", level="WARNING")
        sup = Superlocus(self.transcript, configuration=checked_conf, logger=logger)
        self.assertIn("attributes.tpm", sup.configuration.scoring.scoring)
        sup.get_metrics()
        self.assertIn("attributes.tpm", sup._metrics[self.transcript.id])
        self.assertEqual(sup._metrics[self.transcript.id]["attributes.tpm"], 10)

        sup2 = Superlocus(self.transcript2, configuration=checked_conf)
        sup2.get_metrics()
        self.assertIn("attributes.tpm", sup2._metrics[self.transcript2.id])
        self.assertEqual(sup2._metrics[self.transcript2.id]["attributes.tpm"], 0)

    def test_multiplier(self):

        self.transcript.attributes["tpm"] = 10

        checked_conf = self.conf.copy()
        checked_conf.scoring.scoring["attributes.tpm"] = MinMaxScore.Schema().load(
            {"rescaling": "max", "default": 0, "rtype": "float", 'multiplier': 4})

        self.assertIn('attributes.tpm', checked_conf.scoring.scoring)

        sup = Superlocus(self.transcript, configuration=checked_conf)
        tid = self.transcript.id
        self.assertIn(tid, sup.transcripts)
        sup.filter_and_calculate_scores(check_requirements=False)
        self.assertEqual(sup.scores[tid]['attributes.tpm'], 4)

    def test_default_attribute_score(self):
        self.transcript.attributes["foo"] = True

        checked_conf = self.conf.copy()
        checked_conf.scoring.scoring["attributes.foo"] = MinMaxScore.Schema().load(
            {"rescaling": "max", "rtype": "bool", "default": False})

        self.assertIn('attributes.foo', checked_conf.scoring.scoring)

        sup = Superlocus(self.transcript, configuration=checked_conf)
        sup.get_metrics()
        self.assertIn("attributes.foo", sup._metrics[self.transcript.id])
        self.assertEqual(sup._metrics[self.transcript.id]["attributes.foo"], True)

        sup2 = Superlocus(self.transcript2, configuration=checked_conf)
        sup2.get_metrics()
        self.assertIn("attributes.foo", sup2._metrics[self.transcript2.id])
        self.assertEqual(sup2._metrics[self.transcript2.id]["attributes.foo"], False)

    def test_error_attribute(self):
        self.transcript.attributes["tpm"] = "10a"
        checked_conf = self.conf.copy()
        checked_conf.scoring.scoring["attributes.tpm"] = MinMaxScore.Schema().load(
            {"rescaling": "max", "default": 0, "rtype": "float", 'multiplier': 4})
        self.assertIn('attributes.tpm', checked_conf.scoring.scoring)
        sup = Superlocus(self.transcript, configuration=checked_conf)
        with self.assertRaises(ValueError):
            sup.get_metrics()

    def test_attribute_use_raw(self):

        self.transcript.attributes["tpm"] = 10

        checked_conf = self.conf.copy()
        checked_conf.scoring.scoring["attributes.tpm"] = MinMaxScore.Schema().load(
            {"rescaling": "max", "default": 0, "rtype": "float", 'multiplier': 4,
             'use_raw': True})

        self.assertIn('attributes.tpm', checked_conf.scoring.scoring)

        sup = Superlocus(self.transcript, configuration=checked_conf)
        tid = self.transcript.id
        self.assertIn(tid, sup.transcripts)
        with self.assertRaises(ValueError):
            sup.get_metrics()

    def test_attribute_use_raw_percentage(self):

        self.transcript.attributes["tpm"] = 10

        checked_conf = self.conf.copy()
        checked_conf.scoring.scoring["attributes.tpm"] = MinMaxScore.Schema().load(
            {"rescaling": "max", "default": 0, "rtype": "float", 'multiplier': 4,
             'use_raw': True, 'percentage': True})

        self.assertIn('attributes.tpm', checked_conf.scoring.scoring)

        sup = Superlocus(self.transcript, configuration=checked_conf)
        tid = self.transcript.id
        self.assertIn(tid, sup.transcripts)
        sup.get_metrics()
        self.assertIn("attributes.tpm", sup._metrics[self.transcript.id])
        self.assertEqual(sup._metrics[self.transcript.id]["attributes.tpm"], 0.1)
