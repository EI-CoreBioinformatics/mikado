import unittest
import Mikado.preparation.checking
import Mikado
import multiprocessing as mp
import pkg_resources
import tempfile
import logging
import logging.handlers
import gzip
import pickle
import os
import time
import pyfaidx
import re
from Mikado.tests.test_utils import ProcRunner
from queue import Queue
from sys import version_info


class MiscTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.fasta = pkg_resources.resource_filename("Mikado.tests", "chr5.fas.gz")
        cls.fasta_temp = tempfile.NamedTemporaryFile(suffix=".fa")
        with gzip.open(cls.fasta) as ffile:
            cls.fasta_temp.write(ffile.read())
        cls.fasta_temp.flush()

    @staticmethod
    def create_logger(name, simple=True):
        if simple is True:
            logging_queue = Queue()
            logging_queue.put_nowait = logging_queue.put
        else:
            logging_queue = mp.JoinableQueue(-1)
        log_queue_handler = logging.handlers.QueueHandler(logging_queue)
        log_queue_handler.setLevel(logging.DEBUG)

        logger = Mikado.utilities.log_utils.create_null_logger(name, level="WARNING")
        logger.propagate = False
        listener = logging.handlers.QueueListener(logging_queue, logger)
        listener.propagate = False
        listener.start()
        return logger, listener, logging_queue

    def setUp(self):
        # Create the queues for logging and submission
        self.submission_queue = mp.SimpleQueue()
        self.fasta_out = "temporary.fasta"
        self.gtf_out = "temporary.gtf"

    @unittest.skipUnless((version_info.minor > 5 and version_info.major == 3),
                         "SimpleQueues do not function properly for this in python 3.4 and earlier")
    def test_normal(self):
        logger, listener, logging_queue = self.create_logger("test_normal")

        with self.assertLogs(logger=logger, level="DEBUG") as cmo:
            # FASTA out and GTF out are just the file names, without the temporary directory
            # Moreover they will be complemented by the identifier!

            proc = ProcRunner(Mikado.preparation.checking.CheckingProcess,
                              self.submission_queue,
                              logging_queue,
                              fasta=self.fasta_temp.name,
                              identifier=0,
                              fasta_out=self.fasta_out,
                              gtf_out=self.gtf_out,
                              tmpdir=tempfile.gettempdir(),
                              log_level="DEBUG")
            proc.start()
            time.sleep(0.1)  # Necessary otherwise the check might be too fast for the FileSystem
            self.assertEqual(proc.func.fasta_out, os.path.join(tempfile.gettempdir(), self.fasta_out + "-0"))
            self.assertTrue(os.path.exists(proc.func.fasta_out), proc.func.fasta_out)
            self.assertEqual(proc.func.gtf_out, os.path.join(tempfile.gettempdir(), self.gtf_out + "-0"))
            self.assertTrue(os.path.exists(proc.func.gtf_out), proc.func.gtf_out)
            self.submission_queue.put(("EXIT", None, None, None))
            time.sleep(0.1)
            proc.stop()
            os.remove(proc.func.fasta_out)
            os.remove(proc.func.gtf_out)
            assert not proc.is_alive()

        self.maxDiff = 10000
        self.assertEqual(cmo.output, [
            "DEBUG:Checker-0:Starting Checker-0",
            "DEBUG:Checker-0:Created output FASTA {} and GTF {}".format(proc.func.fasta_out, proc.func.gtf_out),
            "DEBUG:Checker-0:(('GT', 'AG'), ('GC', 'AG'), ('AT', 'AC'))",
            "DEBUG:Checker-0:Finished for Checker-0"])

        self.assertIsInstance(proc.func, mp.Process)

        with self.assertRaises(TypeError):
            picked = pickle.dumps(proc.func)

    # def test_logger_creater(self):

    def test_wrong_initialisation(self):

        if version_info.minor < 6:
            self.submission_queue = mp.Queue(-1)
            logger, listener, logging_queue = self.create_logger("test_wrong_initialisation", simple=False)
        else:
            logger, listener, logging_queue = self.create_logger("test_wrong_initialisation", simple=True)

        kwds = {"submission_queue": self.submission_queue,
                "logging_queue": logging_queue,
                "identifier": 0,
                "fasta_out": self.fasta_out,
                "gtf_out": self.gtf_out,
                "tmpdir": tempfile.gettempdir(),
                "fasta": self.fasta_temp.name,
                "log_level": "WARNING"
                }

        for key in ["submission_queue", "logging_queue", "identifier", "lenient"]:
            with self.subTest(key=key):
                _kwds = kwds.copy()
                _kwds[key] = None
                if "queue" in key:
                    error = TypeError
                else:
                    error = ValueError
                with self.assertRaises(error):
                    Mikado.preparation.checking.CheckingProcess(**_kwds)

        with self.subTest(key="fasta"):
            _kwds = kwds.copy()
            _kwds["fasta"] = None


            with self.assertRaises(AttributeError):
                Mikado.preparation.checking.CheckingProcess(**_kwds)

        for tentative in [None, [], [("A", "G")], [("AG", bytes("GT", encoding="ascii"))]]:
            with self.subTest(tentative=tentative):
                _kwds = kwds.copy()
                _kwds["canonical_splices"] = tentative
                if tentative is None:
                    with self.assertRaises(TypeError):
                        Mikado.preparation.checking.CheckingProcess(**_kwds)
                else:
                    with self.assertRaises(ValueError):
                        Mikado.preparation.checking.CheckingProcess(**_kwds)

        _kwds = kwds.copy()
        _kwds["canonical_splices"] = [("AG", "GT")]
        # just test it does not raise
        _ = Mikado.preparation.checking.CheckingProcess(**_kwds)

    def test_example_model(self):

        fasta = pyfaidx.Fasta(self.fasta_temp.name)
        lines = dict()
        lines["chrom"] = "Chr5"
        lines["strand"] = "+"
        lines["start"] = 208937
        lines["end"] = 210445
        lines["attributes"] = dict()
        lines["tid"], lines["parent"] = "AT5G01530.0", "AT5G01530"
        lines["features"] = dict()
        lines["features"]["exon"] = [(208937, 209593), (209881, 210445)]
        seq = str(fasta[lines["chrom"]][lines["start"] - 1:lines["end"]])

        logger, listener, logging_queue = self.create_logger("test_example_model")

        res = Mikado.preparation.checking.create_transcript(lines, seq, lines["start"], lines["end"],
                                                            logger=logger)
        listener.stop()
        self.assertIsInstance(res, Mikado.transcripts.TranscriptChecker)

        for kwd in lines.keys():
            if kwd in ["end", "start"]:
                continue
            with self.subTest(kwd=kwd, msg="Testing key {}".format(kwd)):
                _lines = lines.copy()
                del _lines[kwd]
                with self.assertLogs("null", level="DEBUG"):
                    res = Mikado.preparation.checking.create_transcript(_lines, seq, lines["start"], lines["end"])
                self.assertIs(res, None)

        _lines = lines.copy()
        _lines["strand"] = "-"
        with self.subTest(msg="Testing an invalid strand"):
            with self.assertLogs("null", level="INFO") as cm:
                res = Mikado.preparation.checking.create_transcript(_lines, seq, lines["start"], lines["end"],
                                                                    strand_specific=True)
            self.assertIsInstance(res, Mikado.transcripts.TranscriptChecker)
            self.assertIn("WARNING:null:Transcript AT5G01530.0 has been assigned to the wrong strand, reversing it.",
                          cm.output)

    @unittest.skipUnless((version_info.minor > 5 and version_info.major == 3),
                         "SimpleQueues do not function properly for this in python 3.4 and earlier")
    def test_example_model_through_process(self):

        logger, listener, logging_queue = self.create_logger("test_example_model_through_process")
        logger.setLevel("DEBUG")
        with self.assertLogs(logger=logger, level="DEBUG") as cmo:
            # FASTA out and GTF out are just the file names, without the temporary directory
            # Moreover they will be complemented by the identifier!

            proc = ProcRunner(Mikado.preparation.checking.CheckingProcess,
                              self.submission_queue,
                              logging_queue,
                              fasta=self.fasta_temp.name,
                              identifier=logger.name,
                              fasta_out=self.fasta_out,
                              gtf_out=self.gtf_out,
                              tmpdir=tempfile.gettempdir(),
                              log_level="DEBUG")
            lines = dict()
            lines["chrom"] = "Chr5"
            lines["strand"] = "+"
            lines["start"] = 208937
            lines["end"] = 210445
            lines["attributes"] = dict()
            lines["tid"], lines["parent"] = "AT5G01530.0", "AT5G01530"
            lines["features"] = dict()
            lines["features"]["exon"] = [(208937, 209593), (209881, 210445)]
            lines["strand_specific"] = True
            lines["is_reference"] = False
            self.submission_queue.put((lines, lines["start"], lines["end"], 0))
            self.submission_queue.put(("EXIT", None, None, None))
            proc.start()
            proc.join()
            time.sleep(0.5)
            self.assertTrue(os.stat(proc.func.fasta_out).st_size > 0, proc.func.fasta_out)
            fasta_lines = []
            with open(proc.func.fasta_out) as f_out:
                for line in f_out:
                    line = line.rstrip()
                    line = re.sub("0/", "", line)
                    fasta_lines.append(line)

            self.assertGreater(len(fasta_lines), 1)
            fasta = pyfaidx.Fasta(self.fasta_temp.name)
            seq = str(fasta[lines["chrom"]][lines["start"] - 1:lines["end"]])
            res = Mikado.preparation.checking.create_transcript(lines, seq, lines["start"], lines["end"])
            self.assertTrue(len(res.cdna), (209593 - 208937 + 1) + (210445 - 209881 + 1))

            with tempfile.NamedTemporaryFile(suffix="fa", delete=True, mode="wt") as faix:

                assert len(fasta_lines) > 1
                assert fasta_lines[0][0] == ">"
                for line in fasta_lines:
                    print(line, file=faix)
                faix.flush()
                try:
                    fa = pyfaidx.Fasta(faix.name)
                except TypeError:
                    raise TypeError([_ for _ in open(faix.name)])
                self.assertEqual(list(fa.keys()), ["AT5G01530.0"])
                self.assertEqual(str(fa["AT5G01530.0"]), str(res.cdna))
                self.assertEqual(len(str(fa["AT5G01530.0"])), res.cdna_length)

                os.remove(faix.name + ".fai")
            fasta.close()

        listener.stop()
