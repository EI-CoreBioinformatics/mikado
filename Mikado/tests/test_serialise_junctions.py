#!/usr/bin/env python3

import logging
import os
import unittest
# import Mikado
from .. import utilities, configuration, parsers, serializers
import tempfile
import sqlalchemy.orm
from sqlalchemy import and_  # , or_
from pkg_resources import resource_stream
import gzip
import re
from pytest import mark


__author__ = 'Luca Venturini'


class TestLoadJunction(unittest.TestCase):

    logger = utilities.log_utils.create_null_logger("test_junction")

    def setUp(self):
        self.dbfile = tempfile.mktemp(suffix=".db")
        self.json_conf = configuration.configurator.to_json(None)
        self.json_conf["db_settings"]["dbtype"] = "sqlite"
        self.json_conf["db_settings"]["db"] = self.dbfile
        self.json_conf["reference"]["genome_fai"] = os.path.join(
            os.path.dirname(__file__),
            "genome.fai")
        self.session = utilities.dbutils.connect(self.json_conf)
        self.junction_file = os.path.join(
            os.path.dirname(__file__),
            "junctions.bed"
        )
        self.junction_serialiser = serializers.junction.JunctionSerializer(
            self.junction_file,
            json_conf=self.json_conf,
            # logger=self.logger
        )

        self.junction_parser = parsers.bed12.Bed12Parser(
            self.junction_file,
            fasta_index=None,
            transcriptomic=False
        )
        self.junction_serialiser()

    def test_first_junc(self):

        line = next(self.junction_parser)
        self.assertIsInstance(line, parsers.bed12.BED12)
        self.assertTrue(line.header)
        line = next(self.junction_parser)
        self.assertFalse(line.header)
        _line = ["Chr1",26510618,26511114, "portcullis_junc_0a",
                 39, "+", 26510751,	26510991, "255,0,0",
                 2,	"133,123", "0,373"]
        _line = "\t".join([str(_) for _ in _line])
        self.assertEqual(line._line, _line)

    def __create_session(self):
        engine = utilities.dbutils.connect(
            self.json_conf, self.logger)
        sessionmaker = sqlalchemy.orm.sessionmaker(bind=engine)
        session = sessionmaker()
        return session

    def test_serialise(self):

        session = self.__create_session()
        self.assertEqual(
            session.query(serializers.junction.Junction).count(),
            372,
            session.query(serializers.junction.Junction).count()
        )

        self.assertEqual(
            session.query(serializers.junction.Junction).filter(
                serializers.junction.Junction.chrom != "Chr5"
            ).count(), 1,
            [_.chrom for _ in
                session.query(serializers.junction.Junction).filter(
                    serializers.junction.Junction.chrom != "Chr5"
            )])

        # It's a BED file translated into 1-based, so add 1 to starts
        self.assertEqual(
            session.query(serializers.junction.Junction).filter(
                and_(
                    serializers.junction.Junction.chrom == "Chr5",
                    serializers.junction.Junction.start == 26510619,
                )
            ).count(), 1,
            [str(_) for _ in
                session.query(serializers.junction.Junction).filter(
                    and_(
                        serializers.junction.Junction.name == "portcullis_junc_0",
                    )
            )])

    def test_serialise_low_maxobject(self):

        self.dbfile = tempfile.mktemp(suffix=".db")
        self.json_conf = configuration.configurator.to_json(None)
        self.json_conf["db_settings"]["dbtype"] = "sqlite"
        self.json_conf["db_settings"]["db"] = self.dbfile
        self.json_conf["reference"]["genome_fai"] = os.path.join(
            os.path.dirname(__file__),
            "genome.fai")
        self.session = utilities.dbutils.connect(self.json_conf)
        self.junction_file = os.path.join(
            os.path.dirname(__file__),
            "junctions.bed"
        )

        self.json_conf["serialise"]["max_objects"] = 100

        self.logger.setLevel("DEBUG")
        self.junction_serialiser = serializers.junction.JunctionSerializer(
            self.junction_file,
            json_conf=self.json_conf,
            logger=self.logger,

        )

        self.junction_parser = parsers.bed12.Bed12Parser(
            self.junction_file,
            fasta_index=None,
            transcriptomic=False
        )

        with self.assertLogs(self.logger, level="DEBUG") as cmo:
            self.junction_serialiser()

        self.assertTrue(any(re.search("DEBUG:{}:Serializing [0-9][0-9][0-9] objects".format(self.logger.name), _)
                            for _ in cmo.output), cmo.output)

    def test_double_thick_end(self):

        session = self.__create_session()
        self.assertEqual(
            session.query(serializers.junction.Junction).filter(
                and_(
                    serializers.junction.Junction.chrom == "Chr5",
                    serializers.junction.Junction.junction_end == 26514549,
                )
            ).count(), 2,
            session.query(serializers.junction.Junction).filter(
                and_(
                    serializers.junction.Junction.chrom == "Chr5",
                    serializers.junction.Junction.junction_end == 26514549,
                )
            )
        )

        first = session.query(serializers.junction.Junction).filter(
                    and_(
                        serializers.junction.Junction.name == "portcullis_junc_0",
                    )).one()
        first_double = session.query(serializers.junction.Junction).filter(
                        and_(
                            serializers.junction.Junction.name == "portcullis_junc_0",
                        )).one()
        second = session.query(serializers.junction.Junction).filter(
                    and_(
                        serializers.junction.Junction.name == "portcullis_junc_1",
                    )).one()

        self.assertTrue(first.is_equal(first_double.chrom,
                                       first_double.start,
                                       first_double.end,
                                       first_double.strand))
        self.assertFalse(first.is_equal(second.chrom,
                                        second.start,
                                        second.end,
                                        second.strand))

    def test_no_logger(self):

        with serializers.junction.JunctionSerializer(
                self.junction_file,
                json_conf=self.json_conf,
                logger=None
                ) as jser:
            self.assertIsInstance(jser.logger, logging.Logger)
            self.assertEqual(jser.logger.name, "junction")

    def test_exiting(self):

        nlogger = logging.getLogger("test_exiting")
        nlogger.setLevel(logging.INFO)
        handler = logging.NullHandler()
        nlogger.addHandler(handler)

        with self.assertLogs("test_exiting", level="WARNING") as cm:
            _ = serializers.junction.JunctionSerializer(
                None,
                json_conf=self.json_conf,
                logger=nlogger
                )
            _()

        self.assertEqual(cm.output,
                         ['WARNING:test_exiting:No input file specified. Exiting.'] * 2,
                         cm.output)

    @mark.slow
    def test_no_fai(self):

        db = tempfile.mktemp(suffix=".db")
        genome_file = tempfile.NamedTemporaryFile("wb", suffix=".fa.gz", prefix="Chr5", dir=".")
        jconf = self.json_conf.copy()
        jconf["db_settings"]["db"] = db
        jconf["reference"]["genome_fai"] = None
        with resource_stream("Mikado.tests", "chr5.fas.gz") as _:
            genome_file.write(_.read())
        genome_file.flush()

        jconf["reference"]["genome"] = genome_file.name

        seri = serializers.junction.JunctionSerializer(
                self.junction_file,
                json_conf=self.json_conf,
                logger=self.logger
                )
        seri()
        genome_file.close()
        os.remove("{}.fai".format(genome_file.name))
        # genome_file.delete()

    def test_invalid_bed12(self):

        with self.assertRaises(TypeError):
            _ = serializers.junction.Junction(None, 0)

    def tearDown(self):
        self.junction_serialiser.close()
        self.junction_parser.close()
        if os.path.exists(self.dbfile):
            os.remove(self.dbfile)


if __name__ == "__main__":
    unittest.main()