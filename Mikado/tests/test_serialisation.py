#!/usr/bin/env python3

import logging
import os
import unittest
import Mikado
import tempfile
import sqlalchemy.orm
from sqlalchemy import and_  # , or_
from pkg_resources import resource_stream
import gzip

__author__ = 'Luca Venturini'


class TestLoadJunction(unittest.TestCase):

    logger = Mikado.utilities.log_utils.create_null_logger("test_junction")
    dbfile = tempfile.mktemp(suffix=".db")

    def setUp(self):

        self.json_conf = Mikado.configuration.configurator.to_json(None)
        self.json_conf["db_settings"]["dbtype"] = "sqlite"
        self.json_conf["db_settings"]["db"] = self.dbfile
        self.json_conf["reference"]["genome_fai"] = os.path.join(
            os.path.dirname(__file__),
            "genome.fai")
        self.session = Mikado.utilities.dbutils.connect(self.json_conf)
        self.junction_file = os.path.join(
            os.path.dirname(__file__),
            "junctions.bed"
        )
        self.junction_serialiser = Mikado.serializers.junction.JunctionSerializer(
            self.junction_file,
            json_conf=self.json_conf,
            logger=self.logger
        )

        self.junction_parser = Mikado.parsers.bed12.Bed12Parser(
            self.junction_file,
            fasta_index=None,
            transcriptomic=False
        )
        self.junction_serialiser()

    def test_first_junc(self):

        line = next(self.junction_parser)
        self.assertIsInstance(line, Mikado.parsers.bed12.BED12)
        self.assertTrue(line.header)
        line = next(self.junction_parser)
        self.assertFalse(line.header)
        _line = ["Chr1",26510618,26511114, "portcullis_junc_0a",
                 39, "+", 26510751,	26510991, "255,0,0",
                 2,	"133,123", "0,373"]
        _line = "\t".join([str(_) for _ in _line])
        self.assertEqual(line._line, _line)

    def __create_session(self):
        engine = Mikado.utilities.dbutils.connect(
            self.json_conf, self.logger)
        sessionmaker = sqlalchemy.orm.sessionmaker(bind=engine)
        session = sessionmaker()
        return session

    def test_serialise(self):

        session = self.__create_session()
        self.assertEqual(
            session.query(Mikado.serializers.junction.Junction).count(),
            372,
            session.query(Mikado.serializers.junction.Junction).count()
        )

        self.assertEqual(
            session.query(Mikado.serializers.junction.Junction).filter(
                Mikado.serializers.junction.Junction.chrom != "Chr5"
            ).count(), 1,
            [_.chrom for _ in
                session.query(Mikado.serializers.junction.Junction).filter(
                    Mikado.serializers.junction.Junction.chrom != "Chr5"
            )])

        # It's a BED file translated into 1-based, so add 1 to starts
        self.assertEqual(
            session.query(Mikado.serializers.junction.Junction).filter(
                and_(
                    Mikado.serializers.junction.Junction.chrom == "Chr5",
                    Mikado.serializers.junction.Junction.start == 26510619,
                )
            ).count(), 1,
            [str(_) for _ in
                session.query(Mikado.serializers.junction.Junction).filter(
                    and_(
                        Mikado.serializers.junction.Junction.name == "portcullis_junc_0",
                    )
            )])

    def test_double_thick_end(self):

        session = self.__create_session()
        self.assertEqual(
            session.query(Mikado.serializers.junction.Junction).filter(
                and_(
                    Mikado.serializers.junction.Junction.chrom == "Chr5",
                    Mikado.serializers.junction.Junction.junction_end == 26514549,
                )
            ).count(), 2,
            session.query(Mikado.serializers.junction.Junction).filter(
                and_(
                    Mikado.serializers.junction.Junction.chrom == "Chr5",
                    Mikado.serializers.junction.Junction.junction_end == 26514549,
                )
            )
        )

        first = session.query(Mikado.serializers.junction.Junction).filter(
                    and_(
                        Mikado.serializers.junction.Junction.name == "portcullis_junc_0",
                    )).one()
        first_double = session.query(Mikado.serializers.junction.Junction).filter(
                        and_(
                            Mikado.serializers.junction.Junction.name == "portcullis_junc_0",
                        )).one()
        second = session.query(Mikado.serializers.junction.Junction).filter(
                    and_(
                        Mikado.serializers.junction.Junction.name == "portcullis_junc_1",
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

        with Mikado.serializers.junction.JunctionSerializer(
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
            _ = Mikado.serializers.junction.JunctionSerializer(
                None,
                json_conf=self.json_conf,
                logger=nlogger
                )
            _()

        self.assertEqual(cm.output,
                         ['WARNING:test_exiting:No input file specified. Exiting.'] * 2,
                         cm.output)

        # with Mikado.py.utilities.log_utils.create_default_logger("test_null") as new_logger:
        #
        #     with  as _:
        #         self.assertLogs(new_logger, "No input file specified. Exiting.")

    def test_no_fai(self):

        db = tempfile.mktemp(suffix=".db")
        genome_file = tempfile.NamedTemporaryFile("wb", suffix=".fa", prefix="Chr5", dir=".")
        jconf = self.json_conf.copy()
        jconf["db_settings"]["db"] = db
        jconf["reference"]["genome_fai"] = None
        with resource_stream("Mikado.tests", "chr5.fas.gz") as _:
            genome_file.write(gzip.decompress(_.read()))
        genome_file.flush()

        jconf["reference"]["genome"] = genome_file.name

        seri = Mikado.serializers.junction.JunctionSerializer(
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
            _ = Mikado.serializers.junction.Junction(None, 0)

    def tearDown(self):
        self.junction_serialiser.close()
        self.junction_parser.close()
        if os.path.exists(self.dbfile):
            os.remove(self.dbfile)


if __name__ == "__main__":
    unittest.main()