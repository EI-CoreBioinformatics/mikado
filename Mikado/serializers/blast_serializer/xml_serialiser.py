import rapidjson as json
from ...parsers.blast_utils import BlastOpener
from .xml_utils import _get_query_for_blast, _get_target_for_blast
from xml.parsers.expat import ExpatError
import xml
from ...utilities.dbutils import connect
from sqlalchemy.orm.session import Session
import msgpack
import tempfile
import os
from . import Query, Target, prepare_hit, InvalidHit
from .xml_utils import get_multipliers, get_off_by_one
from .utils import load_into_db
import multiprocessing as mp
from ...utilities.log_utils import create_null_logger
import pandas as pd
import zlib
import struct


struct_row = struct.Struct(">LLL")


def xml_pickler(json_conf, filename, default_header,
                cache=None,
                max_target_seqs=10):
    valid, _, exc = BlastOpener(filename).sniff(default_header=default_header)
    engine = connect(json_conf, strategy="threadlocal")
    session = Session(bind=engine)

    if not valid:
        err = "Invalid BLAST file: %s" % filename
        raise TypeError(err)
    directory = os.path.dirname(filename)
    try:
        dbname = tempfile.mktemp(suffix=".db", dir=directory)
        dbhandle = open(dbname, "wb")
    except (OSError, PermissionError):
        dbname = tempfile.mktemp(suffix=".db")
        dbhandle = open(dbname, "wb")

    if not isinstance(cache, dict) or set(cache.keys()) != {"query", "target"}:
        cache = dict()
        cache["query"] = dict((item.query_name, item.query_id) for item in session.query(Query))
        cache["target"] = dict((item.target_name, item.target_id) for item in session.query(Target))

    rows = b""
    try:
        with BlastOpener(filename) as opened:
            try:
                qmult, tmult = None, None
                for query_counter, record in enumerate(opened, start=1):
                    if qmult is None:
                        qmult, tmult = get_multipliers(record)
                        off_by_one = get_off_by_one(record)
                    hits, hsps, cache = objectify_record(
                        record, [], [], cache, max_target_seqs=max_target_seqs,
                        qmult=qmult, tmult=tmult, off_by_one=off_by_one)

                    record = {"hits": hits, "hsps": hsps}
                    try:
                        jrecord = json.dumps(record, number_mode=json.NM_NATIVE)
                    except ValueError:
                        jrecord = json.dumps(record)
                    write_start = dbhandle.tell()
                    write_length = dbhandle.write(msgpack.dumps(jrecord))
                    rows += struct_row.pack(query_counter, write_start, write_length)

            except ExpatError as err:
                raise ExpatError("{} is an invalid BLAST file, sending back anything salvageable.\n{}".format(filename,
                                                                                                              err))
    except xml.etree.ElementTree.ParseError as err:
        # logger.error("%s is an invalid BLAST file, sending back anything salvageable", filename)
        raise xml.etree.ElementTree.ParseError(
            "{} is an invalid BLAST file, sending back anything salvageable.\n{}".format(filename, err))
    except ValueError as err:
        # logger.error("Invalid BLAST entry")
        raise ValueError(
            "{} is an invalid BLAST file, sending back anything salvageable.\n{}".format(filename, err))

    return dbname, zlib.compress(rows)


def _serialise_xmls(self):
    # Load sequences in DB, precache IDs

    if isinstance(self.xml, str):
        self.xml = [self.xml]
    else:
        assert isinstance(self.xml, (list, set))

    hits, hsps = [], []
    hit_counter, record_counter = 0, 0

    for filename in self.xml:
        valid, header, exc = BlastOpener(filename).sniff()
        if valid is True:
            self.header = header
        else:
            self.logger.error(exc)
            self.xml.remove(filename)

    cache = {"query": self.queries, "target": self.targets}
    if self._xml_debug is False and (self.single_thread is True or self.procs == 1):
        for filename in self.xml:
            valid, _, exc = BlastOpener(filename).sniff(default_header=self.header)
            if not valid:
                self.logger.error(exc)
                continue
            try:
                self.logger.debug("Analysing %s", filename)
                qmult, tmult = None, None
                with BlastOpener(filename) as opened:
                    for record in opened:
                        if qmult is None:
                            qmult, tmult = get_multipliers(record)
                            off_by_one = get_off_by_one(record)
                        record_counter += 1
                        if record_counter > 0 and record_counter % 10000 == 0:
                            self.logger.info("Parsed %d queries", record_counter)
                        current = len(hits)
                        hits, hsps, cache = objectify_record(
                            record, hits, hsps,
                            qmult=qmult,
                            tmult=tmult,
                            cache=cache,
                            max_target_seqs=self._max_target_seqs, logger=self.logger, off_by_one=off_by_one)
                        hit_counter += len(hits) - current
                        hits, hsps = load_into_db(self, hits, hsps, force=False)
                self.logger.debug("Finished %s", filename)
            except ExpatError:
                self.logger.error("%s is an invalid BLAST file, saving what's available", filename)
        _, _ = load_into_db(self, hits, hsps, force=True)
    elif self._xml_debug is True or self.procs > 1:
        self.logger.debug("Creating a pool with %d processes",
                          min(self.procs, len(self.xml)))
        results = []
        if self._xml_debug is True:
            for num, filename in enumerate(self.xml):
                results.append(xml_pickler(self.json_conf,
                                           filename, self.header,
                                           cache=None,
                                           max_target_seqs=self._max_target_seqs))
        else:
            pool = mp.Pool(self.procs)
            for num, filename in enumerate(self.xml):
                args = (self.json_conf, filename, self.header)
                kwds = {"max_target_seqs": self._max_target_seqs, "cache": None}
                pool.apply_async(xml_pickler, args=args, kwds=kwds, callback=results.append)
            pool.close()
            pool.join()

        for dbfile, rows in results:
            with open(dbfile, "rb") as dbhandle:
                for query_counter, write_start, write_length in struct_row.iter_unpack(zlib.decompress(rows)):
                    dbhandle.seek(write_start)
                    result = json.loads(msgpack.loads(dbhandle.read(write_length), raw=False, use_list=True))
                    __hits = result["hits"]
                    __hsps = result["hsps"]
                    assert (not __hsps and not __hits) or __hsps != __hits, (__hits, __hsps)
                    hit_counter += len(__hits)
                    hits.extend(__hits)
                    hsps.extend(__hsps)
                    hits, hsps = load_into_db(self, hits, hsps, force=False)
                    if record_counter > 0 and record_counter % 10000 == 0:
                        self.logger.debug("Parsed %d queries", record_counter)
            os.remove(dbfile)

        self.logger.debug("Finished sending off the data for serialisation")
        _, _ = load_into_db(self, hits, hsps, force=True)

    self.logger.info("Loaded %d alignments for %d queries",
                     hit_counter, record_counter)

    self.logger.info("Finished loading blast hits")
    if hasattr(self, "logging_queue"):
        self.logging_queue.close()


def objectify_record(record, hits, hsps, cache,
                     max_target_seqs=10000, logger=create_null_logger(),
                     qmult=1, tmult=1, off_by_one=False):
    """
    Private method to serialise a single record into the DB.

    :param record: The BLAST record to load into the DB.
    :type record: Bio.SearchIO.QueryResult
    :param hits: Cache of hits to load into the DB.
    :type hits: list

    :param hsps: Cache of hsps to load into the DB.
    :type hsps: list

    :param cache: cache of match query_name/query_id
    :type cache: dict

    :param logger: logger
    :type logger: logging.Logger

    :param max_target_seqs: Maximum target sequences per query
    :type max_target_seqs: int

    :param qmult: query multiplier (3 for BLASTX or TBLASTX, 1 otherwise)
    :type qmult: int

    :param qmult: query multiplier (1 for TBLASTX or TBLASTN, 1 otherwise)
    :type qmult: int

    :param off_by_one: DIAMOND before 0.9.31 has a bug which causes the query end to be off by 1.
    :type off_by_one: bool

    :returns: hits, hsps, cache
    :rtype: (list, list, dict)
    """

    if len(record.hits) == 0:
        return hits, hsps, cache

    current_query, name, cache["query"] = _get_query_for_blast(record, cache["query"])
    alignments = dict((ccc, (hit.id, max(_.bitscore for _ in hit.hsps))) for ccc, hit in enumerate(record.hits))
    ids, scores = list(zip(*alignments.values()))
    phits = pd.DataFrame({"idx": list(alignments.keys()), "target_name": ids, "bitscore": scores})
    for hit_num, t in enumerate(phits.sort_values(["bitscore", "target_name"],
                                                  ascending=[False, True]).itertuples(name=None)):
        idx, target_name, bitscore = t[1:]
        alignment = record.hits[idx]
        logger.debug("Started the hit %s vs. %s", name, record.hits[idx].id)
        current_target, cache["target"] = _get_target_for_blast(alignment, cache["target"])
        hit_dict_params = dict()
        (hit_dict_params["query_multiplier"],
         hit_dict_params["target_multiplier"]) = (qmult, tmult)
        hit_evalue = min(_.evalue for _ in record.hits[idx].hsps)
        hit_bs = bitscore
        hit_dict_params["hit_number"] = hit_num + 1
        hit_dict_params["evalue"] = hit_evalue
        hit_dict_params["bits"] = hit_bs
        try:
            hit, hit_hsps = prepare_hit(alignment, current_query,
                                        current_target,
                                        query_length=record.seq_len,
                                        off_by_one=off_by_one,
                                        **hit_dict_params)
        except InvalidHit as exc:
            logger.error(exc)
            continue
        hit["query_aligned_length"] = min(record.seq_len, hit["query_aligned_length"])
        hits.append(hit)
        hsps.extend(hit_hsps)

    return hits, hsps, cache
