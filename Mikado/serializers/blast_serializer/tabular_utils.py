from Bio.Align import substitution_matrices
from functools import partial
from .btop_parser import parse_btop
import re
import numpy as np
import pandas as pd
import multiprocessing as mp
from .utils import load_into_db
from collections import defaultdict
import logging
import logging.handlers
from ...utilities.log_utils import create_null_logger, create_queue_logger
from sqlalchemy.orm.session import Session
from ...utilities.dbutils import connect as db_connect
from . import Hit, Hsp
import os
import tempfile
import msgpack
from ...utilities import blast_keys
import shutil
import subprocess as sp


hit_cols = [col.name for col in Hit.__table__.columns]
hsp_cols = [col.name for col in Hsp.__table__.columns]

__author__ = 'Luca Venturini'


class Matrices:

    def __init__(self):

        mnames = substitution_matrices.load()
        self.mnames = dict()

        for mname in mnames:
            self.mnames[mname] = mname
            self.mnames[mname.upper()] = mname
            self.mnames[mname.lower()] = mname

        self.__matrices = dict()

    def __contains__(self, mname):
        return mname in self.mnames

    def keys(self):
        return self.mnames.keys()

    def get(self, key, default=None):
        if key not in self and default is not None:
            return default
        return self.__getitem__(key)

    def __getitem__(self, mname):

        if mname not in self:
            raise KeyError(f"{mname} is not a valid matrix name. Valid names: {self.keys()}")
        if mname not in self.__matrices:
            self.__load_matrix(mname)
        return self.__matrices[mname.lower()]

    def __load_matrix(self, mname):
        matrix = dict()
        orig_mname = self.mnames[mname]
        omatrix = substitution_matrices.load(orig_mname)
        for key, val in omatrix.items():
            if key[::-1] in omatrix and omatrix[key[::-1]] != val:
                raise KeyError((key, val, key[::-1], omatrix[key[::-1]]))
            matrix["".join(key)] = val
            matrix["".join(key[::-1])] = val
        for key in orig_mname, orig_mname.lower(), orig_mname.upper():
            self.__matrices[key] = matrix


matrices = Matrices()


def get_queries(engine):
    queries = pd.read_sql_table("query", engine, index_col="query_name")
    queries.columns = ["qid", "qlength"]
    queries["qid"] = queries["qid"].astype(int)
    queries.index = queries.index.str.replace(id_pattern, "\\1")
    assert queries.qid.drop_duplicates().shape[0] == queries.shape[0]
    return queries


def get_targets(engine):
    targets = pd.read_sql_table("target", engine, index_col="target_name")
    targets.columns = ["sid", "slength"]
    targets["sid"] = targets["sid"].astype(int)
    assert targets.sid.drop_duplicates().shape[0] == targets.shape[0]

    targets.index = targets.index.str.replace(id_pattern, "\\1")
    if targets[targets.slength.isna()].shape[0] > 0:
        raise KeyError("Unbound targets!")
    return targets


def prepare_tab_hsp(key,
                    hsp: pd.Series,
                    columns: dict,
                    qmult=3, tmult=1, matrix_name=None):

    r"""
    Prepare a HSP for loading into the DB.
    The match line will be reworked in the following way:

    - If the position is a match/positive, keep the original value
    - If the position is a gap *for the query*, insert a - (dash)
    - If the position is a gap *for the target*, insert a _ (underscore)
    - If the position is a gap *for both*, insert a \ (backslash)

    :param hsp: A tabular blast row
    :type hsp: tuple
    :param columns: dictionary with the index of the column names
    :return: hsp_dict, numpy array
    :rtype: (tuple, dict, np.array, np.array)
    """

    hsp_dict = dict()
    # We must start from 1, otherwise MySQL crashes as its indices start from 1 not 0
    hsp_dict["query_id"], hsp_dict["target_id"] = key
    try:
        query_array = np.zeros([3, int(hsp[columns["qlength"]])], dtype=np.int)
    except (IndexError,TypeError):
        try:
            raise IndexError((hsp[columns["qlength"]], type(hsp[columns["qlength"]])))
        except IndexError:
            raise IndexError(columns)

    target_array = np.zeros([3, int(hsp[columns["slength"]])], dtype=np.int)
    matrix = matrices.get(matrix_name, matrices["blosum62"])
    if hsp[columns["qstart"]] < 0:
        raise ValueError(hsp.qstart)
    try:
        query_array, target_array, aln_span, match = parse_btop(
        hsp[columns["btop"]],
        query_array=query_array,
        target_array=target_array,
        qpos=int(hsp[columns["qstart"]]),
        spos=int(hsp[columns["sstart"]]),
        qmult=qmult,
        tmult=tmult,
        matrix=matrix)
    except TypeError:
        raise TypeError((hsp[columns["btop"]], hsp[columns["qstart"]],
                         hsp[columns["sstart"]], qmult, tmult))
    hsp_dict["counter"] = hsp[columns["hsp_num"]]
    hsp_dict["query_hsp_start"] = hsp[columns["qstart"]]
    hsp_dict["query_hsp_end"] = hsp[columns["qend"]]
    hsp_dict["query_frame"] = hsp[columns["query_frame"]]
    hsp_dict["target_hsp_start"] = hsp[columns["sstart"]]
    hsp_dict["target_hsp_end"] = hsp[columns["send"]]
    hsp_dict["target_frame"] = hsp[columns["target_frame"]]
    span = np.where(query_array[0] > 0)[0]
    if span.shape[0] == 0:
        raise ValueError((hsp[columns["btop"]], type(hsp[columns["btop"]])))
    if aln_span != hsp[columns["length"]]:
        raise ValueError((aln_span, hsp[columns["length"]]))
    pident = np.where(query_array[1] > 0)[0].shape[0] / (aln_span * qmult) * 100
    # if not np.isclose(pident, hsp[columns["pident"]], atol=.1, rtol=.1):
    #     raise ValueError((pident, hsp[columns["pident"]]))
    hsp_dict["hsp_identity"] = pident
    ppos = np.where(query_array[2] > 0)[0].shape[0] / (aln_span * qmult) * 100
    hsp_dict["hsp_positives"] = ppos
    hsp_dict["match"] = match
    hsp_dict["hsp_length"] = hsp[columns["length"]]
    hsp_dict["hsp_bits"] = hsp[columns["bitscore"]]
    hsp_dict["hsp_evalue"] = hsp[columns["evalue"]]
    return key, hsp_dict, query_array, target_array


def prepare_tab_hit(key: tuple,
                    hit: list,
                    columns: dict,
                    matrix_name: str,
                    qmult=3, tmult=1, **kwargs):
    """
    :param hit:
    :param columns:
    :param matrix_name:
    :param qmult:
    :param tmult:
    :param kwargs:
    :return:
    :rtype: (dict, list[dict])
    """

    hit_dict = dict()
    query_arrays = []
    target_arrays = []
    hsps = []
    hit_row = None
    key = tuple([int(key[0]), int(key[1])])
    for hsp in hit:
        if hit_row is None:
            hit_row = hsp
        key, hsp_dict, query_array, target_array = prepare_tab_hsp(key, hsp,
                                                                   columns=columns,
                                                                   qmult=qmult, tmult=tmult,
                                                                   matrix_name=matrix_name)
        query_arrays.append(query_array)
        target_arrays.append(target_array)
        hsps.append(hsp_dict)
    
    qlength = int(hit_row[columns["qlength"]])
    hit_dict.update(kwargs)
    hit_dict["query_id"] = key[0]
    hit_dict["target_id"] = key[1]
    hit_dict["hit_number"] = int(hit_row[columns["hit_num"]])
    hit_dict["evalue"] = hit_row[columns["min_evalue"]]
    hit_dict["bits"] = hit_row[columns["max_bitscore"]]
    hit_dict["query_multiplier"] = int(qmult)
    hit_dict["target_multiplier"] = int(tmult)

    query_array = sum(query_arrays)
    target_array = sum(target_arrays)
    q_aligned = np.where(query_array[0] > 0)[0]
    hit_dict["query_aligned_length"] = min(qlength, q_aligned.shape[0])
    qstart, qend = q_aligned.min(), q_aligned.max() + 1
    hit_dict["query_start"], hit_dict["query_end"] = int(qstart), int(qend)
    ends, sends = list(zip(*[(int(hsp[columns["qend"]]), int(hsp[columns["send"]])) for hsp in hit]))
    if int(qend) not in ends:
        raise ValueError("Invalid end point: {}, {}".format(qend, ends))
    identical_positions = np.where(query_array[1] > 0)[0].shape[0]
    positives = np.where(query_array[2] > 0)[0].shape[0]
    if identical_positions > q_aligned.shape[0]:
        raise ValueError(
            "Number of identical positions ({}) greater than number of aligned positions ({})!".format(
                len(identical_positions), q_aligned))

    if positives > q_aligned.shape[0]:
        raise ValueError("Number of identical positions ({}) greater than number of aligned positions ({})!".format(
            len(positives), q_aligned))

    t_aligned = np.where(target_array[0] > 0)[0]
    hit_dict["target_aligned_length"] = min(t_aligned.shape[0], hit_row[columns["slength"]])
    hit_dict["target_start"] = int(t_aligned.min())
    hit_dict["target_end"] = int(t_aligned.max() + 1)
    if hit_dict["target_end"] not in sends:
        raise ValueError("Invalid target end point: {}, {}".format(hit_dict["target_end"], sends))
    hit_dict["global_identity"] = identical_positions * 100 / q_aligned.shape[0]
    hit_dict["global_positives"] = positives * 100 / q_aligned.shape[0]
    hit_list = [hit_dict[col] for col in hit_cols]
    hsps = [[hsp[col] for col in hsp_cols] for hsp in hsps]

    return hit_list, hsps


id_pattern = re.compile(r"^[^\|]*\|([^\|]*)\|.*")


def sanitize_blast_data(data: pd.DataFrame, queries: pd.DataFrame, targets: pd.DataFrame,
                        qmult=3, tmult=1):

    if data[data.btop.isna()].shape[0] > 0:
        raise ValueError(data.loc[0])

    data["qseqid"] = data["qseqid"].str.replace(id_pattern, "\\1")
    data["sseqid"] = data["sseqid"].str.replace(id_pattern, "\\1")

    data = data.join(queries, on=["qseqid"]).join(targets, on=["sseqid"]).join(
        data.groupby(["qseqid", "sseqid"]).agg(
            min_evalue=pd.NamedAgg("evalue", np.min),
            max_bitscore=pd.NamedAgg("bitscore", np.max)
        )[["min_evalue", "max_bitscore"]], on=["qseqid", "sseqid"])

    for col in ["qstart", "qend", "sstart", "send", "qlength", "slength"]:
        if col != "slength":
            err_val = (col, data[data[col].isna()].shape[0], data.shape[0])
        else:
            err_val = (col, data[["sseqid"]].head())
        assert ~(data[col].isna().any()), err_val
        try:
            data[col] = data[col].astype(int).values
        except ValueError as exc:
            raise ValueError("{}: {}".format(exc, col))

    for key, multiplier, (start, end), length in [
        ("query_frame", qmult, ("qstart", "qend"), "qlength"),
        ("target_frame", tmult, ("sstart", "send"), "slength")]:
        # Switch start and end when they are not in the correct order
        _ix = (data[start] > data[end])
        if multiplier > 1:
            data.loc[~_ix, key] = data[start] % multiplier
            data.loc[_ix, key] = -((data[length] - data[end] - 1) % multiplier)
            data.loc[(data[key] == 0) & ~_ix, key] = multiplier
            data.loc[(data[key] == 0) & _ix, key] = -multiplier
        else:
            data.loc[:, key] = 0
        data.loc[_ix, [start, end]] = data.loc[_ix, [end, start]].values
        data[start] -= 1
    # Get the minimum evalue for each group
    # data["aln_span"] = data.qend - data.qstart
    # Set the hsp_num
    data["sstart"] = data["sstart"].astype(int).values
    data["hsp_num"] = data.sort_values("bitscore", ascending=False).groupby(["qseqid", "sseqid"]).cumcount() + 1
    temp = data[["qseqid", "sseqid", "max_bitscore"]].drop_duplicates().sort_values(
        ["max_bitscore", "sseqid"], ascending=[False, True])
    temp["hit_num"] = temp.groupby(["qseqid"]).cumcount() + 1
    temp.set_index(["qseqid", "sseqid"], inplace=True)
    data = data.join(temp["hit_num"], on=["qseqid", "sseqid"])
    data = data.sort_values(["qid", "sid"])
    data.set_index(["qid", "sid"], drop=False, inplace=True)
    return data


class Preparer(mp.Process):

    def __init__(self,
                 index_file: str,
                 identifier: int,
                 params_file: str,
                 lock: mp.RLock,
                 conf: dict,
                 maxobjects: int,
                 logging_queue=None,
                 log_level="DEBUG",
                 sql_level="DEBUG",
                 matrix_name=None, qmult=3, tmult=1, **kwargs):

        super().__init__()
        self.matrix_name = matrix_name
        self.qmult, self.tmult = qmult, tmult
        self.identifier = identifier
        self.log_level, self.sql_level = log_level, sql_level
        self.lock = lock
        self.maxobjects = maxobjects
        self.conf = conf
        self.params_file, self.index_file = params_file, index_file
        if logging_queue is None:
            self.logger = create_null_logger("preparer-{}".format(self.identifier))  # create_null_logger
            self.logging_queue = None
        else:
            self.logging_queue = logging_queue

    def run(self):
        with open(self.params_file, "rb") as pfile:
            params = msgpack.loads(pfile.read(), raw=False, strict_map_key=False)
        self.columns = params["columns"]
        prep_hit = partial(prepare_tab_hit,
                           columns=self.columns, qmult=self.qmult, tmult=self.tmult,
                           matrix_name=self.matrix_name)
        if self.logging_queue is not None:
            create_queue_logger(self)
            sql_logger = logging.getLogger("sqlalchemy.engine")
            sql_logger.setLevel(self.sql_level)
            sql_logger.addHandler(self.logger.handlers[0])
        else:
            sql_logger = None
        self.engine = db_connect(self.conf, logger=sql_logger)
        self.logger.debug("Started %s", self.identifier)
        session = Session(bind=self.engine)
        self.session = session
        hits, hsps = [], []
        with open(self.index_file, "rb") as index_handle:
            for key, rows in msgpack.Unpacker(index_handle, raw=False, strict_map_key=False):
                curr_hit, curr_hsps = prep_hit(key, rows)
                hits.append(curr_hit)
                hsps += curr_hsps
                hits, hsps = load_into_db(self, hits, hsps, force=False, raw=True)
        _, _ = load_into_db(self, hits, hsps, force=True, raw=True)
        self.logger.debug("Finished %s", self.identifier)
        os.remove(self.index_file)  # Clean it up
        return True


def parse_tab_blast(self,
                    bname: str,
                    queries: pd.DataFrame,
                    targets: pd.DataFrame,
                    procs: int,
                    matrix_name="blosum62", qmult=3, tmult=1):
    """This function will use `pandas` to quickly parse, subset and analyse tabular BLAST files.
    """

    matrix_name = matrix_name.lower()

    if matrix_name not in matrices:
        raise KeyError("Matrix {} is not valid. Please specify a valid name.".format(matrix_name))

    if procs > 1:
        # We have to set up the processes before the forking.
        lock = mp.RLock()
        conf = dict()
        conf["db_settings"] = self.json_conf["db_settings"].copy()
        params_file = tempfile.mktemp(suffix=".mgp")

        index_files = dict((idx, tempfile.mktemp(suffix=".csv")) for idx in
                           range(procs))
        kwargs = {"conf": conf,
                  "maxobjects": max(int(self.maxobjects / procs), 1),
                  "lock": lock,
                  "matrix_name": matrix_name,
                  "qmult": qmult,
                  "tmult": tmult,
                  "sql_level": self.json_conf["log_settings"]["sql_level"],
                  "log_level": self.json_conf["log_settings"]["log_level"],
                  "logging_queue": self.logging_queue,
                  "params_file": params_file}
        processes = [Preparer(index_files[idx], idx, **kwargs) for idx in range(procs)]

    self.logger.info("Reading %s data", bname)
    # Compatibility with ASN files
    if isinstance(bname, str) and bname.endswith(("asn", "asn.gz", "asn.bz2")):
        bash = shutil.which("bash")
        if bname.endswith("asn"):
            catter = "cat"
        elif bname.endswith(".gz"):
            catter = "gunzip -c"
        else:
            catter = "bunzip2 -c"

        if bash is not None:
            formatter = sp.Popen('blast_formatter -outfmt "6 {blast_keys}" -archive <({catter} {bname})"'.format(
                blast_keys=blast_keys, catter=catter, bname=bname), shell=True, executable=bash, stdout=sp.STDOUT)
        else:
            if catter != "cat":
                temp = tempfile.NamedTemporaryFile(mode="wt", suffix=".asn")
                sp.call("{catter} {bname}".format(bname=bname, catter=catter), shell=True, stdout=temp)
                temp.flush()
            formatter = sp.Popen('blast_formatter -outfmt "6 {blast_keys}" -archive {temp}'.format(
                blast_keys=blast_keys, catter=catter, bname=bname), shell=True, stdout=sp.STDOUT)
        data = pd.read_csv(formatter.stdout, delimiter="\t", names=blast_keys)
    elif isinstance(bname, str) and bname.endswith("daa"):
        formatter = sp.Popen("diamond view -a {bname} --outfmt 6 {blast_keys}".format(
            bname=bname, blast_keys=blast_keys), shell=True, stdout=sp.STDOUT)
        data = pd.read_csv(formatter.stdout, delimiter="\t", names=blast_keys)
    else:
        data = pd.read_csv(bname, delimiter="\t", names=blast_keys)

    data = sanitize_blast_data(data, queries, targets, qmult=qmult, tmult=tmult)
    columns = dict((col, idx) for idx, col in enumerate(data.columns))
    groups = defaultdict(list)
    [groups[val].append(idx) for idx, val in enumerate(data.index)]
    values = data.values

    if procs == 1:
        prep_hit = partial(prepare_tab_hit, columns=columns, qmult=qmult, tmult=tmult,
                           matrix_name=matrix_name)
        hits, hsps = [], []
        self.logger.info("Finished reading %s data, starting serialisation in single-threaded mode", bname)
        for key, group in groups.items():
            curr_hit, curr_hsps = prep_hit(key, values[group, :])
            hits.append(curr_hit)
            hsps += curr_hsps
            hits, hsps = load_into_db(self, hits, hsps, force=False, raw=True)
        hits, hsps = load_into_db(self, hits, hsps, force=True, raw=True)
        assert len(hits) == 0 or isinstance(hits[0], dict), (hits[0], type(hits[0]))
        assert len(hsps) >= len(hits), (len(hits), len(hsps), hits[0], hsps[0])
    else:
        self.logger.info("Finished reading %s data, starting serialisation with %d processors", bname, procs)
        # Now we have to write down everything inside the temporary files.
        # Params_name must contain: shape of the array, dtype of the array, columns
        params = {"columns": columns}
        with open(params_file, "wb") as pfile:
            pfile.write(msgpack.dumps(params))
        assert os.path.exists(params_file)
        # Split the indices
        for idx, split in enumerate(np.array_split(np.array(list(groups.items()),
                                                            dtype=object), procs)):
            with open(index_files[idx], "wb") as index:
                for item in split:
                    vals = (tuple(item[0]), values[item[1], :].tolist())
                    msgpack.dump(vals, index)
            assert os.path.exists(index_files[idx])
            processes[idx].start()

        try:
            res = [proc.join() for proc in processes]
        except KeyboardInterrupt:
            raise KeyboardInterrupt
        except Exception:
            raise
        finally:
            os.remove(params_file)

    return
