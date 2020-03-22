from Bio.SubsMat import MatrixInfo
from functools import partial
from .btop_parser import parse_btop
import re
import numpy as np
import pandas as pd
from ...utilities.log_utils import create_null_logger, create_default_logger
import multiprocessing as mp
from .utils import load_into_db
from collections import namedtuple
from ...utilities.dbutils import connect
from sqlalchemy.orm.session import Session
import logging
import logging.handlers


__author__ = 'Luca Venturini'

# Diamond default: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
# BLASTX default: qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore
blast_keys = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore ppos btop".split()

matrices = dict()
for mname in MatrixInfo.available_matrices:
    matrix = dict()
    omatrix = getattr(MatrixInfo, mname)
    for key, val in omatrix.items():
        if key[::-1] in omatrix and omatrix[key[::-1]] != val:
            raise KeyError((key, val, key[::-1], omatrix[key[::-1]]))
        matrix["".join(key)] = val
        matrix["".join(key[::-1])] = val
    matrices[mname] = matrix


def get_queries(engine):
    queries = pd.read_sql_table("query", engine, index_col="query_name")
    queries.columns = ["qid", "qlength"]
    queries["qid"] = queries["qid"].astype(int)
    assert queries.qid.drop_duplicates().shape[0] == queries.shape[0]
    return queries


def get_targets(engine):
    targets = pd.read_sql_table("target", engine, index_col="target_name")
    targets.columns = ["sid", "slength"]
    targets["sid"] = targets["sid"].astype(int)
    assert targets.sid.drop_duplicates().shape[0] == targets.shape[0]

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
    if not np.isclose(pident, hsp[columns["pident"]], atol=.1, rtol=.1):
        raise ValueError((pident, hsp[columns["pident"]]))
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
    return hit_dict, hsps


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
        assert ~(data[col].isna().any()), (col, data[data[col].isna()].shape[0], data.shape[0])
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
    return data


class Preparer(mp.Process):

    def __init__(self, queue: mp.JoinableQueue,
                 returnqueue: mp.JoinableQueue,
                 identifier: int,
                 matrix_name=None, qmult=3, tmult=1, columns=None, **kwargs):

        super().__init__(**kwargs)
        self.matrix_name = matrix_name
        self.qmult, self.tmult, self.columns = qmult, tmult, columns
        self.queue = queue
        self.returnqueue = returnqueue
        self.identifier = identifier

    def run(self):
        prep_hit = partial(prepare_tab_hit,
                           columns=self.columns, qmult=self.qmult, tmult=self.tmult,
                           matrix_name=self.matrix_name)

        while True:
            data = self.queue.get()
            if data == "EXIT":
                self.queue.put("EXIT")
                self.returnqueue.put((None, self.identifier))
                break
            key, rows = data
            hits, hsps = prep_hit(key, rows)
            self.returnqueue.put((hits, hsps))
        return


def parse_tab_blast(self,
                    bname: str,
                    queries: pd.DataFrame,
                    targets: pd.DataFrame,
                    identifier=None,
                    conf=None,
                    logger=create_null_logger(),
                    level="DEBUG",
                    logging_queue=None,
                    lock=None,
                    matrix_name="blosum62", qmult=3, tmult=1):
    """This function will use `pandas` to quickly parse, subset and analyse tabular BLAST files.
    """

    matrix_name = matrix_name.lower()

    if matrix_name not in matrices:
        raise KeyError("Matrix {} is not valid. Please specify a valid name.".format(matrix_name))

    if logging_queue is not None:
        logger = create_default_logger("parse_tab_blast-{}".format(identifier), level=level)
        [logger.removeHandler(handler) for handler in logger.handlers]
        handler = logging.handlers.QueueHandler(logging_queue)
        logger.addHandler(handler)
        logger.info("Started process %s", identifier)

    if self is None:
        _self = namedtuple("self", ["logger", "session", "engine", "maxobjects"])
        engine = connect(conf, strategy="threadlocal")
        session = Session(bind=engine)
        self = _self(logger, session, engine, int(conf["serialise"]["max_objects"] / conf["threads"]))
        if queries is None:
            queries = get_queries(self.engine)
        if targets is None:
            targets = get_targets(self.engine)

    data = pd.read_csv(bname, delimiter="\t", names=blast_keys)
    data = sanitize_blast_data(data, queries, targets, qmult=qmult, tmult=tmult)
    current_qid, current_sid = None, None
    current_rows = []
    columns = dict((col, idx) for idx, col in enumerate(data.columns))
    prep_hit = partial(prepare_tab_hit, columns=columns, qmult=qmult, tmult=tmult,
                       matrix_name=matrix_name)
    hits, hsps = [], []

    for row in data.itertuples(name=None, index=False):
        qid, sid = row[columns["qid"]], row[columns["sid"]]
        if qid != current_qid or sid != current_sid:
            if current_rows:
                curr_hit, curr_hsps = prep_hit((current_qid, current_sid), current_rows)
                hits.append(curr_hit)
                hsps.extend(curr_hsps)
                hits, hsps = load_into_db(self, hits, hsps, force=False)
            current_qid, current_sid, current_rows = qid, sid, [row]
        else:
            current_rows.append(row)

    if current_rows:
        curr_hit, curr_hsps = prep_hit((current_qid, current_sid), current_rows)
        hits.append(curr_hit)
        hsps += curr_hsps
        hits, hsps = load_into_db(self, hits, hsps, force=True, lock=lock)

    assert len(hits) == 0 or isinstance(hits[0], dict), (hits[0], type(hits[0]))
    assert len(hsps) >= len(hits), (len(hits), len(hsps), hits[0], hsps[0])
    # return hits, hsps
