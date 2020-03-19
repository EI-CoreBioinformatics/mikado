from Bio.SubsMat import MatrixInfo
# from . import Query, Target, Hsp, Hit, InvalidHit
from .btop_parser import parse_btop
import re
import numpy as np
import pandas as pd
from ...utilities.log_utils import create_null_logger
import multiprocessing as mp
from .utils import load_into_db


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


def prepare_tab_hsp(key,
                    hsp: pd.Series,
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
    query_array = np.zeros([3, int(hsp.qlength)], dtype=np.int)
    target_array = np.zeros([3, int(hsp.slength)], dtype=np.int)
    matrix = matrices.get(matrix_name, matrices["blosum62"])
    if hsp.qstart < 0:
        raise ValueError(hsp.qstart)
    try:
        query_array, target_array, aln_span, match = parse_btop(
        hsp.btop,
        query_array=query_array,
        target_array=target_array,
        qpos=int(hsp.qstart),
        spos=int(hsp.sstart),
        qmult=qmult,
        tmult=tmult,
        matrix=matrix)
    except TypeError:
        raise TypeError((hsp.btop, hsp.qstart, hsp.sstart, qmult, tmult))
    hsp_dict["counter"] = hsp.hsp_num
    hsp_dict["query_hsp_start"] = hsp.qstart
    hsp_dict["query_hsp_end"] = hsp.qend
    hsp_dict["query_frame"] = hsp.query_frame
    hsp_dict["target_hsp_start"] = hsp.sstart
    hsp_dict["target_hsp_end"] = hsp.send
    hsp_dict["target_frame"] = hsp.target_frame
    span = np.where(query_array[0] > 0)[0]
    if span.shape[0] == 0:
        raise ValueError((hsp.btop, type(hsp.btop)))
    if aln_span != hsp.length:
        raise ValueError((aln_span, hsp.length))
    pident = np.where(query_array[1] > 0)[0].shape[0] / (aln_span * qmult) * 100
    if not np.isclose(pident, hsp.pident, atol=.1, rtol=.1):
        raise ValueError((pident, hsp.pident))
    hsp_dict["hsp_identity"] = pident
    ppos = np.where(query_array[2] > 0)[0].shape[0] / (aln_span * qmult) * 100
    hsp_dict["hsp_positives"] = ppos
    hsp_dict["match"] = match
    hsp_dict["hsp_length"] = hsp.length
    hsp_dict["hsp_bits"] = hsp.bitscore
    hsp_dict["hsp_evalue"] = hsp.evalue
    return key, hsp_dict, query_array, target_array


def prepare_tab_hit(key: tuple,
                    hit: pd.DataFrame,
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
    for hsp in hit.itertuples(index=False):
        if hit_row is None:
            hit_row = hsp
        key, hsp_dict, query_array, target_array = prepare_tab_hsp(key, hsp, qmult=qmult, tmult=tmult,
                                                                   matrix_name=matrix_name)
        query_arrays.append(query_array)
        target_arrays.append(target_array)
        hsps.append(hsp_dict)
    
    qlength = int(hit_row.qlength)
    hit_dict.update(kwargs)
    hit_dict["query_id"] = key[0]
    hit_dict["target_id"] = key[1]
    hit_dict["hit_number"] = int(hit_row.hit_num)
    hit_dict["evalue"] = hit_row.min_evalue
    hit_dict["bits"] = hit_row.max_bitscore
    hit_dict["query_multiplier"] = int(qmult)
    hit_dict["target_multiplier"] = int(tmult)

    query_array = sum(query_arrays)
    target_array = sum(target_arrays)
    q_aligned = np.where(query_array[0] > 0)[0]
    hit_dict["query_aligned_length"] = min(qlength, q_aligned.shape[0])
    qstart, qend = q_aligned.min(), q_aligned.max() + 1
    hit_dict["query_start"], hit_dict["query_end"] = int(qstart), int(qend)
    if int(qend) not in hit["qend"].astype(int).values:
        raise ValueError("Invalid end point: {}, {}".format(qend, hit["qend"]))
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
    hit_dict["target_aligned_length"] = min(t_aligned.shape[0], hit_row.slength)
    hit_dict["target_start"] = int(t_aligned.min())
    hit_dict["target_end"] = int(t_aligned.max() + 1)
    if hit_dict["target_end"] not in hit["send"].values:
        raise ValueError("Invalid target end point: {}, {}".format(hit_dict["target_end"], hit["send"]))
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
    data["hsp_num"] = data.sort_values("evalue").groupby(["qseqid", "sseqid"]).cumcount() + 1
    temp = data[["qseqid", "sseqid", "min_evalue"]].drop_duplicates().sort_values(
        ["min_evalue", "sseqid"], ascending=True)
    temp["hit_num"] = temp.groupby(["qseqid"]).cumcount() + 1
    temp.set_index(["qseqid", "sseqid"], inplace=True)
    data = data.join(temp["hit_num"], on=["qseqid", "sseqid"])
    return data


def parse_tab_blast(self,
                    bname: str,
                    queries: pd.DataFrame,
                    targets: pd.DataFrame,
                    hits: list,
                    hsps: list,
                    pool: (mp.Pool, None),
                    matrix_name="blosum62", qmult=3, tmult=1, logger=create_null_logger()):
    """This function will use `pandas` to quickly parse, subset and analyse tabular BLAST files.
    """

    matrix_name = matrix_name.lower()
    if matrix_name not in matrices:
        raise KeyError("Matrix {} is not valid. Please specify a valid name.".format(matrix_name))

    data = pd.read_csv(bname, delimiter="\t", names=blast_keys)
    data = sanitize_blast_data(data, queries, targets, qmult=qmult, tmult=tmult)
    groups = data.groupby(["qid", "sid"])
    if pool is not None:
        assert isinstance(pool, mp.pool.Pool)
        for res in (pool.apply_async(
                prepare_tab_hit, args=(key, group), kwds={"qmult": qmult, "tmult": tmult,
                                                          "matrix_name": matrix_name}) for key, group in groups):
            curr_hit, curr_hsps = res.get()
            hits.append(curr_hit)
            hsps.extend(curr_hsps)
            tot = len(hits) + len(hsps)
            if tot >= self.maxobjects:
                hits, hsps = load_into_db(self, hits, hsps, force=False)
    else:
        for key, group in groups:
            curr_hit, curr_hsps = prepare_tab_hit(key, group, qmult=qmult, tmult=tmult, matrix_name=matrix_name)
            hits.append(curr_hit)
            hsps += curr_hsps
            hits, hsps = load_into_db(self, hits, hsps, force=False)

    assert len(hits) == 0 or isinstance(hits[0], dict), (hits[0], type(hits[0]))
    assert len(hsps) >= len(hits), (len(hits), len(hsps), hits[0], hsps[0])
    return hits, hsps
