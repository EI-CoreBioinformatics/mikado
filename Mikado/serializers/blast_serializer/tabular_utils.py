from Bio.SubsMat import MatrixInfo
from . import Query, Target, Hsp, Hit, InvalidHit
from fastnumbers import fast_int, isint
import re
import numpy as np
import pandas as pd
from ...utilities.log_utils import create_null_logger
import multiprocessing as mp
from collections import defaultdict


__author__ = 'Luca Venturini'


blast_keys = "qseqid sseqid pident ppos length mismatch gapopen qstart "\
            "qend sstart send evalue bitscore qseq sseq btop".split()

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


btop_pattern = re.compile(r"(\d+|[A-Z|-]{2,2})")


def _parse_btop(btop, qpos, spos,
                query_array: np.array,
                target_array: np.array,
                matrix, qmult=3, tmult=1):

    """Parse the BTOP lines of tabular BLASTX/DIAMOND output.
    In BLASTX, the alignment *never* skips bases. Ie the relationship is *always* 3 bases to 1 aminoacid,
    even when there are gaps. It is therefore possible to perfectly infer the position of identical vs positive
    matches.

    The function expects the starting position of the match, the

    """

    # 0: match
    # 1: identical
    # 2: positive

    qpos, spos, qmult, tmult = int(qpos), int(spos), int(qmult), int(tmult)

    for pos in re.findall(btop_pattern, btop):
        if isint(pos):
            pos = fast_int(pos)
            try:
                query_array[:, qpos:qpos + pos * qmult] = 1
            except TypeError:
                assert isinstance(qpos, (np.int, int)), (qpos, type(qpos))
                assert isinstance(pos, (np.int, int))
                assert isinstance(qmult, (np.int, int))
            target_array[:, spos:spos + pos * tmult] = 1
            qpos += pos * qmult
            spos += pos * tmult
        elif pos[0] == "-":  # Gap in query
            spos += tmult
        elif pos[1] == "-":  # Gap in target
            qpos += qmult
        else:
            try:
                query_array[0, qpos:qpos + qmult] = 1
                target_array[0, spos:spos + tmult] = 1
            except TypeError:
                raise TypeError((qpos, qpos + qmult, spos, spos + tmult))
            if matrix.get(pos, -1) > 0:
                query_array[2, qpos:qpos + qmult] = 1
                target_array[2, spos:spos + tmult] = 1
            qpos += qmult
            spos += tmult

    return query_array, target_array


def prepare_tab_hsp(hsp, query_array, target_array, qmult=3, tmult=1, matrix_name=None):

    r"""
    Prepare a HSP for loading into the DB.
    The match line will be reworked in the following way:

    - If the position is a match/positive, keep the original value
    - If the position is a gap *for the query*, insert a - (dash)
    - If the position is a gap *for the target*, insert a _ (underscore)
    - If the position is a gap *for both*, insert a \ (backslash)

    :param hsp: A tabular blast row
    :type hsp: pd.Series
    :param counter: a digit that indicates the priority of the HSP in the hit
    :return: hsp_dict, numpy array
    :rtype: (dict, np.array, np.array)
    """

    hsp_dict = dict()
    # We must start from 1, otherwise MySQL crashes as its indices start from 1 not 0
    matrix = matrices.get(matrix_name, matrices["blosum62"])
    query_array, target_array = _parse_btop(hsp.btop,
                                            query_array=query_array,
                                            target_array=target_array,
                                            qpos=hsp.qstart,
                                            spos=hsp.sstart,
                                            qmult=qmult,
                                            tmult=tmult,
                                            matrix=matrix)
    hsp_dict["counter"] = hsp.hsp_num
    hsp_dict["query_hsp_start"] = hsp.qstart
    hsp_dict["query_hsp_end"] = hsp.qend
    hsp_dict["query_frame"] = hsp.query_frame
    hsp_dict["target_hsp_start"] = hsp.sstart
    hsp_dict["target_hsp_end"] = hsp.send
    hsp_dict["target_frame"] = hsp.target_frame
    hsp_dict["hsp_identity"] = hsp.pident
    hsp_dict["hsp_positives"] = hsp.ppos
    hsp_dict["match"] = hsp.btop
    hsp_dict["hsp_length"] = hsp.aln_span
    hsp_dict["hsp_bits"] = hsp.bitscore
    hsp_dict["hsp_evalue"] = hsp.evalue
    return hsp_dict, query_array, target_array


def prepare_tab_hit(hit: pd.DataFrame, qmult=3, tmult=1, matrix_name=None, **kwargs):
    """"""

    hit_dict = dict()
    hsp_dict_list = []

    qlength = hit.qlength.unique()[0]
    assert isinstance(qmult, (int, float)), (qmult, type(qmult))
    assert isinstance(tmult, (int, float)), (tmult, type(tmult))
    hit_dict.update(kwargs)
    hit_dict["query_id"] = int(hit.qid.unique()[0])
    hit_dict["target_id"] = int(hit.sid.unique()[0])
    hit_dict["hit_number"] = int(hit.hit_num.unique()[0])
    query_array = np.zeros([3, hit["qlength"].unique()[0]])
    target_array = np.zeros([3, hit["slength"].unique()[0]])
    hit_dict["evalue"] = hit.min_evalue.unique()[0]
    hit_dict["bits"] = hit.max_bitscore.unique()[0]
    hit_dict["query_multiplier"] = int(qmult)
    hit_dict["target_multiplier"] = int(tmult)

    for idx, hsp in hit.sort_values(
            by=["evalue", "bitscore"], ascending=[True, False]).iterrows():
        hsp_dict, query_array, target_array = prepare_tab_hsp(
            hsp, query_array, target_array, qmult=qmult, matrix_name=matrix_name)
        hsp_dict["query_id"] = hit_dict["query_id"]
        hsp_dict["target_id"] = hit_dict["target_id"]
        hsp_dict_list.append(hsp_dict)

    # q_merged_intervals, q_aligned = merge(q_intervals)
    q_aligned = np.where(query_array[0] > 0)[0]
    hit_dict["query_aligned_length"] = min(qlength, q_aligned.shape[0])
    qstart, qend = q_aligned.min(), q_aligned.max()
    hit_dict["query_start"], hit_dict["query_end"] = int(qstart), int(qend)
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
    hit_dict["target_aligned_length"] = min(t_aligned.shape[0], hit.slength.unique()[0])
    hit_dict["target_start"] = int(t_aligned.min())
    hit_dict["target_end"] = int(t_aligned.max())
    hit_dict["global_identity"] = identical_positions * 100 / q_aligned.shape[0]
    hit_dict["global_positives"] = positives * 100 / q_aligned.shape[0]

    return hit_dict, hsp_dict_list


id_pattern = re.compile(r"^[^\|]*\|([^\|]*)\|.*")


def parse_tab_blast(bname: str,
                    queries: pd.DataFrame,
                    targets: pd.DataFrame,
                    hits: list,
                    hsps: list,
                    pool: (mp.Pool, None),
                    matrix_name="blosum62", qmult=3, tmult=1, logger=create_null_logger()):
    """This function will use `pandas` to quickly parse, subset and analyse tabular BLAST files.
    Files should have been generated (whether with NCBI BLASTX or DIAMOND) with the following keys:

    qseqid sseqid pident ppos length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq btop


    """

    matrix_name = matrix_name.lower()
    if matrix_name not in matrices:
        raise KeyError("Matrix {} is not valid. Please specify a valid name.".format(matrix_name))

    # matrix = matrices[matrix_name]
    data = pd.read_csv(bname, delimiter="\t", names=blast_keys)
    if data[data.btop.isna()].shape[0] > 0:
        raise ValueError(data.loc[0])

    #Correct the names of queries and targets
    data["qseqid"] = data["qseqid"].str.replace(id_pattern, "\\1")
    data["sseqid"] = data["sseqid"].str.replace(id_pattern, "\\1")
    for col in ["qstart", "qend", "sstart", "send"]:
        data[col] = data[col].astype(int)

    # Switch start and env when they are not in the correct order
    _ix = (data.qstart > data.qend)
    data.loc[_ix, ["qstart", "qend"]] = data.loc[_ix, ["qend", "qstart"]].values
    # Get the minimum evalue for each group
    data = data.join(queries, on=["qseqid"]).join(targets, on=["sseqid"]).join(
        data.groupby(["qseqid", "sseqid"]).agg(
            min_evalue=pd.NamedAgg("evalue", np.min),
            max_bitscore=pd.NamedAgg("bitscore", np.max)
        )[["min_evalue", "max_bitscore"]], on=["qseqid", "sseqid"])
    data["query_frame"] = data.qstart % qmult
    data["target_frame"] = data.qstart % tmult
    data["aln_span"] = data.qend - data.qstart + 1

    # Set the hsp_num
    data["hsp_num"] = data.sort_values("evalue").groupby(["qseqid", "sseqid"]).cumcount() + 1
    temp = data[["qseqid", "sseqid", "min_evalue"]].drop_duplicates().sort_values("min_evalue", ascending=True)
    temp["hit_num"] = temp.groupby(["qseqid"]).cumcount() + 1
    temp.set_index(["qseqid", "sseqid"], inplace=True)
    data = data.join(temp["hit_num"], on=["qseqid", "sseqid"])

    if pool is not None:
        assert isinstance(pool, mp.pool.Pool)
        results = []

    previous = len(hits)

    for idx, (qseqid, group) in enumerate(data.groupby(["qseqid", "sseqid"]), 1):
        # Sort by minimum evalue
        # prepare_tab_hit(hit: pd.DataFrame, qmult=3, tmult=1, matrix_name=None, **kwargs):
        if True:
            ghit, ghsps = prepare_tab_hit(group, matrix_name=matrix_name, qmult=qmult, tmult=tmult)
            hits.append(ghit)
            hsps.extend(ghsps)
        else:
            results.append(pool.apply_async(prepare_tab_hit, args=(group, qmult, tmult, matrix_name)))

    if pool is not None:
        for res in results:
            ghit, ghsps = res.get()
            hits.append(ghit)
            hsps.extend(ghsps)

    assert len(hits) - previous == data.groupby(["qseqid", "sseqid"]).ngroups, (
        idx, len(hits) - previous, data.groupby(["qseqid", "sseqid"]).ngroups)

    return hits, hsps
