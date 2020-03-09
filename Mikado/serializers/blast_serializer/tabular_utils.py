from Bio.SubsMat import MatrixInfo
from . import Query, Target, Hsp, Hit, InvalidHit
from fastnumbers import fast_int, isint
import re
import numpy as np
import pandas as pd
from ...utilities.log_utils import create_null_logger
# import pysam


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


def _parse_btop(btop, qpos, array: np.array, matrix, qmult=1):

    """Parse the BTOP lines of tabular BLASTX/DIAMOND output.
    In BLASTX, the alignment *never* skips bases. Ie the relationship is *always* 3 bases to 1 aminoacid,
    even when there are gaps. It is therefore possible to perfectly infer the position of identical vs positive
    matches.

    The function expects the starting position of the match, the

    """

    for pos in re.findall(r"(\d+|[A-Z|-]{2,2})", btop):
        if isint(pos):
            pos = fast_int(pos)
            array[:, qpos:qpos + pos * qmult] = 1
            qpos += pos * qmult
        elif pos[0] == "-":  # Gap in query
            continue
        elif pos[1] == "-":  # Gap in target
            qpos += qmult
        else:
            if matrix.get(pos, -1) > 0:
                array[2, qpos:qpos + pos * qmult] = 1
            qpos += qmult

    return array


def prepare_tab_hit():
    """"""


def parse_tab_blast(bname: str, queries: pd.DataFrame, min_evalue: float,
                    matrix_name="blosum62", qmult=3, tmult=1, logger=create_null_logger()):

    """This function will use `pandas` to quickly parse, subset and analyse tabular BLAST files.
    Files should have been generated (whether with NCBI BLASTX or DIAMOND) with the following keys:

    qseqid sseqid pident ppos length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq btop


    """

    matrix_name = matrix_name.lower()
    if matrix_name not in matrices:
        raise KeyError("Matrix {} is not valid. Please specify a valid name.".format(matrix_name))

    matrix = matrices[matrix_name]
    data = pd.read_csv(bname, delimiter="\t", names=blast_keys)
    # TODO: substitute the FAI with a SQL read, otherwise this will be slow and very memory-heavy
    groups = data.groupby(["qseqid", "sseqid"], as_index=True)
    # Switch start and env when they are not in the correct order
    _ix = (data.qstart > data.qend)
    data.loc[_ix, ["qstart", "qend"]] = data.loc[_ix, ["qend", "qstart"]].values
    # Get the minimum evalue for each group
    data = data.join(queries)
    data = data.join(groups.agg(min_evalue=pd.NamedAgg("evalue", np.min),
                                max_bitscore=pd.NamedAgg("bitscore", np.max())[["min_evalue",
                                                                                "max_bitscore"]],
                                on=["qseqid", "sseqid"]))
    data = data[data.min_evalue <= min_evalue]

    hits, hsps = [], []

    for qseqid, group in data.groupby("qseqid"):
        # Sort by minimum evalue
        ssorted = group.sort_values(["min_evalue", "sseqid"], ascending=True)[["min_evalue",
                                                                               "max_bitscore",
                                                                               "sseqid"]].drop_duplicates()
        for hitnum, (_, row) in enumerate(ssorted.iterrows()):

            hit_dict_params = dict()
            (hit_dict_params["query_multiplier"],
             hit_dict_params["target_multiplier"]) = (qmult, tmult)
            hit_evalue = row.min_evalue
            hit_bs = row.max_bitscore
            hit_dict_params["hit_number"] = hitnum
            hit_dict_params["evalue"] = hit_evalue
            hit_dict_params["bits"] = hit_bs

            # Prepare for bulk load
            rows = group[group.sseqid == row.sseqid]
            try:
                hit, hit_hsps = prepare_tab_hit(rows,
                                                **hit_dict_params)
            except InvalidHit as exc:
                logger.error(exc)
                continue
            # hit["query_aligned_length"] = min(record.seq_len, hit["query_aligned_length"])

            hits.append(hit)
            hsps.extend(hit_hsps)


            ogroup = group[group.sseqid == row.sseqid]


    




__author__ = 'Luca Venturini'
