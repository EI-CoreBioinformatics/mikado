from Bio.SubsMat import MatrixInfo
from sys import maxsize
from fastnumbers import fast_int, isint
import re
import numpy as np

matrices = dict()


def _parse_btop(btop, qpos, matrix, qmult=1, maxpos=maxsize):

    """Parse the BTOP lines of tabular BLASTX/DIAMOND output"""

    ids = set()
    poses = set()

    for pos in re.findall(r"(\d+|[A-Z|-]{2,2})", btop):
        if isint(pos):
            pos = fast_int(pos)
            ids.update(set(range(qpos, min(maxpos, qpos + pos * qmult))))
            poses.update(set(range(qpos,
                                   min(maxpos, qpos + pos * qmult)
                                   )))
            qpos += pos * qmult
        elif pos[0] == "-":  # Gap in query
            continue
        elif pos[1] == "-":  # Gap in target
            qpos += qmult
        else:
            if matrix.get(pos, -1) > 0:
                poses.update(set(range(qpos,
                                       min(maxpos, qpos + pos * qmult)
                                       )))
            qpos += qmult

    return btop, ids, poses



for mname in MatrixInfo.available_matrices:
    matrix = dict()
    omatrix = getattr(MatrixInfo, mname)
    for key, val in omatrix.items():
        if key[::-1] in omatrix and omatrix[key[::-1]] != val:
            raise KeyError((key, val, key[::-1], omatrix[key[::-1]]))
        matrix["".join(key)] = val
        matrix["".join(key[::-1])] = val
    matrices[mname] = matrix


__author__ = 'Luca Venturini'
