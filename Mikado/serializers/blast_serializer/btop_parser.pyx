import re
import numpy as np
cimport numpy as np


btop_pattern = re.compile(r"(\d+|[A-Z|-]{2,2})")


def parse_btop(str btop, long qpos, long spos,
                 np.ndarray[np.long, ndim=2, cast=True] query_array,
                 np.ndarray[np.long, ndim = 2, cast=True] target_array,
                 dict matrix, long qmult=3, long tmult=1):

    """Parse the BTOP lines of tabular BLASTX/DIAMOND output.
    In BLASTX, the alignment *never* skips bases. Ie the relationship is *always* 3 bases to 1 aminoacid,
    even when there are gaps. It is therefore possible to perfectly infer the position of identical vs positive
    matches.

    The function expects the starting position of the match, the

    """

    # 0: match
    # 1: identical
    # 2: positive

    cdef str pos
    cdef long ipos

    for pos in btop_pattern.findall(btop):
        try:
            ipos = int(pos)
        except ValueError:
            ipos = - 1
        if ipos >= 0:
            query_array[:, qpos:qpos + ipos * qmult] = 1
            target_array[:, spos:spos + ipos * tmult] = 1
            qpos += ipos * qmult
            spos += ipos * tmult
        elif pos[0] == "-":  # Gap in query
            spos += tmult
        elif pos[1] == "-":  # Gap in target
            qpos += qmult
        else:
            query_array[0, qpos:qpos + qmult] = 1
            target_array[0, spos:spos + tmult] = 1
            if matrix.get(pos, -1) > 0:
                query_array[2, qpos:qpos + qmult] = 1
                target_array[2, spos:spos + tmult] = 1
            qpos += qmult
            spos += tmult

    return query_array, target_array
