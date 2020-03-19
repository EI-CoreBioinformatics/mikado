# cython: infer_types=True
# cython: language_level=3
# distutils: language=c++
import re
import numpy as np
cimport numpy as np
from libcpp.vector cimport vector
from libcpp.string cimport string
cimport cython


btop_pattern = re.compile(r"(\d+|\D{2,2})")

@cython.wraparound(False)
@cython.cdivision(True)
@cython.boundscheck(False)
cpdef parse_btop(str btop, Py_ssize_t qpos, Py_ssize_t spos,
                 np.ndarray[dtype=np.int, ndim=2, cast=True] query_array,
                 np.ndarray[dtype=np.int, ndim=2, cast=True] target_array,
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
    cdef long qstart = qpos
    cdef long sstart = spos

    cdef long[:,:] query_view = query_array
    cdef long[:,:] target_view = target_array
    cdef long aln_span = 0
    cdef long step = min(qmult, tmult)
    cdef string match
    cdef Py_ssize_t idx

    for pos in btop_pattern.findall(btop):
        try:
            ipos = int(pos)
        except ValueError:
            if len(pos) != 2:
                raise ValueError(pos)
            ipos = - 1
        if ipos >= 0:
            aln_span += step * ipos
            for idx in range(step * ipos):
                match.push_back(b"|")
            query_view[:, qpos:qpos + ipos * qmult] = 1
            target_view[:, spos:spos + ipos * tmult] = 1
            qpos += ipos * qmult
            spos += ipos * tmult
        else:
            aln_span += step
            if pos[0] == "-":  # Gap in query
                target_view[0, spos:spos + tmult] = 1
                spos += tmult
                match.push_back(b"-")
            elif pos[1] == "-":  # Gap in target
                query_view[0, qpos:qpos + qmult] = 1
                qpos += qmult
                if pos[0] == "*":
                    match.push_back(b"*")
                else:
                    match.push_back(b"_")
            else:
                query_view[0, qpos:qpos + qmult] = 1
                target_view[0, spos:spos + tmult] = 1
                if matrix.get(pos, -1) > 0:
                    query_view[2, qpos:qpos + qmult] = 1
                    target_view[2, spos:spos + tmult] = 1
                    match.push_back(b"+")
                else:
                    if pos[0] == "*" or pos[1] == "*":
                        match.push_back(b"*")
                    else:
                        match.push_back(b"/")

                qpos += qmult
                spos += tmult

    fmatch = <bytes>match
    fmatch = fmatch.decode()

    return query_array, target_array, aln_span, fmatch
