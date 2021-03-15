# cython: language_level=3

import numpy as np
cimport numpy as np
from libc.stdlib cimport malloc
from libcpp.string cimport string as cstring
cimport cython

DTYPE = np.int_
ctypedef np.int_t DTYPE_t


@cython.cdivision(True)
cdef _analyze_string(char* qseq, char* sseq, char* mid,
                     long query_start, long query_end, long query_length, long qmult):

    cdef long qpos = -1
    cdef long shape = <long> (query_end - query_start) / qmult
    cdef np.ndarray[DTYPE_t, ndim=1] query_array = np.zeros(shape, dtype=DTYPE)
    cdef DTYPE_t[:] query_view = query_array
    cdef Py_ssize_t match_len = len(mid)
    cdef cstring match
    cdef char qchar, schar, midchar
    cdef Py_ssize_t idx
    if shape > query_length:
        raise ValueError((shape, query_length, query_start, query_end, qmult, mid))

    for idx in range(match_len):
        qchar, schar, midchar = qseq[idx], sseq[idx], mid[idx]
        if qchar == b"-":
            if schar == b"*":
                match.push_back(b"*")
            else:
                match.push_back(b"-")
            continue
        else:
            qpos += 1
            if qchar == b"*" or schar == b"*":
                match.push_back(b"*")
            elif qchar == schar:
                match.push_back(b"|")
                query_view[qpos] = 2
            elif midchar == b"+":
                match.push_back(b"+")
                query_view[qpos] = 1
            elif schar == b"-":
                if qchar == b"*":
                    match.push_back(b"*")
                else:
                    match.push_back(b"_")
                continue
            else:
                match.push_back(b"/")

    return query_array, match


cpdef prepare_aln_strings(hsp, bint off_by_one=0, long qmultiplier=1):

    """This private method calculates the identical positions, the positives, and a re-factored match line
    starting from the HSP.
    :type hsp: Bio.SearchIO.HSP
    """

    cdef bytes qseq = <bytes> str(hsp.query.seq).encode()
    cdef bytes sseq = <bytes> str(hsp.hit.seq).encode()
    cdef bytes mid = <bytes> str(hsp.aln_annotation["similarity"]).encode()
    cdef long match_len = len(mid)

    lett_array = np.array([list(str(hsp.query.seq)),
                           list(str(hsp.aln_annotation["similarity"])),
                           list(str(hsp.hit.seq))], dtype=np.str_)

    if lett_array.shape[0] == 0 or len(lett_array.shape) == 1 or lett_array.shape[1] == 0:  # Empty array!
        raise ValueError("Empty array of matches! {}".format("\n".join(
            [hsp.query.seq, hsp.hit.seq, hsp.aln_annotation["similarity"]])))

    cdef np.ndarray[DTYPE_t, ndim=1] query_array
    cdef long query_start = hsp.query_start
    cdef long query_end = hsp.query_end - off_by_one
    cdef long query_length = len(qseq)
    cdef np.ndarray[DTYPE_t, ndim=2] summer = np.array([[_] for _ in range(qmultiplier)])
    query_array, match = _analyze_string(qseq, sseq, mid, query_start, query_end, query_length, qmultiplier)
    cdef np.ndarray[DTYPE_t, ndim=1] _id_catcher = np.where(query_array >= 2)[0]
    identical_positions = ((_id_catcher * qmultiplier) + summer).flatten()
    cdef np.ndarray[DTYPE_t, ndim=1] _pos_catcher = np.where(query_array >= 1)[0]
    positives = ((_pos_catcher * qmultiplier) + summer).flatten()
    if hsp.query_frame > 0:
        identical_positions = identical_positions + hsp.query_start
        positives = positives + hsp.query_start
    else:
        identical_positions = hsp.query_end - identical_positions - 1
        positives = hsp.query_end - positives - 1

    # identical_positions = set(identical_positions)
    # positives = set(positives)
    identical_positions.sort()
    positives.sort()
    match = <bytes>match
    match = match.decode()

    return match, identical_positions, positives
