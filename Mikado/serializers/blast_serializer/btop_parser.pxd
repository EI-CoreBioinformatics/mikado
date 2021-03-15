cimport numpy as np
import numpy as np


cpdef parse_btop(str btop, Py_ssize_t qpos, Py_ssize_t spos,
                 np.ndarray[dtype=long, ndim=2, cast=True] query_array,
                 np.ndarray[dtype=long, ndim=2, cast=True] target_array,
                 dict matrix, long qmult=?, long tmult=?)