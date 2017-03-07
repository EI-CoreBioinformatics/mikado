cpdef tuple compare(prediction, reference, bint lenient=?, bint strict_strandedness=?)
cdef str __assign_multiexonic_ccode(prediction, reference, long nucl_overlap, double stats[9])
cdef str __assign_monoexonic_ccode(prediction, reference, long nucl_overlap, double stats[9])