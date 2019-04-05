#distutils: language = c++

cdef class IntervalNode:
    cdef float priority
    cdef public object interval
    cdef public int start, end
    cdef int minend, maxend, minstart
    cdef IntervalNode cleft, cright, croot

    cdef IntervalNode rotate_right(IntervalNode self)
    cdef IntervalNode rotate_left(IntervalNode self)
    cdef inline void set_ends(IntervalNode self)
    cdef void _intersect( IntervalNode self, int start, int end, list results)
    cdef void _seek_left(IntervalNode self, int position, list results, int n, int max_dist)
    cdef void _seek_right(IntervalNode self, int position, list results, int n, int max_dist)
    cdef void _traverse(IntervalNode self, object func)
    cpdef right(self, position, int n=?, int max_dist=?)
    cpdef left(self, position, int n=?, int max_dist=?)
    cpdef IntervalNode insert(IntervalNode self, int start, int end, object value=*)
    cpdef bint fuzzy_equal(IntervalNode self, IntervalNode other, int fuzzy)


cdef class Interval:
    cdef public int start, end
    cdef public int begin
    cdef public object value, chrom, strand

    cpdef tuple _as_tuple(self)


cdef class IntervalTree:
    cdef IntervalNode root
    cdef int num_intervals
    cpdef int size(self)
    cdef inline void set_ends(IntervalTree self)