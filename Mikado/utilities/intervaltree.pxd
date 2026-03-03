# distutils: language = c++

from libc.stdint cimport int64_t


cdef class Interval:
    cdef public int64_t start, end
    cdef public object data
    cdef public object value, chrom, strand
    cpdef tuple _as_tuple(Interval self)


cdef class IntervalNode:
    cdef float priority
    cdef public Interval interval
    cdef public int64_t start, end
    cdef int64_t minstop, maxstop, minstart
    cdef IntervalNode cleft, cright, croot

    cpdef tuple _as_tuple(IntervalNode self)
    cdef IntervalNode rotate_right(IntervalNode self)
    cdef IntervalNode rotate_left(IntervalNode self)
    cdef inline void set_stops(IntervalNode self)

    cdef void _intersect(IntervalNode self, int64_t start, int64_t end, list results)
    cdef void _seek_left(IntervalNode self, int64_t position, list results,
                         int64_t n, int64_t max_dist)
    cdef void _seek_right(IntervalNode self, int64_t position, list results,
                          int64_t n, int64_t max_dist)

    cdef void _traverse(IntervalNode self, object func)
    cdef IntervalNode _insert(IntervalNode self, Interval interval)

    cpdef intersect(IntervalNode self, int64_t start, int64_t stop)

    cpdef left(IntervalNode self, Interval f,
               int64_t n=?, int64_t max_dist=?, bint overlap=?)

    cpdef right(IntervalNode self, Interval f,
                int64_t n=?, int64_t max_dist=?, bint overlap=?)

    cpdef bint fuzzy_equal(IntervalNode self, IntervalNode other, int64_t fuzzy)


cdef class IntervalTree:
    cdef IntervalNode root
    cdef int64_t num_intervals

    cpdef int64_t size(IntervalTree self)

    cpdef find(IntervalTree self,
               int64_t start, int64_t end,
               bint strict=?, bint contained_check=?,
               int64_t max_distance=?, int64_t n=?,
               object value=?)

    cpdef insert(IntervalTree self, Interval interval)

    cdef inline void set_stops(IntervalTree self)
