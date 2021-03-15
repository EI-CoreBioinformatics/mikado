#distutils: language = c++

cdef class Interval:
    cdef public int start, end
    cdef public object data
    # cdef public int begin
    cdef public object value, chrom, strand
    cpdef tuple _as_tuple(Interval self)


cdef class IntervalNode:
    cdef float priority
    cdef public Interval interval
    cdef public int start, end
    cdef int minstop, maxstop, minstart
    cdef IntervalNode cleft, cright, croot
    cpdef tuple _as_tuple(IntervalNode self)
    cdef IntervalNode rotate_right(IntervalNode self)
    cdef IntervalNode rotate_left(IntervalNode self)
    cdef inline void set_stops(IntervalNode self)
    cdef void _intersect( IntervalNode self, int start, int end, list results)
    cdef void _seek_left(IntervalNode self, int position, list results, int n, int max_dist)
    cdef void _seek_right(IntervalNode self, int position, list results, int n, int max_dist)
    cdef void _traverse(IntervalNode self, object func)
    cdef IntervalNode _insert(IntervalNode self, Interval interval)
    cpdef intersect(IntervalNode self, int start, int stop)
    # cpdef right(self, Interval f, int n=1, int max_dist=25000)
    cpdef left(IntervalNode self, Interval f, int n=?, int max_dist=?, bint overlap=?)
    cpdef right(IntervalNode self, Interval f, int n=?, int max_dist=?, bint overlap=?)
    # cpdef left(self, Interval f, int n=1, int max_dist=25000)
    # cpdef IntervalNode insert(IntervalNode self, int start, int end, object value=*)
    cpdef bint fuzzy_equal(IntervalNode self, IntervalNode other, int fuzzy)


cdef class IntervalTree:
    cdef IntervalNode root
    cdef int num_intervals
    cpdef int size(IntervalTree self)
    cpdef find(IntervalTree self, int start, int end, bint strict=?,
               bint contained_check=?, int max_distance=?, int n=?, object value=?)
    # cpdef insert(IntervalTree self, int start, int end, object value=?)
    cpdef insert(IntervalTree self, Interval interval)
    cdef inline void set_stops(IntervalTree self)
    # cdef inline void set_stops(IntervalTree self)