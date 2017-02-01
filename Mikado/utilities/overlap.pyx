import cython

@cython.profile(True)
cpdef long overlap(first, second, long flank=0, bint positive=0):

    """This function quickly computes the overlap between two
    ranges, with an optional flank."""

    cdef long start, end, ostart, oend

    if hasattr(first, "start"):
        start, end = first.start, first.end
    else:
        start, end = first[:2]
    if hasattr(second, "start"):
        ostart, oend = second.start, second.end
    else:
        ostart, oend = second[:2]

    return c_overlap(start, end, ostart, oend, flank=flank, positive=positive)


@cython.profile(True)
cdef long c_overlap(long start, long end, long ostart, long oend, long flank, bint positive):
    if start > end:
        start, end = end, start
    if ostart > oend:
        ostart, oend = oend, ostart

    cdef long right
    cdef long left

    cdef long left_one = start - flank
    cdef long left_two = ostart - flank

    if left_one > left_two:
        left = left_one
    else:
        left = left_two

    cdef long right_one = end + flank
    cdef long right_two = oend + flank

    if right_one < right_two:
        right = right_one
    else:
        right = right_two

    cdef long result = right - left
    if positive == 1 and result < 0:
        return 0
    else:
        return result

# def overlap_positive(int start, int end, int ostart, int oend, int flank):
#
#    """This function quickly computes the overlap between two
#    ranges, with an optional flank, and a minimum result of 0."""
#
#    cdef int result = overlap(start, end, ostart, oend, flank)
#    if result < 0:
#        return 0
#    else:
#        return result