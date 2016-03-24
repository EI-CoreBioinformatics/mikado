cpdef int overlap(first, second, int flank=0, bint positive=0):

    """This function quickly computes the overlap between two
    ranges, with an optional flank."""

    cdef int start, end, ostart, oend

    start, end = first[:2]
    ostart, oend = second[:2]

    if start > end:
        start, end = end, start
    if ostart > oend:
        ostart, oend = oend, ostart

    cdef int right
    cdef int left

    cdef int left_one = start - flank
    cdef int left_two = ostart - flank

    if left_one > left_two:
        left = left_one
    else:
        left = left_two

    cdef int right_one = end + flank
    cdef int right_two = oend + flank

    if right_one < right_two:
        right = right_one
    else:
        right = right_two

    cdef int result = right - left
    if positive == 1 and result < 0:
        return 0
    else:
        return result

    # return result


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