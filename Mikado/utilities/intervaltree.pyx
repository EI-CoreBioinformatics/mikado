"""
Data structure for performing intersect queries on a set of intervals which
preserves all information about the intervals (unlike bitset projection methods).

:Authors: James Taylor (james@jamestaylor.org),
          Ian Schenk (ian.schenck@gmail.com),
          Brent Pedersen (bpederse@gmail.com)
"""

# Historical note:
#    This module original contained an implementation based on sorted endpoints
#    and a binary search, using an idea from Scott Schwartz and Piotr Berman.
#    Later an interval tree implementation was implemented by Ian for Galaxy's
#    join tool (see `bx.intervals.operations.quicksect.py`). This was then
#    converted to Cython by Brent, who also added support for
#    upstream/downstream/neighbor queries. This was modified by James to
#    handle half-open intervals strictly, to maintain sort order, and to
#    implement the same interface as the original Intersecter.

#distutils: language = c++
#cython: cdivision=True

import operator
from libcpp.vector cimport vector
from cpython cimport array

cdef extern from "stdlib.h":
    int ceil(float f)
    float log(float f)
    int RAND_MAX
    int rand()
    int strlen(char *)
    int iabs(int)

cdef inline int imax2(int a, int b):
    if b > a: return b
    return a

cdef inline int imax3(int a, int b, int c):
    if b > a:
        if c > b:
            return c
        return b
    if a > c:
        return a
    return c

cdef inline int imin3(int a, int b, int c):
    if b < a:
        if c < b:
            return c
        return b
    if a < c:
        return a
    return c

cdef inline int imin2(int a, int b):
    if b < a: return b
    return a

cdef float nlog = -1.0 / log(0.5)


cdef class IntervalNode:
    """
    A single node of an `IntervalTree`.

    NOTE: Unless you really know what you are doing, you probably should use
          `IntervalTree` rather than using this directly.
    """

    property left_node:
        def __get__(self):
            return self.cleft if self.cleft is not EmptyNode else None
    property right_node:
        def __get__(self):
            return self.cright if self.cright is not EmptyNode else None
    property root_node:
        def __get__(self):
            return self.croot if self.croot is not EmptyNode else None

    property value:
        def __get__(self):
            return self.interval.value

    def __repr__(self):
        return "IntervalNode(%i, %i)" % (self.start, self.end)

    def __cinit__(IntervalNode self, int start, int end, object interval):
        # Python lacks the binomial distribution, so we convert a
        # uniform into a binomial because it naturally scales with
        # tree size.  Also, python's uniform is perfect since the
        # upper limit is not inclusive, which gives us undefined here.
        self.priority = ceil(nlog * log(-1.0/(1.0 * rand()/RAND_MAX - 1)))
        self.start    = start
        self.end      = end
        if not isinstance(interval, (Interval, IntervalNode)):
            interval = Interval(start, end, interval)
        self.interval = interval
        self.maxend   = end
        self.minstart = start
        self.minend   = end
        self.cleft    = EmptyNode
        self.cright   = EmptyNode
        self.croot    = EmptyNode

    def __getstate__(self):
       state = []
       state.append(self.priority)
       state.append(self.start)

    def __setstate__(self, state):
       pass

    cpdef IntervalNode insert(IntervalNode self, int start, int end, object value=None):
        """
        Insert a new IntervalNode into the tree of which this node is
        currently the root. The return value is the new root of the tree (which
        may or may not be this node!)
        """
        cdef IntervalNode croot = self
        # If starts are the same, decide which to add interval to based on
        # end, thus maintaining sortedness relative to start/end
        cdef int decision_endpoint = start
        if isinstance(value, Interval):
            value = value.value

        if start == self.start:
            decision_endpoint = end

        if decision_endpoint > self.start:
            # insert to cright tree
            if self.cright is not EmptyNode:
                self.cright = self.cright.insert( start, end, value )
            else:
                self.cright = IntervalNode( start, end, value)
            # rebalance tree
            if self.priority < self.cright.priority:
                croot = self.rotate_left()
        else:
            # insert to cleft tree
            if self.cleft is not EmptyNode:
                self.cleft = self.cleft.insert( start, end, value)
            else:
                self.cleft = IntervalNode( start, end, value)
            # rebalance tree
            if self.priority < self.cleft.priority:
                croot = self.rotate_right()

        croot.set_ends()
        self.cleft.croot  = croot
        self.cright.croot = croot
        return croot

    cdef IntervalNode rotate_right(IntervalNode self):
        cdef IntervalNode croot = self.cleft
        self.cleft  = self.cleft.cright
        croot.cright = self
        self.set_ends()
        return croot

    cdef IntervalNode rotate_left(IntervalNode self):
        cdef IntervalNode croot = self.cright
        self.cright = self.cright.cleft
        croot.cleft  = self
        self.set_ends()
        return croot

    cdef inline void set_ends(IntervalNode self):
        if self.cright is not EmptyNode and self.cleft is not EmptyNode:
            self.maxend = imax3(self.end, self.cright.maxend, self.cleft.maxend)
            self.minend = imin3(self.end, self.cright.minend, self.cleft.minend)
            self.minstart = imin3(self.start, self.cright.minstart, self.cleft.minstart)
        elif self.cright is not EmptyNode:
            self.maxend = imax2(self.end, self.cright.maxend)
            self.minend = imin2(self.end, self.cright.minend)
            self.minstart = imin2(self.start, self.cright.minstart)
        elif self.cleft is not EmptyNode:
            self.maxend = imax2(self.end, self.cleft.maxend)
            self.minend = imin2(self.end, self.cleft.minend)
            self.minstart = imin2(self.start, self.cleft.minstart)


    def intersect( self, int start, int end, sort=True ):
        """
        given a start and a end, return a list of features
        falling within that range
        """
        cdef list results = []
        self._intersect( start, end, results )
        return results

    find = intersect

    cdef void _intersect( IntervalNode self, int start, int end, list results):
        # Left subtree
        if self.cleft is not EmptyNode and self.cleft.maxend > start:
            self.cleft._intersect( start, end, results )
        # This interval
        if ( self.end >= start ) and ( self.start <= end ):
            results.append( self.interval )
        # Right subtree
        if self.cright is not EmptyNode and self.start < end:
            self.cright._intersect( start, end, results )

    def __getitem__(self, int index):
        self.set_ends()
        if index == 0:
            return self.start
        elif index == 1:
            return self.end
        elif index == 2:
            return self.value
        else:
            return [self.start, self.end, self.value][index]
            # raise IndexError("Intervals only have starts and ends!")

    cdef void _seek_left(IntervalNode self, int position, list results, int n, int max_dist):
        # we know we can bail in these 2 cases.
        if self.maxend + max_dist < position:
            return
        if self.minstart > position:
            return

        # the ordering of these 3 blocks makes it so the results are
        # ordered nearest to farest from the query position
        if self.cright is not EmptyNode:
            self.cright._seek_left(position, results, n, max_dist)

        if -1 < position - self.end < max_dist:
            results.append(self.interval)

        # TODO: can these conditionals be more stringent?
        if self.cleft is not EmptyNode:
                self.cleft._seek_left(position, results, n, max_dist)

    cdef void _seek_right(IntervalNode self, int position, list results, int n, int max_dist):
        # we know we can bail in these 2 cases.
        if self.maxend < position: return
        if self.minstart - max_dist > position: return

        #print "SEEK_RIGHT:",self, self.cleft, self.maxend, self.minstart, position

        # the ordering of these 3 blocks makes it so the results are
        # ordered nearest to farest from the query position
        if self.cleft is not EmptyNode:
                self.cleft._seek_right(position, results, n, max_dist)

        if -1 < self.start - position < max_dist:
            results.append(self.interval)

        if self.cright is not EmptyNode:
                self.cright._seek_right(position, results, n, max_dist)

    cpdef left(self, position, int n=1, int max_dist=2500):
        """
        find n features with a start > than `position`
        f: a Interval object (or anything with an `end` attribute)
        n: the number of features to return
        max_dist: the maximum distance to look before giving up.
        """
        cdef list results = []
        # use start - 1 because .left() assumes strictly left-of
        self._seek_left( position - 1, results, n, max_dist )
        if len(results) == n: return results
        r = results
        r.sort(key=operator.attrgetter('end'), reverse=True)
        return r[:n]

    cpdef right(self, position, int n=1, int max_dist=2500):
        """
        find n features with a end < than position
        f: a Interval object (or anything with a `start` attribute)
        n: the number of features to return
        max_dist: the maximum distance to look before giving up.
        """
        cdef list results = []
        # use end + 1 becuase .right() assumes strictly right-of
        self._seek_right(position + 1, results, n, max_dist)

        if len(results) == n: return results
        r = results
        r.sort(key=operator.attrgetter('start'))
        return r[:n]

    def traverse(self, func):
        self._traverse(func)

    cdef void _traverse(IntervalNode self, object func):
        if self.cleft is not EmptyNode: self.cleft._traverse(func)
        func(self)
        if self.cright is not EmptyNode: self.cright._traverse(func)

    def __iter__(self):
        if self.cleft is not EmptyNode:
            for iv in self.cleft:
                yield iv
        yield self
        if self.cright is not EmptyNode:
            for iv in self.cright:
                yield iv

    def __hash__(self):
        self.set_ends()
        return hash(tuple([self.start, self.end]))

    def __eq__(IntervalNode self, IntervalNode other):

        if not isinstance(other, IntervalNode):
            return False

        return self.start == other.start and self.end == other.end

    cpdef bint fuzzy_equal(IntervalNode self, IntervalNode other, int fuzzy):

        cdef int istart, iend
        cdef int ostart, oend
        if not isinstance(other, IntervalNode):
            return False

        istart, iend = self.start, self.end
        ostart, oend = other.start, other.end
        if abs(ostart - istart) <= fuzzy and abs(oend - iend) <= fuzzy:
            return True
        return False


cdef IntervalNode EmptyNode = IntervalNode( 0, 0, Interval(0, 0))

## ---- Wrappers that retain the old interface -------------------------------

cdef class Interval:
    """
    Basic feature, with required integer start and end properties.
    Also accepts optional strand as +1 or -1 (used for up/downstream queries),
    a name, and any arbitrary data is sent in on the info keyword argument

    >>> from bx.intervals.intersection import Interval

    >>> f1 = Interval(23, 36)
    >>> f2 = Interval(34, 48, value={'chr':12, 'anno':'transposon'})
    >>> f2
    Interval(34, 48, value={'anno': 'transposon', 'chr': 12})

    """

    def __init__(self, int start, int end, object value=None, object chrom=None, object strand=None ):
        assert start <= end, "start must be less than end"
        self.start  = start
        self.begin = self.start  # Necessary for compatibility reasons
        self.end   = end
        self.value = value
        self.chrom = chrom
        self.strand = strand

    def __repr__(self):
        fstr = "Interval(%d, %d" % (self.start, self.end)
        if not self.value is None:
            fstr += ", value=" + str(self.value)
        fstr += ")"
        return fstr

    def __richcmp__(self, other, op):
        if op == 0:
            # <
            return self[0] < other[0] or self[1] < other[1]
        elif op == 1:
            # <=
            return self == other or self < other
        elif op == 2:
            # ==
            return self[0] == other[0] and self[1] == other[1]
        elif op == 3:
            # !=
            return self[0] != other[0] or self[1] != other[1]
        elif op == 4:
            # >
            return self[0] > other[0] or self[1] > other[1]
        elif op == 5:
            # >=
            return self == other or self > other

    def __getitem__(self, int index):
        if index == 0:
            return self.start
        elif index == 1:
            return self.end
        elif index == 2:
            return self.value
        else:
            return [self.start, self.end, self.value][index]
            # raise IndexError("Intervals only have starts and ends!")

    def __iter__(self):

        return iter([self.start, self.end])

    def __getstate__(self):

        return [self.start, self.end, self.value, self.chrom, self.strand]

    def __setstate__(self, state):

        self.start, self.end, self.value, self.chrom, self.strand = state

    cpdef tuple _as_tuple(self):
        return (self.start, self.end)

    def __hash__(self):
        return hash(self._as_tuple())


cdef class IntervalTree:
    """
    Data structure for performing window intersect queries on a set of
    of possibly overlapping 1d intervals.

    Create an empty IntervalTree

    >>> from bx.intervals.intersection import Interval, IntervalTree
    >>> intersecter = IntervalTree()

    An interval is a start and end position and a value (possibly None).
    You can add any object as an interval:

    >>> intersecter.insert( 0, 10, "food" )
    >>> intersecter.insert( 3, 7, dict(foo='bar') )

    >>> intersecter.find( 2, 5 )
    ['food', {'foo': 'bar'}]

    If the object has start and end attributes (like the Interval class) there
    is are some shortcuts:

    >>> intersecter = IntervalTree()
    >>> intersecter.insert_interval( Interval( 0, 10 ) )
    >>> intersecter.insert_interval( Interval( 3, 7 ) )
    >>> intersecter.insert_interval( Interval( 3, 40 ) )
    >>> intersecter.insert_interval( Interval( 13, 50 ) )

    >>> intersecter.find( 30, 50 )
    [Interval(3, 40), Interval(13, 50)]
    >>> intersecter.find( 100, 200 )
    []

    Before/after for intervals

    >>> intersecter.before_interval( Interval( 10, 20 ) )
    [Interval(3, 7)]
    >>> intersecter.before_interval( Interval( 5, 20 ) )
    []

    Upstream/downstream

    >>> intersecter.upstream_of_interval(Interval(11, 12))
    [Interval(0, 10)]
    >>> intersecter.upstream_of_interval(Interval(11, 12, strand="-"))
    [Interval(13, 50)]

    >>> intersecter.upstream_of_interval(Interval(1, 2, strand="-"), num_intervals=3)
    [Interval(3, 7), Interval(3, 40), Interval(13, 50)]


    """

    property start:
        def __get__(self):
            self.set_ends()
            return self.root.minstart

    property end:
        def __get__(self):
            self.set_ends()
            return self.root.maxend

    def __cinit__( self ):
        root = None
        num_intervals = 0

    # ---- Pickling ------

    def __getstate__(self):

        return [self.root, self.num_intervals]

    def __setstate__(self, state):

        self.root, self.num_intervals = state


    # ---- Position based interfaces -----------------------------------------

    def insert( self, int start, int end, object value=None ):
        """
        Insert the interval [start,end) associated with value `value`.
        """
        if self.root is None:
            self.root = IntervalNode( start, end, value )
        else:
            self.root = self.root.insert(start, end, value)
        self.num_intervals += 1

    add = insert

    def find(self, int start, int end, bint strict=0, bint contained_check=0, int max_distance=0,
             int num_intervals=1000, object value=None):
        """
        Return a sorted list of all intervals overlapping [start,end).
        If strict is set to True, only matches which are completely contained will
        be counted as full matches. If set to False, also partial matches will count.
        """

        if self.root is None:
            return []

        if max_distance == 0:
            found=self.root.find( start, end )
        else:
            found = self.root.find( start, end )
            found.extend(self.before( start, num_intervals=num_intervals, max_dist=max_distance))
            found.extend(self.after( end, num_intervals=num_intervals, max_dist=max_distance))
            # Emulate the behaviour of the intervaltree library

        if contained_check is True:
            new_found = []
            for _ in found:
                if (_.start >= start and _.end <= end) or (start >= _.start and end <= _.end ):
                    if value is None or (value is not None and _.value == value):
                        new_found.append(_)
            found = new_found

        elif strict is True:
            # cdef list new_found
            new_found = []
            for _ in found:
                if _.start >= start and _.end <= end:
                    if value is None or (value is not None and _.value == value):
                        new_found.append(_)
            found = new_found

        elif value is not None:
            found = [_ for _ in found if _.value == value]

        return found

    search = find

    def __len__(self):
        """Return the number of intervals in the tree."""

        return self.num_intervals

    cpdef int size(self):
        return self.num_intervals

    def before( self, position, num_intervals=1, max_dist=2000):
        """
        Find `num_intervals` intervals that lie before `position` and are no
        further than `max_dist` positions away
        """
        if self.root is None:
            return []
        return self.root.left( position, num_intervals, max_dist)

    def after( self, position, num_intervals=1, max_dist=25000):
        """
        Find `num_intervals` intervals that lie after `position` and are no
        further than `max_dist` positions away
        """
        if self.root is None:
            return []
        return self.root.right( position, num_intervals, max_dist)

    # ---- Interval-like object based interfaces -----------------------------

    def insert_interval( self, interval ):
        """
        Insert an "interval" like object (one with at least start and end
        attributes)
        """
        self.insert(interval.start, interval.end, interval.value)
        # self.num_intervals += 1

    add_interval = insert_interval

    def before_interval( self, interval, num_intervals=1, max_dist=2500 ):
        """
        Find `num_intervals` intervals that lie completely before `interval`
        and are no further than `max_dist` positions away
        """
        if self.root is None:
            return []
        return self.root.left( interval.start, num_intervals, max_dist )

    def after_interval( self, interval, num_intervals=1, max_dist=2500 ):
        """
        Find `num_intervals` intervals that lie completely after `interval` and
        are no further than `max_dist` positions away
        """
        if self.root is None:
            return []
        return self.root.right( interval.end, num_intervals, max_dist )

    def upstream_of_interval( self, interval, num_intervals=1, max_dist=2500 ):
        """
        Find `num_intervals` intervals that lie completely upstream of
        `interval` and are no further than `max_dist` positions away
        """
        if self.root is None:
            return []
        if interval.strand == -1 or interval.strand == "-":
            return self.root.right( interval.end, num_intervals, max_dist )
        else:
            return self.root.left( interval.start, num_intervals, max_dist )

    def downstream_of_interval( self, interval, num_intervals=1, max_dist=2500 ):
        """
        Find `num_intervals` intervals that lie completely downstream of
        `interval` and are no further than `max_dist` positions away
        """
        if self.root is None:
            return []
        if interval.strand == -1 or interval.strand == "-":
            return self.root.left( interval.start, num_intervals, max_dist )
        else:
            return self.root.right( interval.end, num_intervals, max_dist )

    def traverse(self, fn):
        """
        call fn for each element in the tree
        """
        if self.root is None:
            return None
        return self.root.traverse(fn)

    @classmethod
    def from_tuples(cls, tuples):
        """
        Create a new IntervalTree from an iterable of 2- or 3-tuples,
         where the tuple lists begin, end, and optionally data.
        """

        tree = IntervalTree()
        for iv in tuples:
            tree.insert_interval(Interval(*iv))
        return tree

    @classmethod
    def from_intervals(cls, intervals):
        """
        Create a new IntervalTree from an iterable of Interval instances.
        """

        tree = IntervalTree()
        for iv in intervals:
            tree.insert_interval(iv)
        return tree

    def __iter__(self):
        return self.root.__iter__()

    def __len__(self):
        return self.size()

    def __eq__(self, other):
        if not isinstance(other, IntervalTree):
            return False

        if self.size() != other.size():
            return False
        equal = True
        found = False

        for iv in self:
            found = False
            for oiv in other:
                if iv.start == oiv.start and iv.end == oiv.end:
                    found = True
                    break
            if not found:
                equal = False
                break

        return equal

    cdef inline void set_ends(IntervalTree self):
        self.root.set_ends()

    def fuzzy_equal(self, other, fuzzymatch=0):
        if not isinstance(other, IntervalTree):
            return False

        if self.size() != other.size():
            return False
        equal = True
        found = False

        for iv in self:
            found = False
            for oiv in other:
                if iv.fuzzy_equal(oiv, fuzzymatch):
                    found = True
                    break
            if not found:
                equal = False
                break

        return equal

    def __hash__(self):
        return hash(tuple(list(iter(self))))


# For backward compatibility
Intersecter = IntervalTree