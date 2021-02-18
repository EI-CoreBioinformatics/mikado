#!/usr/bin/python2.5
"""
Intersects ... faster.  Suports GenomicInterval datatype and multiple
chromosomes.
"""
import operator

#@cdef extern from "math.h":
cdef extern from "stdlib.h":
    int ceil(float f)
    float log(float f)
    int RAND_MAX
    int rand()
    int strlen(char *)
    int iabs(int)

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

    def __init__(self, int start, int end, object value=None, object chrom=None, object strand=None):
        assert start <= end, "start must be less than end"
        self.start  = start
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

cpdef int distance(Interval f1, Interval f2):
    """\
    Distance between 2 features. The integer result is always positive or zero.
    If the features overlap or touch, it is zero.
    >>> from quicksect import Interval, distance
    >>> distance(Interval(1, 2), Interval(12, 13))
    10
    >>> distance(Interval(1, 2), Interval(2, 3))
    0
    >>> distance(Interval(1, 100), Interval(20, 30))
    0

    """
    if f1.end < f2.start: return f2.start - f1.end
    if f2.end < f1.start: return f1.start - f2.end
    return 0


cdef class IntervalTree:
    # cdef IntervalNode root

    def __init__(self):
        self.root = None
        self.num_intervals = 0

    # ---- Pickling ------

    def __getstate__(self):

        return [self.root, self.num_intervals]

    def __setstate__(self, state):

        self.root, self.num_intervals = state

    cpdef insert_interval(self, Interval interval):

        if self.root is None:
            self.root = IntervalNode(interval)
        else:
            self.root = self.root.insert(interval)
        self.num_intervals += 1

    add_interval = insert_interval

    cpdef insert(self, int start, int end, value=None):
        return self.insert_interval(Interval(start, end, value))

    add = insert

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
            tree.insert(*iv)
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

    cpdef find(self, int start, int end, bint strict=0, bint contained_check=0, int max_distance=0,
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
            found.extend(self.left(Interval(start, end), num_intervals=num_intervals, max_dist=max_distance))
            found.extend(self.right(Interval(start, end), num_intervals=num_intervals, max_dist=max_distance))
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

    def left(self, Interval f, int num_intervals=1, int max_dist=25000):
        if self.root is None:
            return []
        else:
            return self.root.left(f, num_intervals, max_dist)

    def right(self, Interval f, int num_intervals=1, int max_dist=25000):
        if self.root is None:
            return []
        else:
            return self.root.right(f, num_intervals, max_dist)

    def dump(self, fn):
        try:
            import cPickle
        except ImportError:
            import pickle as cPickle
        l = []
        a = l.append
        self.root.traverse(a)
        fh = open(fn, "wb")
        for f in l:
            cPickle.dump(f, fh)

    def load(self, fn):
        try:
            import cPickle
        except ImportError:
            import pickle as cPickle
        fh = open(fn, "rb")
        while True:
            try:
                feature = cPickle.load(fh)
                self.insert_interval(feature)
            except EOFError:
                break

    cpdef int size(self):
        return self.num_intervals

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
    """\
    Data structure for performing intersect and neighbor queries on a 
    set of intervals. Algorithm uses a segment/interval tree to perform
    efficient queries. 

    Usage
    =====
    >>> from quicksect import IntervalNode, Interval
    >>> tree = IntervalNode(Interval(0, 10))

    Add intervals, the only requirement is that the interval have integer
    start and end attributes. Optional arguments are strand, name, and info.

    >>> Interval(1, 22, info={'chr':12, 'anno': 'anything'})


    >>> tree = tree.insert(Interval(3, 7, 1))
    >>> tree = tree.insert(Interval(3, 40, -1))
    >>> tree = tree.insert(Interval(13, 50, 1))

    Queries
    -------

    find
    ++++

    >>> tree.find(2, 5)
    [Interval(3, 7), Interval(3, 40), Interval(0, 10)]
    >>> tree.find(11, 100)
    [Interval(13, 50), Interval(3, 40)]
    >>> tree.find(100, 200)
    []

    left/right
    ++++++++++
    the left method finds features that are strictly to the left of
    the query feature. overlapping features are not considered:

    >>> tree.left(Interval(0, 1))
    []
    >>> tree.left(Interval(11, 12))
    [Interval(0, 10)]

    """
    # cdef float priority
    # cdef public Interval interval
    # cdef public int start, end
    # cdef int minstop, maxstop, minstart
    # cdef IntervalNode cleft, cright, croot

    property left_node:
        def __get__(self):
            return self.cleft if self.cleft is not EmptyNode else None
    property right_node:
        def __get__(self):
            return self.cright if self.cright is not EmptyNode else None
    property root_node:
        def __get__(self):
            return self.croot if self.croot is not EmptyNode else None

    def __repr__(self):
        return "IntervalNode(%i, %i)" % (self.start, self.end)

    def __cinit__(IntervalNode self, Interval interval):
        # Python lacks the binomial distribution, so we convert a
        # uniform into a binomial because it naturally scales with
        # tree size.  Also, python's uniform is perfect since the
        # upper limit is not inclusive, which gives us undefined here.
        self.priority   = ceil(nlog * log(-1.0/(1.0 * rand()/RAND_MAX - 1)))
        self.start      = interval.start
        self.end       = interval.end
        self.interval   = interval
        self.maxstop    = interval.end
        self.minstart   = interval.start
        self.minstop    = interval.end
        self.cleft       = EmptyNode
        self.cright      = EmptyNode
        self.croot       = EmptyNode

    def insert(self, interval):
        return self._insert(interval)

    cdef IntervalNode _insert(IntervalNode self, Interval interval):
        cdef IntervalNode croot = self
        if interval.start > self.start:

            # insert to cright tree
            if self.cright is not EmptyNode:
                self.cright = self.cright._insert(interval )
            else:
                self.cright = IntervalNode(interval)
            # rebalance tree
            if self.priority < self.cright.priority:
                croot = self.rotate_left()

        else:
            # insert to cleft tree
            if self.cleft is not EmptyNode:
                self.cleft = self.cleft._insert(interval)
            else:
                self.cleft = IntervalNode(interval)
            # rebalance tree
            if self.priority < self.cleft.priority:
                croot = self.rotate_right()
    
        croot.set_stops()
        self.cleft.croot  = croot
        self.cright.croot = croot
        return croot

    cdef IntervalNode rotate_right(IntervalNode self):
        cdef IntervalNode croot = self.cleft
        self.cleft  = self.cleft.cright
        croot.cright = self
        self.set_stops()
        return croot

    cdef IntervalNode rotate_left(IntervalNode self):
        cdef IntervalNode croot = self.cright
        self.cright = self.cright.cleft
        croot.cleft  = self
        self.set_stops()
        return croot

    cdef inline void set_stops(IntervalNode self):
        if self.cright is not EmptyNode and self.cleft is not EmptyNode: 
            self.maxstop = imax3(self.end, self.cright.maxstop, self.cleft.maxstop)
            self.minstop = imin3(self.end, self.cright.minstop, self.cleft.minstop)
            self.minstart = imin3(self.start, self.cright.minstart, self.cleft.minstart)
        elif self.cright is not EmptyNode:
            self.maxstop = imax2(self.end, self.cright.maxstop)
            self.minstop = imin2(self.end, self.cright.minstop)
            self.minstart = imin2(self.start, self.cright.minstart)
        elif self.cleft is not EmptyNode:
            self.maxstop = imax2(self.end, self.cleft.maxstop)
            self.minstop = imin2(self.end, self.cleft.minstop)
            self.minstart = imin2(self.start, self.cleft.minstart)
        

    cpdef intersect(self, int start, int stop):
        """
        given a start and a stop, return a list of features
        falling within that range
        """
        cdef list results = []
        self._intersect(start, stop, results)
        return results

    find = intersect
        
    cdef void _intersect(IntervalNode self, int start, int stop, list results):
        # to have starts, stops be non-inclusive, replace <= with <  and >= with >
        #if start <= self.end and stop >= self.start: results.append(self.interval)
        if (not self.end < start) and (not self.start > stop): results.append(self.interval)
        #if self.cleft is not EmptyNode and start <= self.cleft.maxstop:
        if self.cleft is not EmptyNode and not self.cleft.maxstop < start:
            self.cleft._intersect(start, stop, results)
        #if self.cright is not EmptyNode and stop >= self.start:
        if self.cright is not EmptyNode and not self.start > stop:
            self.cright._intersect(start, stop, results)

    cdef void _seek_left(IntervalNode self, int position, list results, int n, int max_dist):
        # we know we can bail in these 2 cases.
        if self.maxstop + max_dist < position: return
        if self.minstart > position: return

        #import sys;sys.stderr.write( " ".join(map(str, ["SEEK_LEFT:", self, self.cleft, self.maxstop, self.minstart,  position])))

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
        if self.maxstop < position: return
        if self.minstart - max_dist > position: return

        #print "SEEK_RIGHT:",self, self.cleft, self.maxstop, self.minstart, position

        # the ordering of these 3 blocks makes it so the results are
        # ordered nearest to farest from the query position
        if self.cleft is not EmptyNode: 
                self.cleft._seek_right(position, results, n, max_dist)

        if -1 < self.start - position < max_dist:
            results.append(self.interval)

        if self.cright is not EmptyNode:
                self.cright._seek_right(position, results, n, max_dist)

    def neighbors(self, Interval f, int n=1, int max_dist=25000):
        cdef list neighbors = []

        cdef IntervalNode right = self.cright
        while right.cleft is not EmptyNode:
            right = right.cleft

        cdef IntervalNode left = self.cleft
        while left.cright is not EmptyNode:
            left = left.cright
        return [left, right]

    cpdef left(self, Interval f, int num_intervals=1, int max_dist=25000):
        """find n features with a start > than f.end
        f: a Interval object
        num_intervals: the number of features to return
        max_dist: the maximum distance to look before giving up.
        """
        cdef list results = []
        # use start - 1 becuase .left() assumes strictly left-of
        self._seek_left(f.start - 1, results, num_intervals, max_dist)
        if len(results) <= num_intervals: return results
        r = results
        r.sort(key=operator.attrgetter('end'), reverse=True)
        if distance(f, r[num_intervals]) != distance(f, r[num_intervals-1]):
            return r[:num_intervals]
        while num_intervals < len(r) and distance(r[num_intervals], f) == distance(r[num_intervals - 1], f):
            num_intervals += 1
        return r[:num_intervals]

    cpdef right(self, Interval f, int num_intervals=1, int max_dist=25000):
        """find n features with a stop < than f.start
        f: a Interval object
        num_intervals: the number of features to return
        max_dist: the maximum distance to look before giving up.
        """
        cdef list results = []
        # use stop + 1 becuase .right() assumes strictly right-of
        self._seek_right(f.end + 1, results, num_intervals, max_dist)
        if len(results) <= num_intervals: return results
        r = results
        r.sort(key=operator.attrgetter('start'))
        if distance(f, r[num_intervals]) != distance(f, r[num_intervals-1]):
            return r[:num_intervals]
        while num_intervals < len(r) and distance(r[num_intervals], f) == distance(r[num_intervals - 1], f):
            num_intervals += 1
        return r[:num_intervals]

    def __iter__(self):
            
        if self.cleft is not EmptyNode:
            yield self.cleft

        yield self.interval

        if self.cright is not EmptyNode:
            yield self.cright

    def traverse(self, func):
        self._traverse(func)

    cdef void _traverse(IntervalNode self, object func):
        if self.cleft is not EmptyNode: self.cleft._traverse(func)
        func(self.interval)
        if self.cright is not EmptyNode: self.cright._traverse(func)

    cpdef tuple _as_tuple(self):
        return (self.start, self.end)

    def __hash__(self):
        return hash(self._as_tuple())

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



cdef IntervalNode EmptyNode = IntervalNode(Interval(0, 0))
