# coding: utf-8

"""This module is used to serialise BLAST objects into a database.
It consists of various different classes:

- Query
- Target
- Hit
- HSP

The module also contains helper functions such as mean().
"""
import collections
import queue
import os
import sqlalchemy
import gzip
import subprocess
from Bio import SeqIO
import operator
import sqlalchemy.exc
from sqlalchemy.ext.hybrid import hybrid_property  # hybrid_method
from sqlalchemy import Column, String, Integer, Float, ForeignKey, Index
from sqlalchemy.sql.schema import PrimaryKeyConstraint, UniqueConstraint
from sqlalchemy.orm import relationship, backref
from Bio.Blast.NCBIXML import parse as xparser
import io
from sqlalchemy import create_engine
from sqlalchemy.orm.session import sessionmaker
from mikado_lib.serializers.dbutils import dbBase
from mikado_lib.serializers.dbutils import connect
import logging
import threading
import multiprocessing
import time


def create_opener(filename):

    """
    Function to create the appropriate opener for a BLAST file.
    If a handle is given instead of a filename, the function returns the input immediately.

    :param filename
    :return:
    """

    if type(filename) in (gzip.GzipFile, io.TextIOWrapper):
        return filename
    elif type(filename) is not str or not os.path.exists(filename):
        raise OSError("Non-existent file: {0}".format(filename))

    if filename.endswith(".gz"):
        if filename.endswith(".xml.gz"):
            return gzip.open(filename, "rt")
        elif filename.endswith(".asn.gz"):
            zcat = subprocess.Popen(["zcat", filename], shell=False,
                                    stdout=subprocess.PIPE)  # I cannot seem to make it work with gzip.open
            blast_formatter = subprocess.Popen(['blast_formatter', '-outfmt', '5',
                                                '-archive', '-'], shell=False,
                                               stdin=zcat.stdout,
                                               stdout=subprocess.PIPE)
            return io.TextIOWrapper(blast_formatter.stdout, encoding="UTF-8")
    elif filename.endswith(".xml"):
        return open(filename)
    else:
        raise ValueError("Unrecognized file format: {0}".format(filename))


def check_beginning(handle, filename, previous_header):

    """
    Static method to check that the beginning of the XML file is actually correct.

    :return
    """

    exc = None
    header = []
    position = 0
    try:
        first_line = next(handle)
        position += len(first_line)
        first_line = first_line.strip()
        if first_line != '<?xml version="1.0"?>':
            exc = ValueError("Invalid header for {0}!\n\t{1}".format(filename, first_line))
            raise exc
        second_line = next(handle)
        position += len(second_line)
        second_line = second_line.strip()
        valid_schemas = [
            " ".join(['<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN"',
                     '"http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">']),
            '<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "NCBI_BlastOutput.dtd">'
        ]

        if second_line not in valid_schemas:
            exc = ValueError("Invalid XML type for {0}!\n\t{1}".format(filename, second_line))
            raise exc
        header = [first_line, second_line]
    except StopIteration:
        exc = OSError("Empty file: {0}".format(filename))
    except ValueError:
        pass

    if exc is not None:
        return handle, header, exc

    while True:
        line = next(handle)

        if "<Iteration>" in line:
            # handle.seek(position)
            break
        position += len(line)
        line = line.rstrip()
        if not line:
            exc = ValueError("Invalid header for {0}:\n\n{1}".format(
                filename,
                "\n".join(header)
            ))
            break
        if len(header) > 10**3:
            exc = ValueError("Abnormally long header ({0}) for {1}:\n\n{2}".format(
                len(header),
                filename,
                "\n".join(header)
            ))
            break
        header.append(line)

    if not any(iter(True if "BlastOutput" in x else False for x in header)):
        exc = ValueError("Invalid header for {0}:\n\n{1}".format(filename, "\n".join(header)))

    if previous_header is not None:
        checker = list(filter(lambda x: "BlastOutput_query" not in x, header))
        previous_header = list(filter(lambda x: "BlastOutput_query" not in x, previous_header))
        if checker != previous_header:
            exc = ValueError("BLAST XML header does not match for {0}".format(filename))

    return handle, header, exc


class _Merger(multiprocessing.Process):

    def __init__(self, filenames, header, other_queue, logger=None):
        multiprocessing.Process.__init__(self)
        self.queue = other_queue
        self.filenames = filenames
        self.header = header
        self.logger = logger

    def run(self):
        # self.logger.info("Merger running")
        print_header = True

        for filename in self.filenames:
            self.logger.debug("Begun {0}".format(filename))
            if print_header is True:
                # self.logger.info("Printing header")
                self.queue.put_nowait(self.header)
                print_header = False
            try:
                handle = create_opener(filename)
            except Exception as exc:
                self.logger.exception(exc)
                continue
            handle, header, exc = check_beginning(handle,
                                                  filename,
                                                  self.header)
            if exc is not None:
                self.logger.exception(exc)
                self.logger.error("Skipped {0}".format(filename))

                continue
            self.logger.debug("Finished parsing header for {0}".format(filename))
            lines = collections.deque()
            lines.append("<Iteration>")
            bo_found = False
            for line in handle:
                if line.strip() == "":
                    continue
                if "BlastOutput" in line:
                    bo_found = True
                    break
                lines.append(line.rstrip())

            if bo_found is False:
                exc = ValueError("{0} is an invalid XML file".format(filename))
                self.logger.exception(exc)
                continue

            self.logger.debug("Finished parsing lines for {0}".format(filename))
            self.queue.put_nowait(lines)
            # self.logger.info("Dispatched {0} lines".format(counter))
            self.logger.debug("Sent {0} lines for {1}".format(len(lines), filename))
            # We HAVE  to wait some seconds, otherwise the XML parser might miss the end of the file.
            time.sleep(2)

        # We HAVE  to wait some seconds, otherwise the XML parser might miss the end of the file.
        time.sleep(1)
        self.queue.put_nowait(["</BlastOutput_iterations>\n</BlastOutput>"])

        self.queue.put("Finished")
        return


class XMLMerger(threading.Thread):

    logger = logging.getLogger("blast_merger")
    logger.propagate = False
    logger.setLevel(logging.INFO)
    stream_handler = logging.StreamHandler()
    formatter = logging.Formatter("{asctime} - {levelname} - {message}", style='{')
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    def __init__(self, filenames, log_level=logging.WARNING, log=None):

        threading.Thread.__init__(self)
        self.logger.setLevel(log_level)
        if log is not None:
            self.file_handler = logging.FileHandler(log, "w")
            self.file_handler.setFormatter(self.formatter)
            self.logger.addHandler(self.file_handler)
            self.logger.removeHandler(self.stream_handler)

        self.filenames = collections.deque(filenames)
        self.event = threading.Event()
        manager = multiprocessing.Manager()
        self.queue = manager.Queue()
        self.lines = collections.deque()
        self.started = False
        self.finished = False
        if len(self.filenames) > 0:
            while True:
                handle, header, exc = check_beginning(create_opener(self.filenames[0]),
                                                      self.filenames[0],
                                                      None
                                                      )
                if exc is not None:
                    self.logger.exception(exc)
                    _ = self.filenames.popleft()
                    continue
                self.header = header
                break
        self.logger.debug("Header is {0} long".format(len(self.header)))
        if len(self.filenames) == 0:
            raise IndexError("No files left!")

        self.merger = _Merger(self.filenames, self.header, self.queue, logger=self.logger)
        self.start()

    def run(self):
        self.merger.start()

        self.logger.info("Reader started")
        while True:
            try:
                lines = self.queue.get_nowait()
            except queue.Empty:
                self.logger.debug("Queue was empty, waiting")
                time.sleep(0.01)
                continue
            if lines == "Finished":
                break
            self.logger.debug("Received {0} lines".format(len(lines)))
            self.lines.extend(lines)

        self.finished = True
        self.merger.join()
        return

    def read(self, size=None):
        # self.logger.info("Requested {0} bytes".format(size))
        total = ""
        while len(self.lines) == 0:
            if self.finished is True:
                break
        if len(self.lines) == 0:
            raise StopIteration
        if size is None:
            total = self.lines
        else:
            previous_val = None
            import sys
            while sys.getsizeof(total) < size and len(self.lines) > 0:
                val = self.lines.popleft()
                if val == "</BlastOutput_iterations>\n</BlastOutput>":
                    self.logger.info("Received termination lines")
                    if previous_val is not None:
                        self.logger.info("Previous val: {0}".format(previous_val[-2000:]))
                total += val
                if previous_val is not None:
                    previous_val += val
                else:
                    previous_val = val

        total = total.rstrip()
        return total

    def __next__(self):

        while len(self.lines) == 0:
            if self.finished is True:
                break
        if len(self.lines) == 0:
            # Thread is finished and deque is exhausted. Finished
            raise StopIteration
        return self.lines.popleft()

    def __iter__(self):
        return self

    def __exit__(self, *args):
        _ = args
        self.join()

    def __enter__(self):
        pass


def mean(l: list):
    """
    :param l: a list of integers
    :type l: list(int)

    Very simple function to calculate the mean of a numeric array.
    """

    if len(l) == 0:
        raise ZeroDivisionError
    return sum(l) / len(l)


def merge(intervals: [(int, int)]):
    """
    :param intervals: a list of integer duplexes
    :type intervals: list

    This function is used to merge together intervals, which have to be supplied as a list
    of duplexes - (start,stop). The function will then merge together overlapping tuples and
    return a list of non-overlapping tuples.
    If the list is composed by only one element, the function returns immediately.
    """

    # Assume tuple of the form (start,end)
    # And return 0- and 1-length intervals
    new_intervals = []
    for interval in intervals:
        new_intervals.append(tuple(sorted(interval)))

    intervals = new_intervals[:]
    if len(intervals) < 2:
        return intervals

    # Sort according to start, end
    intervals = sorted(intervals, key=operator.itemgetter(0, 1))
    final_list = [intervals[0]]

    for start, end in intervals[1:]:
        if start > final_list[-1][1]:
            final_list.append(tuple([start, end]))
        elif end > final_list[-1][1]:
            final_list[-1] = tuple([final_list[-1][0], end])
    return final_list


class Query(dbBase):
    """
    Very simple serialization class for Query objects.

    :return id: integer key
    :rtype id: int

    :return name: name of the queries
    :rtype name: str

    :return length: length of the queries
    :rtype length: int
    """

    __tablename__ = "query"
    query_id = Column(Integer, primary_key=True)
    query_name = Column(String(200), unique=True, index=True)
    query_length = Column(Integer, nullable=True)  # This so we can load data also from the orf class

    named_tup = collections.namedtuple("Query", ["query_id", "query_name", "query_length"])

    def __init__(self, name, length):
        self.query_name = name
        self.query_length = length

    def as_tuple(self):
        """Quick function to convert the SQLalchemy object into a named tuple with the same fields"""
        return self.named_tup(self.query_id, self.query_name, self.query_length)


class Target(dbBase):

    """
    Very simple serialization class for Target objects.
    """

    __tablename__ = "target"

    target_id = Column(Integer, primary_key=True)
    target_name = Column(String(200), unique=True, index=True)
    target_length = Column(Integer)
    named_tup = collections.namedtuple("Target", ["target_id", "target_name", "target_length"])

    def __init__(self, target_name, target_length):
        """
        Constructor method.
        :param target_name: name of the targets
        :type target_name: str

        :param target_length: length of the targets
        :type target_length: int
        """

        self.target_name = target_name
        self.target_length = target_length

    def as_tuple(self):
        """Quick function to convert the SQLalchemy object into a named tuple with the same fields"""
        return self.named_tup(self.target_id, self.target_name, self.target_length)


class Hsp(dbBase):
    """
    This class serializes and stores into the DB the various HSPs.
    It is directly connected to the Hit table, through the "hit_id" 
    reference key.
    The Hit reference can be accessed through the hit_object attribute;
    back-reference (Hit to Hsps) is given by the "hsps" attribute.
    
    Keys:

    :return hit_id: Reference key for the Hit table
    :rtype hit_id: int

    :return counter: It indicates the progressive number of the HSP for the hit
    :rtype counter: int

    :return query_hsp_start: Start position on the query
    :rtype query_hsp_start; int

    :return query_hsp_end: End position on the query
    :rtype query_hsp_end: int

    :return target_hsp_start: Start position on the target
    :rtype target_hsp_start: int

    :return target_hsp_end: End position on the target
    :rtype target_hsp_end: int

    :return hsp_evalue: Evalue of the HSP
    :rtype hsp_evalue: float

    :return hsp_bits: Bit-score of the HSP
    :rtype hsp_bits: float

    :return hsp_identity: Identity (in %) of the alignment
    :rtype  hsp_identity: float

    :return hsp_length: Length of the HSP
    :rtype hsp_length: int
    
    An HSP row has the following constraints:
    - Counter,hit_id must be unique (and are primary keys)
    - The combination ("Hit_id","query_hsp_start","query_hsp_end", "target_hsp_start", "target_hsp_end") must be unique

    Moreover, the following properties are also present:

    :return query_object: The referenced Query
    :rtype query_object: Query

    :return target_object: The reference Target
    :rtype target_object: Target

    """

    __tablename__ = "hsp"
    counter = Column(Integer)  # Indicates the number of the HSP inside the hit
    query_id = Column(Integer, ForeignKey(Query.query_id), unique=False)
    target_id = Column(Integer, ForeignKey(Target.target_id), unique=False)
    pk_constraint = PrimaryKeyConstraint("counter", "query_id", "target_id", name="hsp_constraint")
    query_index = Index("hsp_query_idx", "query_id", unique=False)
    target_index = Index("hsp_target_idx", "target_id", unique=False)
    combined_index = Index("hsp_combined_idx", "query_id", "target_id", unique=False)
    full_index = Index("hsp_full_idx", "counter", "query_id", "target_id", unique=True)
    query_hsp_start = Column(Integer)
    query_hsp_end = Column(Integer)
    target_hsp_start = Column(Integer)
    target_hsp_end = Column(Integer)
    uni_constraint = UniqueConstraint("query_id", "target_id", "query_hsp_start", "query_hsp_end", "target_hsp_start",
                                      "target_hsp_end")
    hsp_evalue = Column(Float)
    hsp_bits = Column(Float)
    hsp_identity = Column(Float)
    hsp_length = Column(Integer)

    query_object = relationship(Query, uselist=False)
    target_object = relationship(Target, uselist=False)

    __table_args__ = (pk_constraint, query_index, target_index, combined_index)

    def __init__(self, hsp, counter, query_id, target_id):

        """

        :param hsp: an hsp object from the serialized XML

        :param counter: integer which indicates the HSP position in the HSP list for the hit
        :type counter: int

        :param query_id: Foreign key for the Query table
        :type query_id: int

        :param target_id: Foreign key for the Target table
        :type target_id: int
        """

        self.counter = counter
        self.query_hsp_start = hsp.query_start
        self.query_hsp_end = hsp.query_end
        self.target_hsp_start = hsp.sbjct_start
        self.target_hsp_end = hsp.sbjct_end
        self.hsp_identity = float(hsp.identities) / hsp.align_length * 100
        self.hsp_length = hsp.align_length
        self.hsp_bits = hsp.bits
        self.hsp_evalue = hsp.expect
        self.query_id = query_id
        self.target_id = target_id

    def __str__(self):
        """Simple printing function."""
        line = [self.query, self.target, self.query_hsp_start, self.query_hsp_end, self.target_hsp_start,
                self.target_hsp_end, self.hsp_evalue]
        return "\t".join([str(x) for x in line])

    # @profile

    @classmethod
    def as_dict_static(cls, state_obj):
        """Method to return a dict representation of the object. Necessary for storing.
        This method returns a dictionary *without any attribute that requires joined data*.
        As a static method, it is useful to be used outside of the class.
        :param state_obj: an instance of the HSP class or a collections.namedtuple object from a direct query

        :rtype : dict
        """

        keys = [
            "query_hsp_start",
            "query_hsp_end",
            "target_hsp_start",
            "target_hsp_end",
            "hsp_evalue",
            "hsp_bits"
        ]

        state = dict().fromkeys(keys)
        for key in keys:
            state[key] = getattr(state_obj, key)
        return state

    def as_dict(self):
        """Method to return a dict representation of the object. Necessary for storing.
        This method returns a dictionary *without any attribute that requires joined data*.
        It is meant to be used only by the method as_dict of the Hit class."""

        return self.as_dict_static(self)

    def as_full_dict(self):
        """Method to return a dict representation of the object.
        This method also checks query name and hit name, so it is slower than as_dict and used
        when it is necessary to retrieve data independently from Hit.
        """

        state = self.as_dict()
        state["query"] = self.query
        state["target"] = self.target
        state["query_hsp_cov"] = self.query_hsp_cov
        state["target_hsp_cov"] = self.target_hsp_cov

        return state

    @hybrid_property
    def query(self):
        """Returns the name of the query sequence, through a nested SQL query."""
        return self.query_object.query_name

    @hybrid_property
    def target(self):
        """Returns the name of the target sequence, through a nested SQL query."""
        return self.target_object.target_name

    @hybrid_property
    def query_hsp_cov(self):
        """This property returns the percentage of the query which is covered by the HSP."""
        return (self.query_hsp_end - self.query_hsp_start + 1) / self.query_object.query_length

    @hybrid_property
    def target_hsp_cov(self):
        """This property returns the percentage of the target which is covered by the HSP."""
        return (self.target_hsp_end - self.target_hsp_start + 1) / self.target_object.target_length


class Hit(dbBase):
    """This class is used to serialise and store in a DB a BLAST hit.
    Stored attributes:
    
    - id                Indexing key
    - query_id            Foreign ID key for the query table
    - target_id            Foreign ID key for the target table
    - qt_constrating    
    
    """

    __tablename__ = "hit"
    query_id = Column(Integer, ForeignKey(Query.query_id), unique=False)
    target_id = Column(Integer, ForeignKey(Target.target_id), unique=False)
    qt_constraint = PrimaryKeyConstraint("query_id", "target_id", name="hit_id")
    qt_index = Index("qt_index", "query_id", "target_id", unique=True)
    query_index = Index("hit_query_idx", "query_id", unique=False)
    target_index = Index("hit_target_idx", "target_id", unique=False)
    evalue = Column(Float)
    bits = Column(Float)
    global_identity = Column(Float)
    query_start = Column(Integer)
    query_end = Column(Integer)
    target_start = Column(Integer)
    target_end = Column(Integer)
    hit_number = Column(Integer)
    query_multiplier = Column(Float)  # Probably I should move this to a separate table!
    target_multiplier = Column(Float)
    query_aligned_length = Column(Integer)
    target_aligned_length = Column(Integer)

    query_object = relationship(Query, uselist=False,
                                lazy="immediate",
                                backref=backref("hits"))
    target_object = relationship(Target,
                                 uselist=False,
                                 lazy="immediate",
                                 backref=backref("hits"))

    hsps = relationship(Hsp, uselist=True,
                        # lazy="immediate",
                        backref=backref("hit_object", uselist=False),
                        foreign_keys=[query_id, target_id],
                        primaryjoin="and_(Hit.query_id==Hsp.query_id, Hit.target_id==Hsp.target_id)")

    __table_args__ = (qt_constraint,qt_index,query_index,target_index)

    def __init__(self, query_id, target_id, alignment, evalue, bits, hit_number=1, query_multiplier=1,
                 target_multiplier=1):
        """This function takes as input the id of a target, the id of the query, and a hit-object from the XML.
        The multiplier keyword is used to calculate the ratio between the query and the target.

        :param query_id: reference key for the Query table
        :type query_id: int

        :param target_id: reference key for the Target table
        :type target_id: int

        :param alignment: The BLAST Hit object

        :param evalue: Evalue of the hit (recovered from the "description" object)
        :type evalue: float

        :param bits: BitScore of the hit (recovered from the "description" object)
        :type bits: float

        :param hit_number: progressive index that indicates the priority of the hit in the database
        :type hit_number: int

        :param query_multiplier: either 1 or 3, it depends on the type of BLAST
        :type query_multiplier: int

        :param target_multiplier: either 1 or 3, it depends on the type of BLAST
        :type target_multiplier: int
        """

        self.query_id = query_id
        self.target_id = target_id
        self.query_multiplier = query_multiplier
        self.target_multiplier = target_multiplier
        self.hit_number = hit_number
        self.evalue = evalue
        self.bits = bits

        self.global_identity = mean([hsp.identities / hsp.align_length * 100 for hsp in alignment.hsps])
        q_intervals = [tuple([hsp.query_start, hsp.query_end]) for hsp in alignment.hsps]
        q_merged_intervals = sorted(merge(q_intervals), key=operator.itemgetter(0, 1))
        q_aligned = sum([tup[1] - tup[0] + 1 for tup in q_merged_intervals])
        self.query_aligned_length = q_aligned
        self.query_start = q_merged_intervals[0][0]
        self.query_end = q_merged_intervals[-1][1]

        t_intervals = [tuple([hsp.sbjct_start, hsp.sbjct_end]) for hsp in alignment.hsps]
        t_merged_intervals = sorted(merge(t_intervals), key=operator.itemgetter(0, 1))
        t_aligned = sum([tup[1] - tup[0] + 1 for tup in t_merged_intervals])
        self.target_aligned_length = t_aligned
        self.target_start = t_merged_intervals[0][0]
        self.target_end = t_merged_intervals[-1][1]

    def __str__(self):
        line = [self.query, self.target, self.evalue, self.bits, self.query_start, self.query_end, self.target_start,
                self.target_end, self.query_len, self.target_len]

        return "\t".join(str(x) for x in line)

    # @profile

    @classmethod
    def as_dict_static(cls, state_obj):
        """Method to return a dict representation of the object.
        Static method to be called from outside the class.

        :param state_obj: a namedtuple or an instance of this class

        :rtype: dict
        """

        keys = [
            "evalue",
            "bits",
            "global_identity",
            "query_start",
            "query_end",
            "target_start",
            "target_end",
            "hit_number",
            "query_multiplier",
            "target_multiplier",
            "query_aligned_length",
            "target_aligned_length",
        ]

        state = dict().fromkeys(keys)

        for key in keys:
            state[key] = getattr(state_obj, key)

        return state

    @classmethod
    def as_full_dict_static(cls, hit_tuple,
                            hsp_list,
                            query_tuple,
                            target_tuple):
        """
        :param hit_tuple: Hit namedtuple (from direct query to the DB)
        :type hit_tuple: collections.namedtuple
        :param hsp_list: list of hsp dictionaries from Hsp.as_dict_static
        :type hsp_list: list(collections.namedtuple)

        :param query_tuple: Query namedtuple
        :type query_tuple: collections.namedtuple
        :param target_tuple: Target namedtuple
        :type target_tuple: collections.namedtuple
        :rtype: dict
        """

        state = cls.as_dict_static(hit_tuple)
        hsps = [Hsp.as_dict_static(h) for h in hsp_list]

        state["query"] = query_tuple.query_name
        state["target"] = target_tuple.target_name
        state["query_len"] = query_tuple.query_length
        state["target_len"] = target_tuple.target_length
        state["query_hit_ratio"] = (query_tuple.query_length * hit_tuple.query_multiplier) / \
                                   (target_tuple.target_length * hit_tuple.target_multiplier)
        state["hit_query_ratio"] = (target_tuple.target_length * hit_tuple.target_multiplier) / \
                                   (query_tuple.query_length * hit_tuple.query_multiplier)
        state["hsps"] = []
        for h in hsps:
            h["query_hsp_cov"] = (h["query_hsp_end"] - h["query_hsp_start"] + 1) / (state["query_len"])
            h["target_hsp_cov"] = (h["target_hsp_end"] - h["target_hsp_start"] + 1) / (state["target_len"])
            state["hsps"].append(h)

        return state

    def as_dict(self):
        """Method to return a dict representation of the object.
        Necessary for storing.

        :rtype: dict
        """

        state = self.as_dict_static(self)

        # Retrieving the values ONCE
        query_object = self.query_object.as_tuple()
        target_object = self.target_object.as_tuple()

        state["query"] = query_object.query_name
        state["target"] = target_object.target_name
        state["query_len"] = query_object.query_length
        state["target_len"] = target_object.target_length
        state["query_hit_ratio"] = state["query_len"] * state["query_multiplier"] /\
            (state["target_len"] * state["target_multiplier"])

        state["hit_query_ratio"] = state["target_len"] * state["target_multiplier"] /\
            (state["query_len"] * state["query_multiplier"])

        state["hsps"] = []
        for hsp in self.hsps:
            h = hsp.as_dict()
            h["query_hsp_cov"] = (h["query_hsp_end"] - h["query_hsp_start"] + 1) / (state["query_len"])
            h["target_hsp_cov"] = (h["target_hsp_end"] - h["target_hsp_start"] + 1) / (state["target_len"])
            state["hsps"].append(h)

        return state

    @hybrid_property
    def query_len(self):
        """
        This method retrieves the length of the query from the entry in the Query table.
        :rtype int
        """
        return self.query_object.query_length

    @hybrid_property
    def query(self):
        """
        This method retrieves the name of the query from the entry in the Query table.
        :rtype str
        """
        return self.query_object.query_name

    @hybrid_property
    def target(self):
        """
        This method retrieves the name of the target from the entry in the Target table.
        :rtype str
        """

        return self.target_object.target_name

    @hybrid_property
    def target_len(self):
        """
        This method retrieves the length of the target from the entry in the Target table.
        :rtype int
        """

        return self.target_object.target_length

    @hybrid_property
    def query_hit_ratio(self):
        """
        This property returns the quotient (Query Length)/(Target Length)
        """

        return self.query_len * self.query_multiplier / (self.target_len * self.target_multiplier)

    @hybrid_property
    def hit_query_ratio(self):
        """
        This property returns the quotient (Target Length)/(Query Length)
        """

        return self.target_len * self.target_multiplier / (self.query_len * self.query_multiplier)


class XmlSerializer:
    """This class has the role of taking in input a blast XML file and (partially) serialise it into
    a database. We are using SQLalchemy, so the database type could be any of SQLite, MySQL, PSQL, etc."""

    logger = logging.getLogger("blast_serialiser")
    logger.propagate = False
    logger.setLevel(logging.INFO)
    handler = logging.StreamHandler()
    formatter = logging.Formatter("{asctime} - {levelname} - {message}", style='{')
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    def __init__(self, xml, max_target_seqs=float("Inf"),
                 db=None,
                 target_seqs=None,
                 query_seqs=None,
                 discard_definition=True, maxobjects=10000,
                 json_conf=None
                 ):
        """Initializing method. Arguments:

        :param db: the name of the SQLite database
        :type db: str | None

        :param xml: The XML to parse.

        Optional arguments:

        :param target_seqs: either a BioPython index of a FASTA file or the file itself.

        :param query_seqs: either a BioPython index of a FASTA file or the file itself.

        :param discard_definition: flag. If set to True, the "id" field in the XML
        instead of "definition" will be used for serializing.
        :type discard_definition: bool

        :param maxobjects: maximum number of objects to keep in memory before bulk loading. Default: 10^4
        :type maxobjects: int

        :param json_conf: a configuration dictionary.
        :type json_conf: dict


        """

        if xml is None:
            self.logger.warning("No BLAST XML provided. Exiting.")
            return

        if json_conf is not None:
            self.engine = connect(json_conf)
        else:
            if db is None:
                db = ":memory:"
            self.engine = create_engine("sqlite:///{0}".format(db))

        session = sessionmaker()
        session.configure(bind=self.engine)
        dbBase.metadata.create_all(self.engine)  # @UndefinedVariable
        self.session = session()
        self.logger.debug("Created the session")
        if type(xml) is str:
            self.xml_parser = xparser(create_opener(xml))
        else:
            assert type(xml) in (list, set)
            self.xml_parser = XMLMerger(xml)  # Merge in memory

        # Runtime arguments
        self.discard_definition = discard_definition

        if type(query_seqs) is str:
            assert os.path.exists(query_seqs)
            self.query_seqs = SeqIO.index(query_seqs, "fasta")
        elif query_seqs is None:
            self.query_seqs = None
        else:
            assert "SeqIO.index" in repr(query_seqs)
            self.query_seqs = query_seqs

        if type(target_seqs) is str:
            assert os.path.exists(target_seqs)
            self.target_seqs = SeqIO.index(target_seqs, "fasta")
        elif target_seqs is None:
            self.target_seqs = None
        else:
            assert "SeqIO.index" in repr(target_seqs)
            self.target_seqs = target_seqs

        self.max_target_seqs = max_target_seqs
        self.maxobjects = maxobjects

    def serialize(self):

        """Method to serialize the BLAST XML file into a database provided with the __init__ method """

        q_mult = 1  # Query multiplier: 3 for BLASTX (nucleotide vs. protein), 1 otherwise
        h_mult = 1  # Target multiplier: 3 for TBLASTN (protein vs. nucleotide), 1 otherwise

        objects = []

        targets = dict()
        queries = dict()
        self.logger.info("Loading previous IDs")
        for query in self.session.query(Query):
            queries[query.query_name] = query.query_id
        for query in self.session.query(Target):
            targets[query.target_name] = query.target_id
        self.logger.info("Loaded previous IDs")

        self.logger.info("Started the serialisation")
        if self.target_seqs is not None:
            self.logger.info("Started to serialise the targets")
            for record in self.target_seqs:
                if record in targets:
                    continue
                objects.append(Target(record, len(self.target_seqs[record])))
                if len(objects) >= self.maxobjects:
                    self.logger.info("Loading {0} objects into the \"target\" table".format(len(objects)))
                    self.engine.execute(Target.__table__.insert(),
                                        [{"target_name": obj.target_name,
                                          "target_length": obj.target_length} for obj in objects]
                                        )
                    #                         self.session.bulk_save_objects(objects, return_defaults=False)
                    self.logger.info("Loaded {0} objects into the \"target\" table".format(len(objects)))
                    objects = []
            self.logger.info("Loading {0} objects into the \"target\" table".format(len(objects)))
            self.engine.execute(Target.__table__.insert(),
                                [{"target_name": obj.target_name,
                                  "target_length": obj.target_length} for obj in objects]
                                )
            #             self.session.bulk_save_objects(objects)
            self.logger.info("Loaded {0} objects into the \"target\" table".format(len(objects)))
            objects = []
            self.logger.info("Loaded targets")
            self.session.commit()

        if self.query_seqs is not None:
            self.logger.info("Started to serialise the queries")
            for record in self.query_seqs:
                if record in queries:
                    continue
                objects.append(Query(record, len(self.query_seqs[record])))
                if len(objects) >= self.maxobjects:
                    self.logger.info("Loading {0} objects into the \"query\" table".format(len(objects)))
                    self.engine.execute(Query.__table__.insert(),
                                        [{"query_name": obj.query_name,
                                          "query_length": obj.query_length} for obj in objects]
                                        )
                    #                     self.session.bulk_save_objects(objects, return_defaults=False)
                    self.logger.info("Loaded {0} objects into the \"query\" table".format(len(objects)))
                    objects = []
            self.logger.info("Loading {0} objects into the \"query\" table".format(len(objects)))
            self.engine.execute(Query.__table__.insert(),
                                [{"query_name": obj.query_name, "query_length": obj.query_length} for obj in objects]
                                )
            #             self.session.bulk_save_objects(objects)
            self.session.commit()
            self.logger.info("Loaded {0} objects into the \"query\" table".format(len(objects)))
            objects = []
            self.logger.info("Queries serialised")

        self.logger.info("Loading all IDs")
        for query in self.session.query(Query):
            queries[query.query_name] = (query.query_id, query.query_length is not None)
        for query in self.session.query(Target):
            targets[query.target_name] = (query.target_id, query.target_length is not None)
        self.logger.info("Loaded all IDs")

        # Memorize the mapping ... it is just faster

        query_counter = 0
        for record in self.xml_parser:
            if record.application in ("BLASTN", "TBLASTX", "BLASTP"):
                q_mult = 1
                h_mult = 1
            elif record.application == "BLASTX":
                q_mult = 3
                h_mult = 1
            elif record.application == "TBLASTN":
                q_mult = 1
                h_mult = 3
            if len(record.descriptions) == 0:
                continue
            query_counter += 1
            if self.discard_definition is False:
                name = record.query.split()[0]
            else:
                name = record.query_id
            self.logger.debug("Started with {0}".format(name))

            if name in queries:
                current_query = queries[name][0]
                if queries[name][1] is False:
                    self.session.query(Query).filter(Query.query_name == name).update(
                        {"query_length": record.query_length})
                    self.session.commit()
            else:
                self.logger.warn("Adding {0} to the db".format(name))
                current_query = Query(name, record.query_length)
                self.session.add(current_query)
                self.session.commit()
                queries[name] = (current_query.query_id, True)
                current_query = current_query.query_id

            for ccc, alignment in filter(lambda x: x[0] <= self.max_target_seqs, enumerate(record.alignments)):

                self.logger.debug("Started the hit {0}-{1}".format(name, record.alignments[ccc].accession))
                evalue = record.descriptions[ccc].e
                bits = record.descriptions[ccc].bits
                alignment = record.alignments[ccc]
                hit_num = ccc + 1
                if record.alignments[ccc].accession in targets:
                    current_target = targets[record.alignments[ccc].accession][0]
                    if targets[record.alignments[ccc].accession][1] is False:
                        self.session.query(Target).filter(Target.target_name == record.alignments[ccc].accession).\
                            update({"target_length": record.query_length})

                else:
                    current_target = Target(record.alignments[ccc].accession, record.alignments[ccc].length)
                    self.session.add(current_target)
                    try:
                        self.session.commit()
                        assert type(current_target.target_id) is int
                        targets[record.alignments[ccc].accession] = (current_target.target_id, True)
                        current_target = current_target.target_id
                    except sqlalchemy.exc.IntegrityError:
                        self.session.rollback()
                        continue

                current_hit = Hit(current_query, current_target, alignment, evalue, bits,
                                  hit_number=hit_num,
                                  query_multiplier=q_mult,
                                  target_multiplier=h_mult)
                objects.append(current_hit)

                for counter, hsp in enumerate(alignment.hsps):
                    current_hsp = Hsp(hsp, counter, current_query, current_target)
                    objects.append(current_hsp)

            if len(objects) >= self.maxobjects:
                self.logger.info("Loading {0} objects into the hit, hsp tables; done {1} queries".format(len(objects),
                                                                                                         query_counter))
                self.load_into_db(objects)
                self.logger.info(
                    "Loaded {0} objects into the hit, hsp tables; done {1} queries".format(len(objects), query_counter))
                objects = []

        self.logger.info("Loading {0} objects into the hit, hsp tables".format(len(objects)))
        self.load_into_db(objects)

        self.logger.info("Loaded {0} objects into the hit, hsp tables".format(len(objects)))
        self.session.commit()
        self.logger.info("Finished loading blast hits")

    def __call__(self):
        """
        Alias for serialize
        """
        self.serialize()

    def load_into_db(self, objects):
        """
        :param objects: Objects to be loaded into the database
        :type objects: list

        Method to perform the bulk loading of objects into the SQL database.

        """

        try:
            self.session.bulk_save_objects(objects, return_defaults=False)
            self.session.commit()
        except Exception as err:
            self.logger.error('Database corrupted'.format(err))
            self.logger.error(err)
            self.logger.error('Dropping and reloading')
            self.session.rollback()
            self.session.query(Hsp).delete()
            self.session.query(Hit).delete()
            self.session.bulk_save_objects(objects, return_defaults=False)
            self.session.commit()
