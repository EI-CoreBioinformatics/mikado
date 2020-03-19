"""
Generic utilities used for BLAST serialising into a DB.
"""

from ...parsers.blast_utils import merge
import numpy as np
from .aln_string_parser import prepare_aln_strings

# valid_matches = set([chr(x) for x in range(65, 91)] + [chr(x) for x in range(97, 123)] +
#                     ["|", "*"])


def get_off_by_one(record):
    """"""
    ref = getattr(record, "reference", None)
    if ref is None:
        ref = record.__dict__.get("reference", None)
    if ref is None:
        return False
    elif isinstance(ref, str) and "diamond":
        version = record.version
        try:
            macro, minor, micro = [int(_) for _ in version.split(".")]
        except (ValueError, TypeError):
            return False
        if macro > 0 or minor > 9 or micro > 30:
            return False
        else:
            return True


def get_multipliers(record, application=None):
    """
    Private quick method to determine the multipliers for a BLAST alignment
    according to the application present in the record.
    :param record:
    :type record: Bio.SearchIO.Record
    :return:
    """

    q_mult, h_mult = 1, 1

    if record is not None:
        application = record.program.upper()
    else:
        application = application.upper()
    if application in ("BLASTN", "TBLASTX", "BLASTP"):
        q_mult = 1
        h_mult = 1
    elif application == "BLASTX":
        q_mult = 3
        h_mult = 1
    elif application == "TBLASTN":
        q_mult = 1
        h_mult = 3

    return q_mult, h_mult


def _get_query_for_blast(record, cache):
    """ This private method formats the name of the query
    recovered from the BLAST hit. It will cause an exception if the target is not
    present in the dictionary.
    :param record: Bio.SearchIO.Record
    :return: current_query (ID in the database), name
    """

    if record.id in cache:
        return cache[record.id], record.id, cache
    elif record.id.split()[0] in cache:
        return cache[record.id.split()[0]], record.id.split()[0], cache
    else:
        raise ValueError("{} not found (Accession: {})".format(record.id, record.__dict__))


def _get_target_for_blast(alignment, cache):
    """ This private method retrieves the correct target_id
    key for the target of the BLAST. If the entry is not present
    in the database, it will be created on the fly.
    The method returns the index of the current target and
    and an updated target dictionary.
    :param alignment: an alignment child of a BLAST record object
    :type alignment: Bio.SearchIO.Hit
    :return: current_target (ID in the database), targets
    """

    if alignment.accession in cache:
        return cache[alignment.accession], cache
    elif alignment.id in cache:
        return cache[alignment.id], cache
    else:
        raise ValueError("{} not found (Accession: {})".format(alignment.id, alignment.__dict__))


def prepare_hsp(hsp, counter, off_by_one=False, qmultiplier=1, tmultiplier=1):

    r"""
    Prepare a HSP for loading into the DB.
    The match line will be reworked in the following way:

    - If the position is a match/positive, keep the original value
    - If the position is a gap *for the query*, insert a - (dash)
    - If the position is a gap *for the target*, insert a _ (underscore)
    - If the position is a gap *for both*, insert a \ (backslash)

    :param hsp: An HSP object from Bio.Blast.NCBIXML
    :type hsp: Bio.SearchIO.HSP
    :param counter: a digit that indicates the priority of the HSP in the hit
    :return: hsp_dict, identical_positions, positives
    :rtype: (dict, set, set)
    """

    hsp_dict = dict()
    # We must start from 1, otherwise MySQL crashes as its indices start from 1 not 0
    hsp_dict["counter"] = counter + 1
    hsp_dict["query_hsp_start"] = hsp.query_start
    match, identical_positions, positives = prepare_aln_strings(hsp, off_by_one=off_by_one,
                                                                qmultiplier=qmultiplier)
    hsp_dict["query_hsp_end"] = hsp.query_end + off_by_one
    hsp_dict["query_frame"] = hsp.query_frame
    hsp_dict["target_hsp_start"] = hsp.hit_start
    hsp_dict["target_hsp_end"] = hsp.hit_end
    hsp_dict["target_frame"] = hsp.hit_frame
    hsp_dict["hsp_identity"] = hsp.ident_num / hsp.aln_span * 100
    hsp_dict["hsp_positives"] = (match.count("+") + match.count("|")) / hsp.aln_span * 100
    hsp_dict["match"] = match
    hsp_dict["hsp_length"] = hsp.aln_span
    hsp_dict["hsp_bits"] = hsp.bitscore
    hsp_dict["hsp_evalue"] = hsp.evalue
    return hsp_dict, identical_positions, positives


def _np_grouper(data):
    return np.array(np.split(data, np.where(np.diff(data) != 1)[0] + 1))


def prepare_hit(hit, query_id, target_id, off_by_one=False, **kwargs):
    """Prepare the dictionary for fast loading of Hit and Hsp objects.
    global_positives: the similarity rate for the global hit *using the query perspective*
    global_identity: the identity rate for the global hit *using the query perspective*

    :param hit: the hit to parse.
    :type hit: Bio.SearchIO.Hit

    :param query_id: the numeric ID of the query in the database. Necessary for serialisation.
    :type query_id: int

    :param target_id: the numeric ID of the target in the database. Necessary for serialisation.
    :type target_id: int

    :param kwargs: additional properties to give to the hit_dict. Retrieved e.g. from descriptions.
    :type kwargs: dict
    """

    hit_dict = dict()
    hsp_dict_list = []
    q_intervals = []
    t_intervals = []

    qmulti = kwargs["query_multiplier"]
    tmulti = kwargs["target_multiplier"]
    qlength = kwargs["query_length"]
    assert isinstance(qmulti, (int, float)), type(qmulti)
    assert isinstance(tmulti, (int, float)), type(tmulti)
    hit_dict.update(kwargs)
    hit_dict["query_id"] = query_id
    hit_dict["target_id"] = target_id

    query_array = np.zeros([2, int(qlength)], dtype=np.int)

    for counter, hsp in enumerate(hit.hsps):
        if hsp.query_start + off_by_one - 1 > qlength:
            raise ValueError("Invalid length: {}, {}", hsp.query_start + off_by_one - 1, qlength)

        hsp_dict, ident, posit = prepare_hsp(hsp, counter, qmultiplier=qmulti, tmultiplier=tmulti,
                                             off_by_one=off_by_one)
        if ident.max() > query_array.shape[1] or ident.min() < 0:
            raise IndexError("Invalid indexing values (max {}; frame {}; hsp: {})!"\
"Too low: {}\nToo high: {}".format(
                query_array.shape[1],
                hsp.query_frame, hsp.__dict__,
                list(ident[(ident < 0)]),
                list(ident[ident > query_array.shape[1]])))
        try:
            query_array[0, ident] = 1
        except IndexError as exc:
            raise IndexError("{}, off by one: {}; min, max: {}, {}; hsp {}".format(
                exc, off_by_one, ident.min(), ident.max(), hsp._items[0].__dict__))
        try:
            query_array[1, posit] = 1
        except IndexError as exc:
            raise IndexError("{}, off by one: {}; min, max: {}, {}; frame: {}".format(
                exc, off_by_one, posit.min(), posit.max(), hsp.query_frame))
        hsp_dict["query_id"] = query_id
        hsp_dict["target_id"] = target_id
        hsp_dict_list.append(hsp_dict)
        q_intervals.append((hsp.query_start, hsp.query_end - 1))
        t_intervals.append((hsp.hit_start, hsp.hit_end - 1))

    q_merged_intervals, q_aligned = merge(q_intervals)
    hit_dict["query_aligned_length"] = min(qlength, q_aligned)
    qstart, qend = q_merged_intervals[0][0], q_merged_intervals[-1][1]
    hit_dict["query_start"], hit_dict["query_end"] = qstart, qend
    identical = np.where(query_array[0] == 1)[0]
    positives = np.where(query_array[1] == 1)[0]

    if identical.shape[0] > q_aligned:
        raise ValueError(
            "Number of identical positions ({}) greater than number of aligned positions ({})!\n{}\n{}".format(
            identical.shape[0], q_aligned, q_intervals, q_merged_intervals))

    if positives.shape[0] > q_aligned:
        raise ValueError("Number of identical positions ({}) greater than number of aligned positions ({})!".format(
            positives.shape[0], q_aligned))

    t_merged_intervals, t_aligned = merge(t_intervals)
    hit_dict["target_aligned_length"] = min(t_aligned, hit.seq_len)
    hit_dict["target_start"] = t_merged_intervals[0][0]
    hit_dict["target_end"] = t_merged_intervals[-1][1] + 1
    hit_dict["global_identity"] = identical.shape[0] * 100 / q_aligned
    hit_dict["global_positives"] = positives.shape[0] * 100 / q_aligned

    return hit_dict, hsp_dict_list
# pylint: enable=too-many-locals
