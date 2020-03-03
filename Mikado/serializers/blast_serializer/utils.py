"""
Generic utilities used for BLAST serialising into a DB.
"""

from ...parsers.blast_utils import merge
import numpy as np

valid_letters = 'ACDEFGHIKLMNPQRSTVWYBXZJUO'
letters = np.array(list(valid_letters), dtype=np.str_)


__author__ = 'Luca Venturini'

# valid_matches = set([chr(x) for x in range(65, 91)] + [chr(x) for x in range(97, 123)] +
#                     ["|", "*"])


def prepare_hsp(hsp, counter, qmultiplier=1, tmultiplier=1):

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
    match, identical_positions, positives = _prepare_aln_strings(hsp, qmultiplier=qmultiplier)
    hsp_dict["counter"] = counter + 1
    hsp_dict["query_hsp_start"] = hsp.query_start
    hsp_dict["query_hsp_end"] = hsp.query_end
    hsp_dict["query_frame"] = hsp.query_frame
    hsp_dict["target_hsp_start"] = hsp.hit_start
    hsp_dict["target_hsp_end"] = hsp.hit_end
    hsp_dict["target_frame"] = hsp.hit_frame
    hsp_dict["hsp_identity"] = hsp.ident_num / hsp.aln_span * 100
    hsp_dict["hsp_positives"] = hsp.pos_num / hsp.aln_span * 100
    hsp_dict["match"] = match
    hsp_dict["hsp_length"] = hsp.aln_span
    hsp_dict["hsp_bits"] = hsp.bitscore
    hsp_dict["hsp_evalue"] = hsp.evalue
    return hsp_dict, identical_positions, positives


def _np_grouper(data):
    return np.array(np.split(data, np.where(np.diff(data) != 1)[0] + 1))


def _prepare_aln_strings(hsp, qmultiplier=1):

    """This private method calculates the identical positions, the positives, and a re-factored match line
    starting from the HSP.
    :type hsp: Bio.SearchIO.HSP
    """

    lett_array = np.array([
        list(str(hsp.query.seq)),
        list(hsp.aln_annotation["similarity"]),
        list(str(hsp.hit.seq))], dtype=np.str_)

    match = lett_array[1]
    match[np.where(~((lett_array[1] == "+") | (np.isin(lett_array[1], letters))))] = " "
    match[np.where(
        (np.isin(lett_array[0], letters)) &
        (np.isin(lett_array[2], letters)) &
        (lett_array[0] != lett_array[2]) &
        (lett_array[1] != "+"))] = "X"
    match[np.where((lett_array[0] == "-") & (lett_array[2] == "*"))] = "*"
    match[np.where((lett_array[0] == "-") & ~(lett_array[2] == "*"))] = "-"
    match[np.where((lett_array[2] == "-") & (lett_array[0] == "*"))] = "*"
    match[np.where((lett_array[2] == "-") & ~(lett_array[0] == "*"))] = "_"

    summer = np.array([[_] for _ in range(qmultiplier)])
    v = np.array([[1]] * qmultiplier)
    identical_positions = np.where(np.isin(match, valid_letters)) * v
    identical_positions = set(np.array(
        [identical_positions[_] * 3 + summer for _ in range(identical_positions.shape[0])]).flatten())
    positives = np.where(~np.isin(match, np.array(["*", "-", "_", " "]))) * v
    positives = set(np.array(
        [positives[_] * 3 + summer for _ in range(positives.shape[0])]).flatten())
    str_match = "".join(match)

    return str_match, identical_positions, positives


def prepare_hit(hit, query_id, target_id, **kwargs):
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

    identical_positions, positives = set(), set()

    best_hsp = (float("inf"), float("-inf"))  # E-Value, BitS

    def hsp_sorter(val):
        """
        :param val: Evalue, Bit-Score
        :return:
        """

        evalue, bits = val
        return -evalue, bits

    qmulti = kwargs["query_multiplier"]
    tmulti = kwargs["target_multiplier"]
    qlength = kwargs["query_length"]
    assert isinstance(qmulti, (int, float)), type(qmulti)
    assert isinstance(tmulti, (int, float)), type(tmulti)
    hit_dict.update(kwargs)
    hit_dict["query_id"] = query_id
    hit_dict["target_id"] = target_id

    for counter, hsp in enumerate(hit.hsps):
        hsp_dict, ident, posit = prepare_hsp(hsp, counter, qmultiplier=qmulti, tmultiplier=tmulti)
        identical_positions.update(ident)
        positives.update(posit)
        best_hsp = sorted([best_hsp,
                           (hsp_dict["hsp_evalue"], hsp_dict["hsp_bits"])],
                          key=hsp_sorter, reverse=True)[0]

        # if hsp_dict["hsp_evalue"] < best_hsp[0] and hsp_dict["hsp_bits"] > best_hsp[1]:
        #     best_hsp = (hsp_dict["hsp_evalue"], hsp_dict["hsp_bits"])
        hsp_dict["query_id"] = query_id
        hsp_dict["target_id"] = target_id
        hsp_dict_list.append(hsp_dict)
        q_intervals.append((hsp.query_start, hsp.query_end))
        # t_intervals.append((hsp.sbjct_start, hsp.sbjct_end))
        t_intervals.append((hsp.hit_start, hsp.hit_end))

    q_merged_intervals, q_aligned = merge(q_intervals)
    assert isinstance(q_aligned, np.int), (q_merged_intervals, q_aligned, type(q_aligned))
    hit_dict["query_aligned_length"] = min(qlength, q_aligned)
    qstart, qend = q_merged_intervals[0][0], q_merged_intervals[-1][1]
    assert isinstance(qstart, np.int), (q_merged_intervals, type(qstart))
    assert isinstance(qend, np.int), (q_merged_intervals, type(qend))

    hit_dict["query_start"], hit_dict["query_end"] = qstart, qend

    if len(identical_positions) > q_aligned:
        raise ValueError("Number of identical positions ({}) greater than number of aligned positions ({})!".format(
            len(identical_positions), q_aligned))

    if len(positives) > q_aligned:
        raise ValueError("Number of identical positions ({}) greater than number of aligned positions ({})!".format(
            len(positives), q_aligned))

    t_merged_intervals, t_aligned = merge(t_intervals)
    hit_dict["target_aligned_length"] = min(t_aligned, hit.seq_len)
    hit_dict["target_start"] = t_merged_intervals[0][0]
    hit_dict["target_end"] = t_merged_intervals[-1][1]
    hit_dict["global_identity"] = len(identical_positions) * 100 / q_aligned
    hit_dict["global_positives"] = len(positives) * 100 / q_aligned

    return hit_dict, hsp_dict_list
# pylint: enable=too-many-locals
