"""
Generic utilities used for BLAST serialising into a DB.
"""

from ...parsers.blast_utils import merge
import numpy as np
import functools


__author__ = 'Luca Venturini'

valid_matches = set([chr(x) for x in range(65, 91)] + [chr(x) for x in range(97, 123)] +
                    ["|", "*"])


def prepare_hsp(hsp, counter, qmultiplier=1, tmultiplier=1):

    r"""
    Prepare a HSP for loading into the DB.
    The match line will be reworked in the following way:

    - If the position is a match/positive, keep the original value
    - If the position is a gap *for the query*, insert a - (dash)
    - If the position is a gap *for the target*, insert a _ (underscore)
    - If the position is a gap *for both*, insert a \ (backslash)

    :param hsp: An HSP object from Bio.Blast.NCBIXML
    # :type hsp: Bio.Blast.Record.HSP
    :type hsp: Bio.Blast.Record.HSP
    :param counter: a digit that indicates the priority of the HSP in the hit
    :return: hsp_dict, identical_positions, positives
    :rtype: (dict, set, set)
    """

    hsp_dict = dict()
    # We must start from 1, otherwise MySQL crashes as its indices start from 1 not 0
    match, identical_positions, positives = _prepare_aln_strings(hsp, qmultiplier=qmultiplier, tmultiplier=tmultiplier)
    hsp_dict["counter"] = counter + 1
    hsp_dict["query_hsp_start"] = hsp.query_start
    hsp_dict["query_hsp_end"] = hsp.query_end
    hsp_dict["query_frame"] = hsp.frame[0]
    hsp_dict["target_hsp_start"] = hsp.sbjct_start
    hsp_dict["target_hsp_end"] = hsp.sbjct_end
    hsp_dict["target_frame"] = hsp.frame[1]
    hsp_dict["hsp_identity"] = hsp.identities / hsp.align_length * 100
    hsp_dict["hsp_positives"] = hsp.positives / hsp.align_length * 100
    hsp_dict["match"] = match
    hsp_dict["hsp_length"] = hsp.align_length
    hsp_dict["hsp_bits"] = hsp.score
    hsp_dict["hsp_evalue"] = hsp.expect
    return hsp_dict, identical_positions, positives


def _prepare_aln_strings(hsp, qmultiplier=1, tmultiplier=1):

    """This private method calculates the identical positions, the positives, and a re-factored match line
    starting from the HSP."""

    # for query_aa, middle_aa, target_aa in zip(hsp.query, hsp.match, hsp.sbjct):
    query_pos, target_pos = hsp.query_start - 1, hsp.sbjct_start - 1

    def categoriser(middle_aa, query_aa, target_aa, qmultiplier, tmultiplier):
        qpos = 0
        tpos = 0
        identical = set()
        positives = set()
        matched = ""
        if query_aa == target_aa == "-":
            matched = "\\"
        elif query_aa == "-":
            tpos = tmultiplier
            if target_aa == "*":
                matched = "*"
            else:
                matched = "-"
        elif target_aa == "-":
            qpos = qmultiplier
            if query_aa == "*":
                matched = "*"
            else:
                matched = "_"
        elif middle_aa == " ":
            matched = " "
            qpos = qmultiplier
            tpos += tmultiplier
        elif middle_aa == "+" or middle_aa in valid_matches:
            if middle_aa != "+":
                identical = set(range(query_pos, query_pos + qmultiplier))
            positives = set(range(query_pos, query_pos + + qmultiplier))
            qpos = qmultiplier
            tpos += tmultiplier
            matched = middle_aa
        return qpos, tpos, matched, identical, positives

    partial_categorizer = functools.partial(categoriser,
                                            qmultiplier=qmultiplier, tmultiplier=tmultiplier)

    results = [partial_categorizer(middle_aa, query_aa, target_aa) for
               middle_aa, query_aa, target_aa in zip(hsp.match, hsp.query, hsp.sbjct)]

    qposes, tposes, matches, identicals, posis = list(zip(*results))
    query_pos += sum(qposes)
    target_pos += sum(tposes)
    match = "".join(matches)
    identical_positions = set.union(*identicals)
    positives = set.union(*posis)

    assert query_pos <= hsp.query_end and target_pos <= hsp.sbjct_end, ((query_pos, hsp.query_end),
                                                                        (target_pos, hsp.sbjct_end),
                                                                        hsp.match, hsp.query, hsp.sbjct)

    return match, identical_positions, positives


def prepare_hit(hit, query_id, target_id, **kwargs):
    """Prepare the dictionary for fast loading of Hit and Hsp objects.
    global_positives: the similarity rate for the global hit *using the query perspective*
    global_identity: the identity rate for the global hit *using the query perspective*

    :param hit: the hit to parse.
    :type hit: Bio.Blast.Record.Alignment

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
        t_intervals.append((hsp.sbjct_start, hsp.sbjct_end))

    q_merged_intervals, q_aligned = merge(q_intervals)
    assert isinstance(q_aligned, np.int), (q_merged_intervals, q_aligned, type(q_aligned))
    hit_dict["query_aligned_length"] = min(q_aligned, hit.query_length)
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
    hit_dict["target_aligned_length"] = min(t_aligned, hit.length)
    hit_dict["target_start"] = t_merged_intervals[0][0]
    hit_dict["target_end"] = t_merged_intervals[-1][1]
    hit_dict["global_identity"] = len(identical_positions) * 100 / q_aligned
    hit_dict["global_positives"] = len(positives) * 100 / q_aligned

    return hit_dict, hsp_dict_list
# pylint: enable=too-many-locals
