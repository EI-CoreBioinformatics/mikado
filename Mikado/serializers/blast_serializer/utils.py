"""
Generic utilities used for BLAST serialising into a DB.
"""

from ...parsers.blast_utils import merge
from ...exceptions import InvalidHit
import operator

__author__ = 'Luca Venturini'

valid_matches = set([chr(x) for x in range(65, 91)] + [chr(x) for x in range(97, 123)] +
                    ["|", "*"])


def prepare_hsp(hsp, counter):

    r"""
    Prepare a HSP for loading into the DB.
    The match line will be reworked in the following way:

    - If the position is a match/positive, keep the original value
    - If the position is a gap *for the query*, insert a - (dash)
    - If the position is a gap *for the target*, insert a _ (underscore)
    - If the position is a gap *for both*, insert a \ (backslash)

    :param hsp: An HSP object from Bio.Blast.NCBIXML
    # :type hsp: Bio.Blast.Record.HSP
    :type hsp: Bio.SearchIO._model.hsp.HSP
    :param counter: a digit that indicates the priority of the HSP in the hit
    :return: hsp_dict, identical_positions, positives
    :rtype: (dict, set, set)
    """

    identical_positions, positives = set(), set()

    hsp_dict = dict()
    # We must start from 1, otherwise MySQL crashes
    # as its indices start from 1 not 0
    hsp_dict["counter"] = counter + 1
    hsp_dict["query_hsp_start"] = hsp.query_start
    hsp_dict["query_hsp_end"] = hsp.query_end
    # Prepare the list for later calculation
    # q_intervals.append((hsp.query_start, hsp.query_end))

    # hsp_dict["target_hsp_start"] = hsp.sbjct_start
    # hsp_dict["target_hsp_end"] = hsp.sbjct_end
    hsp_dict["target_hsp_start"] = hsp.hit_start
    hsp_dict["target_hsp_end"] = hsp.hit_end

    # Prepare the list for later calculation
    # t_intervals.append((hsp.sbjct_start, hsp.sbjct_end))

    # hsp_dict["hsp_identity"] = hsp.identities / hsp.align_length * 100
    # hsp_dict["hsp_positives"] = hsp.positives / hsp.align_length * 100
    hsp_dict["hsp_identity"] = hsp.ident_num / hsp.aln_span * 100
    hsp_dict["hsp_positives"] = hsp.pos_num / hsp.aln_span * 100

    # Prepare the list for later calculation
    # hit_dict["global_identity"].append(hsp_dict["hsp_identity"])
    match = ""
    query_pos, target_pos = hsp.query_start - 1, hsp.hit_start - 1
    positive_count, iden_count = 0, 0
    # for query_aa, middle_aa, target_aa in zip(hsp.query, hsp.match, hsp.sbjct):
    for middle_aa, query_aa, target_aa in zip(hsp.aln_annotation["similarity"],
                                              *hsp.aln.get_all_seqs()):
        if middle_aa in valid_matches or middle_aa == "+":
            query_pos += 1
            target_pos += 1
            match += middle_aa
            positives.add(query_pos)
            positive_count += 1
            if middle_aa != "+":
                iden_count += 1
                identical_positions.add(query_pos)
        elif query_aa == target_aa == "-":
            match += "\\"
        elif query_aa == "-":
            target_pos += 1
            match += "-"
        elif target_aa == "-":
            query_pos += 1
            match += "_"

    assert query_pos <= hsp.query_end and target_pos <= hsp.hit_end
    # assert positive_count == hsp.positives and iden_count == hsp.identities
    # assert positive_count == hsp.pos_num and iden_count == hsp.ident_num

    hsp_dict["match"] = match

    hsp_dict["hsp_length"] = hsp.aln_span
    hsp_dict["hsp_bits"] = hsp.bitscore
    hsp_dict["hsp_evalue"] = hsp.evalue

    return hsp_dict, identical_positions, positives


# Splitting this function would make it less clear, and the internal variables
# are all necessary for a correct serialisation.
# pylint: disable=too-many-locals

# class InvalidHit(ValueError):
#
#     pass


def prepare_hit(hit, query_id, target_id, **kwargs):
    """Prepare the dictionary for fast loading of Hit and Hsp objects.
    global_positives: the similarity rate for the global hit *using the query perspective*
    global_identity: the identity rate for the global hit *using the query perspective*

    :param hit: the hit to parse.
    :type hit: Bio.SearchIO._model.hit.Hit

    :param query_id: the numeric ID of the query in the database. Necessary for serialisation.
    :type query_id: int

    :param target_id: the numeric ID of the target in the database. Necessary for serialisation.
    :type target_id: int

    :param kwargs: additional properties to give to the hit_dict. Retrieved e.g. from descriptions.
    :type kwargs: dict
    """

    hit_dict = dict()
    hsp_dict_list = []
    # hit_dict["global_identity"] = []
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

    for counter, hsp in enumerate(hit.hsps):
        hsp_dict, ident, posit = prepare_hsp(hsp, counter)
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

    hit_dict.update(kwargs)
    hit_dict["query_id"] = query_id
    hit_dict["target_id"] = target_id

    q_merged_intervals = sorted(merge(q_intervals), key=operator.itemgetter(0, 1))
    q_aligned = sum([tup[1] - tup[0] + 1 for tup in q_merged_intervals])
    hit_dict["query_aligned_length"] = q_aligned
    hit_dict["query_start"] = q_merged_intervals[0][0]
    hit_dict["query_end"] = q_merged_intervals[-1][1]

    t_merged_intervals = sorted(merge(t_intervals), key=operator.itemgetter(0, 1))
    t_aligned = sum([tup[1] - tup[0] + 1 for tup in t_merged_intervals])
    hit_dict["target_aligned_length"] = t_aligned
    hit_dict["target_start"] = t_merged_intervals[0][0]
    hit_dict["target_end"] = t_merged_intervals[-1][1]
    hit_dict["global_identity"] = len(identical_positions) * 100 / q_aligned
    hit_dict["global_positives"] = len(positives) * 100 / q_aligned
    # if hit_dict["evalue"] != best_hsp[0] or hit_dict["bits"] != best_hsp[1]:
    #     raise InvalidHit("Discrepant evalue/bits for hsps and hit for {0} vs. {1}; \
    #     best: {2}, reported {3}".format(
    #         hit.id,
    #         query_id,
    #         best_hsp,
    #         (hit_dict["evalue"], hit_dict["bits"])
    #     ))

    return hit_dict, hsp_dict_list
# pylint: enable=too-many-locals
