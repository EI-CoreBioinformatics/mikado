"""
Generic utilities used for BLAST serialising into a DB.
"""

from ...parsers.blast_utils import merge
import numpy as np

# valid_matches = set([chr(x) for x in range(65, 91)] + [chr(x) for x in range(97, 123)] +
#                     ["|", "*"])


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
    match[np.where((lett_array[0] == "-") & (lett_array[2] == "*"))] = "*"
    match[np.where((lett_array[0] == "-") & ~(lett_array[2] == "*"))] = "-"
    match[np.where((lett_array[2] == "-") & (lett_array[0] == "*"))] = "*"
    match[np.where((lett_array[2] == "-") & ~(lett_array[0] == "*"))] = "_"
    match[np.where(
        (lett_array[0] != lett_array[2]) & (lett_array[1] != "+") &
        (match != "*") & (match != "_") & (match != "-")
    )] = "/"

    query_array = np.zeros(hsp.query_end - hsp.query_start)
    qpos = 0
    for idx in range(lett_array.shape[1]):
        pos = lett_array[:, idx]
        if pos[0] == pos[2]:
            query_array[qpos] = 2
            qpos += 1
        elif pos[1] == "+":
            query_array[qpos] = 1
            qpos += 1
        elif pos[0] == "-":
            continue
        elif pos[2] == "-":
            qpos += 1
            continue

    summer = np.array([[_] for _ in range(qmultiplier)])
    _id_catcher = np.where(query_array >= 2)
    assert hsp.ident_num == _id_catcher[0].shape[0]
    identical_positions = ((_id_catcher[0] * qmultiplier) + summer).flatten()
    _pos_catcher = np.where(query_array >= 1)
    assert hsp.pos_num == _pos_catcher[0].shape[0], (hsp.pos_num, _pos_catcher[0].shape[0])
    positives = ((_pos_catcher[0] * qmultiplier) + summer).flatten()
    if hsp.query_frame > 0:
        identical_positions = identical_positions + hsp.query_start
        positives = positives + hsp.query_start
    else:
        identical_positions = hsp.query_end - identical_positions
        positives = hsp.query_end - positives

    assert hsp.query_start <= positives.min() <= positives.max() <= hsp.query_end, (
        hsp.query_frame, hsp.query_start, positives.min(), positives.max(), hsp.query_end
    )
    assert hsp.query_start <= identical_positions.min() <= identical_positions.max() <= hsp.query_end
    identical_positions = set(identical_positions)
    positives = set(positives)
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

    qmulti = kwargs["query_multiplier"]
    tmulti = kwargs["target_multiplier"]
    qlength = kwargs["query_length"]
    assert isinstance(qmulti, (int, float)), type(qmulti)
    assert isinstance(tmulti, (int, float)), type(tmulti)
    hit_dict.update(kwargs)
    hit_dict["query_id"] = query_id
    hit_dict["target_id"] = target_id

    query_array = np.zeros([3, qlength])
    target_array = np.zeros([3, hit.seq_len])

    for counter, hsp in enumerate(hit.hsps):
        hsp_dict, ident, posit = prepare_hsp(hsp, counter, qmultiplier=qmulti, tmultiplier=tmulti)
        identical_positions.update(ident)
        positives.update(posit)
        hsp_dict["query_id"] = query_id
        hsp_dict["target_id"] = target_id
        hsp_dict_list.append(hsp_dict)
        q_intervals.append((hsp.query_start, hsp.query_end))
        t_intervals.append((hsp.hit_start, hsp.hit_end))

    q_merged_intervals, q_aligned = merge(q_intervals)
    hit_dict["query_aligned_length"] = min(qlength, q_aligned)
    qstart, qend = q_merged_intervals[0][0], q_merged_intervals[-1][1]
    hit_dict["query_start"], hit_dict["query_end"] = qstart, qend

    if len(identical_positions) > q_aligned:
        raise ValueError(
            "Number of identical positions ({}) greater than number of aligned positions ({})!\n{}\n{}".format(
            len(identical_positions), q_aligned, q_intervals, q_merged_intervals))

    if len(positives) > q_aligned:
        raise ValueError("Number of identical positions ({}) greater than number of aligned positions ({})!".format(
            len(positives), q_aligned))

    t_merged_intervals, t_aligned = merge(t_intervals)
    hit_dict["target_aligned_length"] = min(t_aligned, hit.seq_len)
    hit_dict["target_start"] = t_merged_intervals[0][0]
    hit_dict["target_end"] = t_merged_intervals[-1][1]
    hit_dict["global_identity"] = len(identical_positions) * 100 / qlength  # q_aligned
    hit_dict["global_positives"] = len(positives) * 100 / qlength   # q_aligned

    return hit_dict, hsp_dict_list
# pylint: enable=too-many-locals
