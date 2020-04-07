"""
Generic utilities used for BLAST serialising into a DB.
"""

import numpy as np
from ...parsers.blast_utils import merge
from .aln_string_parser import prepare_aln_strings


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
        if macro > 0 or minor > 9 or (micro > 30 or micro < 25):
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


def _np_grouper(data):
    return np.array(np.split(data, np.where(np.diff(data) != 1)[0] + 1))
