# coding: utf-8

"""
This module defines the classes needed for the "compare" script, namely:

- result_storer:    a glorified struct to
                    hold the comparison results of a prediction vs. the reference
- reference_gene:   data structure which holds the transcripts of a gene.
                    No check is performed - transcripts are grouped according
                    to their parent information.
- assigner:         Main workhorse of the sublibrary. This class
                    assigns each transcript to its best match(es) in the reference annotation.
- accountant:       This class calculates the final summary statistics for the comparison.
"""

from ..utilities import f1

from . import resultstorer

# def calc_f1(recall, precision):
#     """
#     Static method to calculate the F1 statistic given precision
#     and recall (order is unimportant). Definition:
#     F1 = (2 * precision * recall) / (precision + recall)
#
#     :type recall: float
#     :type precision: float
#     :rtype: float
#     """
#     if max(precision, recall) == 0:
#         return 0
#     else:
#         product = 2 * precision * recall
#         summa = precision + recall
#         return product / summa


# noinspection PyPep8
# from . import accountant
# noinspection PyPep8
# from .assignment import assigner
from . import compare
# from .assignment.assigner import Assigner
# from .accountant import Accountant
