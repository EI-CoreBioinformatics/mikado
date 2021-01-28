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

from . import resultstorer
from . import compare
from .contrast import compare as c_compare
