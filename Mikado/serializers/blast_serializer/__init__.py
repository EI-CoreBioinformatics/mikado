"""
This package contains all the modules necessary for BLAST serialisation and analysis.
"""

from ...exceptions import InvalidHit
from .query import Query
from .target import Target
from .hsp import Hsp, prepare_hsp
from .hit import Hit, prepare_hit
from .blast_serialiser import BlastSerializer

__author__ = 'Luca Venturini'
