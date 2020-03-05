"""
This package contains all the modules necessary for BLAST serialisation and analysis.
"""

from .xml_utils import prepare_hit, prepare_hsp
from ...exceptions import InvalidHit
from .query import Query
from .target import Target
from .hsp import Hsp
from .hit import Hit
from .blast_serialiser import BlastSerializer

__author__ = 'Luca Venturini'
