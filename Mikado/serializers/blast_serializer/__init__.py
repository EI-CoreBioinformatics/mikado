"""
This package contains all the modules necessary for BLAST serialisation and analysis.
"""

from .utils import prepare_hit, prepare_hsp
from ...exceptions import InvalidHit
from .query import Query
from .target import Target
from .hsp import Hsp
from .hit import Hit
from .xml_serialiser import XmlSerializer

__author__ = 'Luca Venturini'
