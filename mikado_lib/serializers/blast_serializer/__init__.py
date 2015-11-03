"""
This package contains all the modules necessary for BLAST serialisation and analysis.
"""

from mikado_lib.serializers.blast_serializer.utils import prepare_hit, prepare_hsp
from mikado_lib.serializers.blast_serializer.query import Query
from mikado_lib.serializers.blast_serializer.target import Target
from mikado_lib.serializers.blast_serializer.hsp import Hsp
from mikado_lib.serializers.blast_serializer.hit import Hit
from mikado_lib.serializers.blast_serializer.xml_serialiser import XmlSerializer

__author__ = 'Luca Venturini'
