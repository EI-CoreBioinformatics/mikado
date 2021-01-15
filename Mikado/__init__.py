#!/usr/bin/env python3
# coding: utf_8

"""
Mikado is a Python suite whose purpose is to find and resolve genic loci in a
genomic annotation. This is the library it relies onto.
"""

from Mikado.version import __version__

__title__ = "Mikado"
__author__ = 'Luca Venturini'
__license__ = 'GPL3'
__copyright__ = 'Copyright 2015-2020 Luca Venturini'

__all__ = ["configuration",
           "exceptions",
           "loci",
           "parsers",
           "picking",
           "preparation",
           "scales",
           "serializers",
           "subprograms",
           "utilities",
           "__version__"]


from .utilities.log_utils import create_default_logger
from numpy._pytesttester import PytestTester
test = PytestTester(__name__)
