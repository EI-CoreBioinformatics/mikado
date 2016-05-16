#!/usr/bin/env python3
# coding: utf_8

"""
Mikado is a Python suite whose purpose is to find and resolve genic loci in a
genomic annotation. This is the library it relies onto.
"""

__title__ = "Mikado"
__author__ = 'Luca Venturini'
__license__ = 'GPL3'
__copyright__ = 'Copyright 2015-2016 Luca Venturini'
__version__ = "0.22.1"

__all__ = ["configuration",
           "exceptions",
           "loci",
           "parsers",
           "picking",
           "preparation",
           "scales",
           "serializers",
           "subprograms",
           "utilities"]

from . import configuration
from . import exceptions
from . import loci
from . import parsers
from . import picking
from . import preparation
from . import scales
from . import serializers
from . import subprograms
from . import utilities
