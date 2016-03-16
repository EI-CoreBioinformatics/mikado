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
__version__ = "0.19.0"

__all__ = ["exceptions",
           "utilities",
           "parsers",
           "serializers",
           "loci_objects",
           "configuration",
           "scales"]

from . import exceptions
from . import utilities
from . import configuration
from . import parsers
from . import serializers
from . import loci
from . import scales
