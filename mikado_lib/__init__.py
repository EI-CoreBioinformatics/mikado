#!/usr/bin/env python3
# coding: utf_8

"""
mikado_lib is a Python suite whose purpose is to find and resolve genic loci in a
genomic annotation. This is the library it relies onto.
"""

__title__ = "mikado_lib"
__version__ = '0.9.6'
__author__ = 'Luca Venturini'
__license__ = 'GPL3'
__copyright__ = 'Copyright 2015 Luca Venturini'

import mikado_lib.exceptions
import mikado_lib.utilities
import mikado_lib.parsers
import mikado_lib.serializers
import mikado_lib.loci_objects
import mikado_lib.configuration
import mikado_lib.scales
