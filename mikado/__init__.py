#!/usr/bin/env python3
# coding: utf_8

"""
mikado is a Python suite whose purpose is to find and resolve genic loci in a
genomic annotation. This is the library it relies onto.
"""

__title__ = "mikado"
__version__ = '0.9.6'
__author__ = 'Luca Venturini'
__license__ = 'GPL3'
__copyright__ = 'Copyright 2015 Luca Venturini'

import mikado.exceptions
import mikado.utilities
import mikado.parsers
import mikado.serializers
import mikado.loci_objects
import mikado.configuration
import mikado.scales
