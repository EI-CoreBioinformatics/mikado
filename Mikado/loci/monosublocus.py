#!/usr/bin/env python3
# coding: utf-8


"""
This module defines monosubloci, an intermediate step between the definition of subloci
and the final definition of loci.
A Monosublocus is characterized by its containing one and only one transcript.
"""

from .abstractlocus import Abstractlocus


# pylint: disable=too-many-instance-attributes,too-many-public-methods
class Monosublocus(Abstractlocus):
    """Very basic class which holds a single transcript."""

    __name__ = "Monosublocus"

    # ########## Special methods ############

    def __init__(self,
                 transcript_instance=None,
                 json_conf=None,
                 logger=None,
                 verified_introns=None,
                 **kwargs):

        self.counter = 0  # simple tag to avoid collisions
        Abstractlocus.__init__(self, json_conf=json_conf, logger=logger, verified_introns=verified_introns, **kwargs)
        # this must be defined straight away
        if transcript_instance is not None:
            Abstractlocus.add_transcript_to_locus(self, transcript_instance)
            self.tid = transcript_instance.id
            self.monoexonic = transcript_instance.monoexonic
        self.feature = "Monosublocus"
        self.parent = None
        self.attributes = dict()

    # pylint: disable=arguments-differ
    def __str__(self, print_cds=True, source_in_name=True):

        raise NotImplementedError(
            """This is a container used for computational purposes only,
            it should not be printed out directly!""")

    # pylint: enable=arguments-differ

    # ########## Class instance methods ##############

    def add_transcript_to_locus(self, transcript, check_in_locus=False, **kwargs):
        """For this basic class, this method raises a NotImplementedError -
        as this container should hold only one transcript.

        :param transcript
        :param check_in_locus: flag. Ignored.
        :type check_in_locus: bool
        """

        raise NotImplementedError("In a Monosublocus there should be one and only one transcript!")

    def is_intersecting(self):
        """Not implemented: this function makes no sense for a single-transcript container."""
        raise NotImplementedError(
            "Monosubloci hold a single transcript, so intersections are not calculated.")

    # ######### Properties ############

    @property
    def id(self):
        """
        Override of the Abstractlocus method, to set the name appropriately.
        :rtype : str
        """
        if self.monoexonic is True:
            addendum = "mono"
        else:
            addendum = "multi"
        if self.counter > 0:
            addendum = "{0}.{1}".format(addendum, self.counter)
        return "{0}.{1}".format(super().id, addendum)
