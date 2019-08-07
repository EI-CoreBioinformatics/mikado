# coding:utf-8

"""
Custom exceptions for Mikado.
"""


class NotInLocusError(AssertionError):
    """
    Error to be raised when a method tries to add a transcript to a Locus it does not belong to.
    """
    pass


class NoJsonConfigError(ValueError):
    """
    Exception to be raised if no/an invalid configuration dictionary has been provided.
    """
    pass


class InvalidLocusError(ValueError):
    """
    Exception to be raised when something has made a Locus object invalid.
    """
    pass


class ModificationError(RuntimeError):
    """This exception is raised when something tries to modify a finalized object."""
    pass


class UnrecognizedOperator(ValueError):
    """
    Exception to be raised when the configuration file contains an unsupported operator
    for a scoring/filtering function.
    """
    pass


class UnrecognizedRescaler(ValueError):
    """
    Exception to be raised when the configuration file contains an unsupported rescaling
    function for a scoring operation.
    """

    pass


class InvalidJson(KeyError):
    """
    Exception to be raised when the JSON/YAML is invalid.
    """
    pass


class InvalidTranscript(ValueError):
    """
    Exception to be raised when a transcript contains corrupted data
    (e.g. overlapping or missing exons).
    """

    pass


class InvalidCDS(InvalidTranscript):
    """
    Exception to be raised when a transcript contains an invalid CDS
    (e.g. with a UTR in the middle).
    """

    pass


class IncorrectStrandError(InvalidTranscript):
    """
    Exception to be raised when a transcript contains an invalid intron
    (e.g. an AG-GT intron assigned to the minus strand).
    """

    pass


class UnsortedInput(ValueError):
    """
    Exception to be raised when the input for pick is not properly sorted.
    """

    pass


class InvalidAssembly(ValueError):
    """
    Exception to be raised when the input for prepare is not properly formatted.
    """

    pass


class RedundantNames(KeyError):
    pass


class CorruptIndex(ValueError):
    """
    Exception to be raised when the index for Mikado compare is corrupt.
    """

class InvalidHit(ValueError):
    """
    Exception to be raised when a Hit has discrepant declared and real best e-value.
    """


class InvalidSerialization(KeyError):
    """
    Exception to be raised when trying to include in the Mikado database incongruent data.
    """
