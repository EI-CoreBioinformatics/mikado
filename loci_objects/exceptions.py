class NotInLocusError(AssertionError):
    pass

class NoJsonConfigError(ValueError):
    pass

class InvalidLocusError(ValueError):
    pass

class ModificationError(RuntimeError):
    '''This exception is raised when something tries to modify a finalized object.'''
    pass

class UnrecognizedOperator(ValueError):
    pass

class UnrecognizedRescaler(ValueError):
    pass

class InvalidJson(KeyError):
    pass

class InvalidTranscript(ValueError):
    pass

class IncorrectStrandError(Exception):
    pass