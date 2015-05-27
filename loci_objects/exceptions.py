class NotInLocusError(AssertionError):
    pass

class NoJsonConfigError(ValueError):
    pass

class InvalidLocusError(ValueError):
    pass

class ModificationError(RuntimeError):
    '''This exception is raised when something tries to modify a finalized object.'''
    pass