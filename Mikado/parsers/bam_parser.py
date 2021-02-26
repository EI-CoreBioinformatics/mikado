from .parser import Parser
import pysam

from ..exceptions import InvalidParsingFormat


class BamParser(Parser):

    __annot_type__ = "bam"

    def __init__(self, handle):
        self.__closed = False
        try:
            self._handle = pysam.AlignmentFile(handle, "rb", check_sq=False)
        except KeyboardInterrupt:
            raise
        except (ValueError, TypeError):
            raise InvalidParsingFormat("This is not a valid BAM file: {}.".format(handle))

    def __iter__(self):
        return self

    def __next__(self):
        return next(self._handle)

    @property
    def name(self):
        return self._handle.filename.decode()
