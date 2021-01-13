from . import Parser
import pysam


class BamParser(Parser):

    __annot_type__ = "bam"

    def __init__(self, handle):
        self.__closed = False
        self._handle = pysam.AlignmentFile(handle, "rb", check_sq=False)

    def __iter__(self):
        return self

    def __next__(self):
        return next(self._handle)

    @property
    def name(self):
        return self._handle.filename.decode()
