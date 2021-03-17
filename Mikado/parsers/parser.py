import abc
from ..utilities.file_type import filetype
import io
import gzip
import bz2
from functools import partial


class HeaderError(Exception):
    """
    Mock exception which is raised when a header/comment line (e.g. starting with "#") is found.
    """
    pass


class Parser(metaclass=abc.ABCMeta):
    """Generic parser iterator. Base parser class."""

    def __init__(self, handle):
        self.__closed = False
        self._handle = self.__get_handle(handle)
        self.closed = False

    def __iter__(self):
        return self

    def __get_handle(self, handle, position=None):
        if not isinstance(handle, io.IOBase):
            if handle.endswith(".gz") or filetype(handle) == b"application/gzip":
                opener = gzip.open
            elif handle.endswith(".bz2") or filetype(handle) == b"application/x-bzip2":
                opener = bz2.open
            else:
                opener = partial(open, **{"buffering": 1})
            try:
                handle = opener(handle, "rt")
            except FileNotFoundError:
                raise FileNotFoundError("File not found: {0}".format(handle))
        if position is not None:
            handle.seek(position)
        return handle

    def __next__(self):
        return next(self._handle)

    def __enter__(self):
        if self.closed is True:
            raise ValueError('I/O operation on closed file.')
        return self

    def __exit__(self, *args):
        _ = args
        self._handle.close()
        self.closed = True

    def close(self):
        """
        Alias for __exit__
        """
        self.__exit__()

    @property
    def name(self):
        """
        Return the filename.
        """
        if hasattr(self._handle, "name"):
            return self._handle.name
        return None

    @property
    def closed(self):
        """
        Boolean flag. If True, the file has been closed already.
        """
        return self.__closed

    @closed.setter
    def closed(self, *args):
        """
        :param args: boolean flag

        This sets the closed flag of the file.

        """
        if not isinstance(args[0], bool):
            raise TypeError("Invalid value: {0}".format(args[0]))

        self.__closed = args[0]

    def __getstate__(self):
        try:
            position = self._handle.tell()
        except:
            position = None
        state = dict()
        state.update(self.__dict__)
        state["position"] = position
        if hasattr(self._handle, "filename"):
            _handle = self._handle.filename
            if isinstance(_handle, bytes):
                _handle = _handle.decode()
            state["_handle"] = _handle
        elif hasattr(self._handle, "name"):
            _handle = self._handle.name
            if isinstance(_handle, bytes):
                _handle = _handle.decode()
            state["_handle"] = _handle
        else:
            raise TypeError("Unknown handle: {}".format(self._handle))
        state.pop("logger", None)
        return state

    def __setstate__(self, state):
        position = state.get("position")
        del state["position"]
        self.__dict__.update(state)
        self._handle = self.__get_handle(state["_handle"], position=position)

