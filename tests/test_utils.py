import threading
import multiprocessing as mp


class ProcRunner(threading.Thread):

    def __init__(self, function: [mp.Process], *args, **kwargs):

        self.__function = function
        self.args = args
        self.kwargs = kwargs
        self._func = self.__function(*self.args, **self.kwargs)
        super().__init__()

    def run(self):
        self._func.run()

    def __enter__(self):
        return self

    def __exit__(self):
        self.stop()

    @property
    def func(self):
        return self._func

    def join(self, *args, **kwargs):
        if self.func._popen is not None:
            self.func.join()
            self.func.terminate()
        super().join(timeout=0.1)
        if self._tstate_lock is not None:
            assert hasattr(self._tstate_lock, "release")
            if self._tstate_lock.locked():
                self._tstate_lock.release()
        self._stop()

    def stop(self):
        self.join()
        self._stop()