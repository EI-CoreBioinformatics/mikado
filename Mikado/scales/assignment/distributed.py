import multiprocessing as mp
from ...utilities.namespace import Namespace
from .assigner import Assigner
from ...transcripts import Transcript
import os
import msgpack


class Assigners(mp.Process):

    def __init__(self, index, args: Namespace, queue, returnqueue, log_queue, counter, dump_dbname):
        super().__init__()
        self._dump_dbname = dump_dbname
        # self.accountant_instance = Accountant(genes, args, counter=counter)
        if hasattr(args, "fuzzymatch"):
            self.__fuzzymatch = args.fuzzymatch
        else:
            self.__fuzzymatch = 0
        self.__counter = counter
        self._index = index
        self.queue = queue
        self.returnqueue = returnqueue
        self._args = args
        self.log_queue = log_queue

    def run(self):
        self.__connection = open(self._dump_dbname, mode="rb")
        self._args.__dict__["log_queue"] = self.log_queue
        self.assigner_instance = Assigner(self._index, self._args,
                                          printout_tmap=False,
                                          counter=self.__counter,
                                          fuzzymatch=self.__fuzzymatch)
        while True:
            transcr = self.queue.get()
            if transcr == "EXIT":
                self.queue.put_nowait("EXIT")
                self.assigner_instance.dump()
                self.returnqueue.put(self.assigner_instance.db.name)
                break
            else:
                assert os.path.exists(self._dump_dbname), self._dump_dbname
                start, end = transcr
                if end - start == 0:
                    continue
                try:
                    assert os.stat(self._dump_dbname).st_size > 0, self._dump_dbname
                except AssertionError:
                    raise AssertionError((start, end, self._dump_dbname))
                self.__connection.seek(start)
                dumped = self.__connection.read(end - start)
                # dumped = self.__cursor.execute("SELECT json FROM dump WHERE idx=?", (transcr,)).fetchone()
                dumped = msgpack.loads(dumped)
                transcr = Transcript()
                transcr.load_dict(dumped, trust_orf=True, accept_undefined_multi=True)
                self.assigner_instance.get_best(transcr)

        self.__connection.close()


class FinalAssigner(mp.Process):

    def __init__(self, index: str, args: Namespace, queue, log_queue):

        super().__init__()
        self.index = index
        self.queue = queue
        failed = set()
        self.args = args
        self.log_queue = log_queue

    def run(self):
        try:
            self.args.__dict__["log_queue"] = self.log_queue
        except AttributeError:
            raise AttributeError(self.args)
        self.assigner = Assigner(self.index, self.args, printout_tmap=True)
        while True:
            dbname = self.queue.get()
            if dbname == "EXIT":
                break
            self.assigner.load_result(dbname)
        self.assigner.finish()


