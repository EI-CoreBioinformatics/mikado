import multiprocessing as mp
from ...utilities.namespace import Namespace
from .assigner import Assigner
from ...transcripts import Transcript
from ..resultstorer import ResultStorer
import msgpack


def encode_result(obj):
    if isinstance(obj, ResultStorer):
        return obj.as_dict()
    else:
        return obj


class Assigners(mp.Process):

    def __init__(self, index, args: Namespace, queue, returnqueue, log_queue, counter):
        super().__init__()
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
        self._args.__dict__["log_queue"] = self.log_queue
        self.assigner_instance = Assigner(self._index, self._args,
                                          printout_tmap=False,
                                          counter=self.__counter,
                                          fuzzymatch=self.__fuzzymatch)
        results = []
        while True:
            transcripts = self.queue.get()
            if transcripts == "EXIT":
                self.queue.put("EXIT")
                self.queue.task_done()
                result = msgpack.dumps([res.as_dict() for res in results], strict_types=True)
                self.returnqueue.put(("tmap", result))
                refmap, stats = self.assigner_instance.dump()
                self.returnqueue.put(("refmap", refmap, stats))
                self.returnqueue.put("EXIT")
                break
            else:
                dumped = msgpack.loads(transcripts)
                for dump in dumped:
                    transcr = Transcript()
                    transcr.load_dict(dump, trust_orf=True, accept_undefined_multi=True)
                    result = self.assigner_instance.get_best(transcr)
                    if isinstance(result, ResultStorer):
                        result = [result]
                    results.extend(result)
                    if len(results) >= 1000:
                        result = msgpack.dumps([res.as_dict() for res in results], strict_types=True)
                        self.returnqueue.put(("tmap", result))
                        results = []
                self.queue.task_done()


class FinalAssigner(mp.Process):

    def __init__(self, index: str, args: Namespace, queue, log_queue, nprocs):

        super().__init__()
        self.index = index
        self.queue = queue
        # failed = set()
        self.args = args
        self.log_queue = log_queue
        self.nprocs = nprocs

    def run(self):
        try:
            self.args.__dict__["log_queue"] = self.log_queue
        except AttributeError:
            raise AttributeError(self.args)
        self.assigner = Assigner(self.index, self.args, printout_tmap=True)
        finished_children = 0
        while finished_children < self.nprocs:
            tmap_row = self.queue.get()
            if tmap_row == "EXIT":
                finished_children += 1
                self.queue.task_done()
                continue
            elif tmap_row[0] == "tmap":
                rows = [ResultStorer(state=row) for row in msgpack.loads(tmap_row[1])]
                for row in rows:
                    self.assigner.print_tmap(row)
            elif tmap_row[0] == "refmap":
                self.assigner.load_result(tmap_row[1], tmap_row[2])
            self.queue.task_done()

        self.queue.join()
        self.assigner.finish()
        self.assigner.logger.info("Finished everything, shutting down")
        return
