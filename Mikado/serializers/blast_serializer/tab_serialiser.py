from .tabular_utils import parse_tab_blast, get_queries, get_targets
import pandas as pd
import multiprocessing as mp
from .utils import load_into_db
import functools


def _serialise_tabular(self):
    if isinstance(self.xml, str):
        self.xml = [self.xml]
    else:
        assert isinstance(self.xml, (list, set))

    matrix_name = self.json_conf["serialise"]["substitution_matrix"]
    program = self.json_conf["serialise"].get("blast_flavour", "blastx")
    qmult, tmult = self.get_multipliers(None, program)

    if self._xml_debug is False and (self.single_thread is True or self.procs == 1):
        queries = get_queries(self.engine)
        targets = get_targets(self.engine)
        parser = functools.partial(parse_tab_blast,
                                   self=self,
                                   queries=queries,
                                   targets=targets,
                                   conf=None,
                                   logging_queue=None,
                                   logger=self.logger,
                                   matrix_name=matrix_name,
                                   qmult=qmult, tmult=tmult)
        for fname in self.xml:
            parser(bname=fname, identifier=None)
            # _, _ = load_into_db(self, hits, hsps, force=True)
            self.logger.debug("Finished %s", fname)
    else:
        self.logger.info("Creating a pool with %d workers for analysing BLAST results",
                         self.procs)
        lock = mp.RLock()
        parser = functools.partial(parse_tab_blast,
                                   self=None,
                                   queries=None,
                                   targets=None,
                                   conf=self.json_conf,
                                   logger=None,
                                   level=self.logger.level,
                                   logging_queue=self.logging_queue,
                                   matrix_name=matrix_name,
                                   qmult=qmult, tmult=tmult)
        pool = mp.Pool(self.procs)
        for idx, fname in enumerate(self.xml, 1):
            parser(bname=fname, identifier=idx, lock=lock)
            # pool.apply_async(parser, args=(fname,), kwds={"identifier": idx})
        pool.close()
        pool.join()

    self.logger.info("Finished loading blast hits")
