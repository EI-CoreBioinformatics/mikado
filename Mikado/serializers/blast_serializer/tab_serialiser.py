from .tabular_utils import parse_tab_blast
import pandas as pd
import multiprocessing as mp
from .utils import load_into_db
import functools


def _serialise_tabular(self):
    if isinstance(self.xml, str):
        self.xml = [self.xml]
    else:
        assert isinstance(self.xml, (list, set))

    queries = pd.read_sql_table("query", self.engine, index_col="query_name")
    queries.columns = ["qid", "qlength"]
    queries["qid"] = queries["qid"].astype(int)
    assert queries.qid.drop_duplicates().shape[0] == queries.shape[0]

    targets = pd.read_sql_table("target", self.engine, index_col="target_name")
    targets.columns = ["sid", "slength"]
    targets["sid"] = targets["sid"].astype(int)
    assert targets.sid.drop_duplicates().shape[0] == targets.shape[0]

    if targets[targets.slength.isna()].shape[0] > 0:
        raise KeyError("Unbound targets!")

    # cache = {"query": self.queries, "target": self.targets}
    matrix_name = self.json_conf["serialise"]["substitution_matrix"]
    program = self.json_conf["serialise"].get("blast_flavour", "blastx")
    qmult, tmult = self.get_multipliers(None, program)
    hits, hsps = [], []

    if self._xml_debug is False and (self.single_thread is True or self.procs == 1):
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
                                   queries=queries,
                                   targets=targets,
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
