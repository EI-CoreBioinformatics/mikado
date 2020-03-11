from .tabular_utils import parse_tab_blast
import pandas as pd
import multiprocessing as mp
from .utils import load_into_db


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

    if self._xml_debug is False and (self.single_thread is True or self.procs == 1):
        pool = None
    else:
        self.logger.info("Creating a pool with %d workers for analysing BLAST results",
                         self.procs)
        pool = mp.Pool(self.procs)

    hits, hsps = [], []

    for fname in self.xml:
        self.logger.warning("Analysing %s", fname)
        hits, hsps = parse_tab_blast(fname,
                                     queries,
                                     targets,
                                     hits=hits,
                                     hsps=hsps,
                                     pool=pool,
                                     matrix_name=matrix_name,
                                     qmult=qmult, tmult=tmult,
                                     logger=self.logger)
        from collections import defaultdict
        checker = defaultdict(list)
        bad = set()
        for hit in hits:
            key = (hit["query_id"], hit["target_id"])
            checker[key].append(hit)
            if len(checker[key]) > 1:
                bad.add(key)

        if bad:
            self.logger.error("Duplicated keys:\n{}".format("\n".join([str(_) for _ in bad])))
            raise KeyError

        hits, hsps = load_into_db(self, hits, hsps, force=False)
        self.logger.debug("Finished %s", fname)
    _, _ = load_into_db(self, hits, hsps, force=True)
    self.logger.info("Finished loading blast hits")
