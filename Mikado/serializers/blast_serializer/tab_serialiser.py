from .tabular_utils import parse_tab_blast
import pandas as pd
import multiprocessing as mp
from .utils import load_into_db


def _serialise_tabular(self):
    if isinstance(self.xml, str):
        self.xml = [self.xml]
    else:
        assert isinstance(self.xml, (list, set))

    queries = pd.read_sql_table("query", self.engine, index_col="query_name",
                                columns=["query_name", "query_id", "query_length"])
    queries.columns = ["qid", "qlength"]
    targets = pd.read_sql_table("target", self.engine, index_col="target_name",
                                columns=["target_name", "target_id", "target_length"])
    targets.columns = ["sid", "slength"]

    # cache = {"query": self.queries, "target": self.targets}
    matrix_name = self.json_conf["serialise"]["substitution_matrix"]
    program = self.json_conf["serialise"].get("blast_flavour", "blastx")
    qmult, tmult = self.get_multipliers(None, program)

    if self._xml_debug is False and (self.single_thread is True or self.procs == 1):
        pool = None
    else:
        pool = mp.Pool(self.procs)

    hits, hsps = [], []

    for fname in self.xml:
        self.logger.debug("Analysing %s", fname)
        hits, hsps = parse_tab_blast(fname,
                                       queries,
                                       targets,
                                       hits=hits,
                                       hsps=hsps,
                                       pool=pool,
                                       matrix_name=matrix_name,
                                       qmult=qmult, tmult=tmult,
                                       logger=self.logger)
        hits, hsps = load_into_db(self, hits, hsps, force=False)
        self.logger.debug("Finished %s", fname)
    _, _ = load_into_db(self, hits, hsps, force=True)
    self.logger.info("Finished loading blast hits")
