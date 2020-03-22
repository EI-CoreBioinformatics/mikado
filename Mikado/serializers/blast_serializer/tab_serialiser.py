from .tabular_utils import parse_tab_blast, get_queries, get_targets, TabParserWrapper
import multiprocessing as mp
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
                                   matrix_name=matrix_name,
                                   qmult=qmult, tmult=tmult)
        for fname in self.xml:
            parser(bname=fname)
            self.logger.debug("Finished %s", fname)
    else:
        self.logger.info("Creating a pool with %d workers for analysing BLAST results",
                         self.procs)
        # lock = mp.RLock()
        # queue = mp.JoinableQueue(-1)
        # conf = dict()
        pool = mp.Pool(self.procs)
        queries = get_queries(self.engine)
        targets = get_targets(self.engine)
        parser = functools.partial(parse_tab_blast,
                                   self=self,
                                   queries=queries,
                                   targets=targets,
                                   pool=pool,
                                   matrix_name=matrix_name,
                                   qmult=qmult, tmult=tmult)
        # conf["serialise"] = dict()
        # conf["serialise"]["max_objects"] = self.json_conf["serialise"]["max_objects"]
        # conf["threads"] = self.json_conf["threads"]
        # conf["db_settings"] = self.json_conf["db_settings"].copy()
        #
        # processes = [TabParserWrapper(
        #     conf=conf, idx=idx,
        #     queue=queue, lock=lock,
        #     logging_queue=self.logging_queue, level=self.logger.level,
        #     matrix_name=matrix_name, qmult=qmult, tmult=tmult) for idx in range(1, self.procs + 1)]
        # [proc.start() for proc in processes]
        for fname in self.xml:
            parser(bname=fname)
        # queue.put("EXIT")
        # [proc.join() for proc in processes]

    self.logger.info("Finished loading blast hits")
