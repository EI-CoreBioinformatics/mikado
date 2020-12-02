from .tabular_utils import parse_tab_blast, get_queries, get_targets
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
                                   procs=1,
                                   matrix_name=matrix_name,
                                   qmult=qmult, tmult=tmult)
        for fname in self.xml:
            parser(bname=fname)
            self.logger.debug("Finished %s", fname)
    else:
        self.logger.info("Creating a pool with %d workers for analysing BLAST results",
                         self.procs)
        queries = get_queries(self.engine)
        targets = get_targets(self.engine)
        parser = functools.partial(parse_tab_blast,
                                   self=self,
                                   queries=queries,
                                   targets=targets,
                                   procs=self.procs,
                                   matrix_name=matrix_name,
                                   qmult=qmult, tmult=tmult)
        for fname in self.xml:
            parser(bname=fname)

    self.logger.info("Finished loading blast hits")
