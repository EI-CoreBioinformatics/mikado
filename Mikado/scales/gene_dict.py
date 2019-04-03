from ..loci.reference_gene import Gene
import sqlite3
import json


class GeneDict:

    def __init__(self, cursor: sqlite3.Cursor, logger):

        self.__cursor = cursor
        self.logger = logger

    def __getitem__(self, item):

        res = self.__cursor.execute("SELECT json from genes where gid=?", (item,)).fetchmany()
        if res:
            try:
                return self.__load_gene(res[0][0])
            except IndexError:
                raise IndexError(res)
        else:
            return None

    def __load_gene(self, jdict):
        gene = Gene(None, logger=self.logger)
        gene.load_dict(json.loads(jdict))
        gene.finalize()
        return gene

    def load_all(self):

        for row in self.__cursor.execute("SELECT gid, json from genes"):
            gid, gene = row[0], self.__load_gene(row[1])
            yield (gid, gene)

