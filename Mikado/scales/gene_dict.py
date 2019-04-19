from ..loci.reference_gene import Gene
import sqlite3
import json
from time import sleep
import msgpack


class GeneDict:

    def __init__(self, dbname: str, logger=None):

        self.__dbname = dbname
        self.__db = sqlite3.connect("file:{}?mode=ro".format(self.__dbname), uri=True)
        self.__cursor = self.__db.cursor()
        self.logger = logger

    def __getitem__(self, item):

        failed = 0
        while True:
            try:
                res = self.__cursor.execute("SELECT json from genes where gid=?", (item,)).fetchmany()
                break
            except sqlite3.OperationalError as exc:
                failed += 1
                if failed > 10:
                    raise sqlite3.OperationalError(exc)
                sleep(3)

        if res:
            try:
                return self.__load_gene(res[0][0])
            except IndexError:
                raise IndexError(res)
        else:
            return None

    def __load_gene(self, jdict):
        gene = Gene(None, logger=self.logger)
        gene.load_dict(msgpack.loads(jdict, raw=False))
        gene.finalize()
        return gene

    def load_all(self):

        for row in self.__cursor.execute("SELECT gid, json from genes"):
            gid, gene = row[0], self.__load_gene(row[1])
            yield (gid, gene)

    def __iter__(self):
        self.__db.row_factory = lambda cursor, row: row[0]
        cursor = self.__db.cursor()
        for gid in cursor.execute("SELECT gid from genes"):
            yield gid
        self.__db.row_factory = None

    def items(self):

        for gid in self:
            yield (gid, self[gid])
