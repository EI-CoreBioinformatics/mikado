# coding: utf-8

"""
This class defines the results of the Assigner.compare method.
"""


class RestultStorer:
    """This class stores the results in pre-defined slots, to reduce memory usage."""

    __slots__ = ["RefId", "RefGene", "ccode",
                 "TID", "GID", "n_prec",
                 "n_recall", "n_f1", "j_prec",
                 "j_recall", "j_f1", "distance"]

    def __init__(self, *args):

        """
        :param args: a list/tuple
        :type args: list | tuple

        """

        if len(args) != len(self.__slots__):
            raise ValueError("Result_storer expected {0} but only received {1}".format(len(self.__slots__), len(args)))

        self.RefId, self.RefGene, self.ccode, self.TID, self.GID, \
            self.n_prec, self.n_recall, self.n_f1, self.j_prec, self.j_recall, \
            self.j_f1, self.distance = args

        for index, key in enumerate(self.__slots__):
            if index < 3:
                if type(getattr(self, self.__slots__[index])) is str:
                    setattr(self, key, tuple([getattr(self, self.__slots__[index])]))
            elif 4 < index < len(self.__slots__):
                if type(getattr(self, self.__slots__[index])) in (float, int):
                    setattr(self, key, tuple( [getattr(self, self.__slots__[index])]))

    def _asdict(self):

        """
        :return: a dictionary containing the items of the class
        :rtype : dict
        """
        d = dict().fromkeys(self.__slots__)

        for attr in self.__slots__[:3]:
            try:
                d[attr] = ",".join(list(getattr(self, attr)))
            except TypeError as exc:
                raise TypeError("{0}; {1}".format(exc, getattr(self, attr)))
        for attr in self.__slots__[3:5]:
            d[attr] = getattr(self, attr)
        for attr in self.__slots__[5:-1]:
            d[attr] = ",".join("{0:,.2f}".format(x) for x in getattr(self, attr))
        d["distance"] = self.distance[0]  # Last attribute
        return d

    def __str__(self):

        r = self._asdict()
        line = []
        for key in self.__slots__:
            line.append(str(r[key]))
        return "\t".join(line)

    def __repr__(self):

        t = "result( "
        for key in self.__slots__:
            t += "{0}={1}, ".format(key, getattr(self, key))
        t += ")"
        return t
