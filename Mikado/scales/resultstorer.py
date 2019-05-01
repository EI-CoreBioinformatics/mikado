# coding: utf-8

"""
This class defines the results of the Assigner.compare method.
"""


class ResultStorer:
    """This class stores the results in pre-defined slots, to reduce memory usage."""

    __slots__ = ["ref_id", "ref_gene", "ccode",
                 "tid", "gid",
                 "tid_num_exons", "ref_num_exons",
                 "n_prec", "n_recall", "n_f1",
                 "j_prec", "j_recall", "j_f1",
                 "e_prec", "e_recall", "e_f1",
                 "distance",
                 "location"]

    @staticmethod
    def to_float_tuple(tup):
        if isinstance(tup, (tuple, list)):
            return tuple([float(_) for _ in tup])
        elif isinstance(tup, (str, float, int, bytes)):
            return tuple([float(tup)])
        else:
            raise TypeError(type(tup))

    @staticmethod
    def to_int_tuple(tup):
        if isinstance(tup, (tuple, list)):
            return tuple([int(_) for _ in tup])
        elif isinstance(tup, (str, float, int, bytes)):
            return tuple([int(tup)])
        else:
            raise TypeError(type(tup))

    @staticmethod
    def to_tuple(tup):
        if isinstance(tup, (tuple, list)):
            return tuple([_ for _ in tup])
        elif isinstance(tup, (str, float, int)):
            return tuple([str(tup)])
        elif isinstance(tup, bytes):
            return tuple([tup.decode()])
        else:
            raise TypeError(type(tup))

    @property
    def types(self):
        return [self.to_tuple, self.to_tuple, self.to_tuple,
                str, str, int, int,
                self.to_float_tuple, self.to_float_tuple, self.to_float_tuple,
                self.to_float_tuple, self.to_float_tuple, self.to_float_tuple,
                self.to_float_tuple, self.to_float_tuple, self.to_float_tuple,
                self.to_int_tuple, str]


    def __init__(self, *args, state=None):

        """
        :param args: a list/tuple
        :type args: list | tuple

        """

        if state is not None and len(args) == 0:
            self._load_dict(state)
            return
        
        if len(args) != len(self.__slots__):
            err_msg = "Result_storer expected {0} but only received {1}".format(
                len(self.__slots__), len(args))
            raise ValueError(err_msg)

        self.ref_id, self.ref_gene, self.ccode, self.tid, self.gid, \
            self.tid_num_exons, self.ref_num_exons, \
            self.n_prec, self.n_recall, self.n_f1,\
            self.j_prec, self.j_recall, self.j_f1, \
            self.e_prec, self.e_recall, self.e_f1, \
            self.distance, self.location = args

        for index, key in enumerate(self.__slots__):
            setattr(self, key, self.types[index](getattr(self, self.__slots__[index])))
            #
            # if index < 3:
            #     if isinstance(getattr(self, self.__slots__[index]), str):
            #         setattr(self, key, tuple([getattr(self, self.__slots__[index])]))
            # elif 6 <= index < len(self.__slots__):
            #     if isinstance(getattr(self, self.__slots__[index]), (float, int)):
            #         setattr(self, key, tuple([]))
            #     elif isinstance(getattr(self, self.__slots__[index]), str):
            #         __val = tuple(getattr(self, self.__slots__[index]).split(","))
            #         if str.isdecimal(__val[0]) or str.isdigit(__val[0]):
            #             __val = tuple([float(_) for _ in __val])
            #         setattr(self, key, __val)
            #     elif not isinstance(getattr(self, self.__slots__[index]), tuple):
            #         setattr(self, self.__slots__[index],
            #                 tuple([getattr(self, self.__slots__[index])]))

    def _load_dict(self, state):

        for index, key in enumerate(self.__slots__):
            setattr(self, key, self.types[index](state[key]))
        return
                    
    def _asdict(self):

        """
        :return: a dictionary containing the items of the class
        :rtype : dict
        """
        result_dict = dict().fromkeys(self.__slots__)
        for attr in self.__slots__:
            result_dict[attr] = getattr(self, attr)
            if isinstance(result_dict[attr], (tuple, list)):
                result_dict[attr] = ",".join(["{}".format(_) for _ in result_dict[attr]])
        return result_dict

    def as_dict(self):
        """
        Wrapper for the protected method _asdict
        :return: dictionary
        """

        return self._asdict()

    def __str__(self):

        result_dict = self._asdict()
        line = []
        for key in self.__slots__:
            line.append(str(result_dict[key]))
        return "\t".join(line)

    def __repr__(self):

        represent = "result( "
        for key in self.__slots__:
            represent += "{0}={1}, ".format(key, getattr(self, key))
        represent += ")"
        return represent

    def __getitem__(self, item):
        if item in self.__slots__:
            return getattr(self, item)
        else:
            raise KeyError(item)
