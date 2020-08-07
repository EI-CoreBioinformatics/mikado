class Namespace:

    __name__ = "Namespace"

    def __init__(self, default=0, **kwargs):
        self.__default = default
        self.__values = dict()
        self.update(kwargs)

    @property
    def default(self):
        return self.__default

    def __getitem__(self, item):
        self.__dict__.setdefault(item, self.default)
        self.__values[item] = self.__dict__[item]
        return self.__values[item]

    def __getstate__(self):

        default = self.default
        del self.__default
        state = dict()
        state.update(self.__dict__)
        state["default"] = default
        self.__default = default
        return state

    def __setstate__(self, state):
        self.__default = state["default"]
        del state["default"]
        self.__dict__.update(state)

    def __getattr__(self, item):
        self.__dict__.setdefault(item, self.default)
        self.__values[item] = self.__dict__[item]
        return self.__values[item]

    def update(self, dictionary):
        self.__dict__.update(dictionary)
        self.__values.update(dictionary)

    def get(self, item):
        return self.__getitem__(item)

    def __iter__(self):
        return iter(_ for _ in self.__values)

    def __deepcopy__(self, memodict={}):
        new = Namespace(default=self.default)
        new.update(self.__values)
        return new

    def keys(self):
        return self.__values.keys()

    def items(self):
        return self.__values.items()


