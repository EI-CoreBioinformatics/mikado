class result_storer:
    
    '''This class is used by shangaicompare to store the results of a comparison.'''
    
    __slots__ = [ "RefId", "RefGene", "ccode", "TID", "GID", "n_prec", "n_recall", "n_f1", "j_prec", "j_recall", "j_f1", "distance"]
    
    def __init__(self, *args):
        
        if len(args)!=len(self.__slots__):
            raise ValueError("Result_storer expected {0} but only received {1}".format(len(self.__slots__), len(args)))
        
        for index,key in enumerate(self.__slots__):
            if index<3:
                if type(args[index]) is str:
                    setattr(self, key, tuple([args[index]]) )
                else:
                    setattr(self, key, args[index])
            elif index in (3,4, len(self.__slots__)):
                setattr(self, key, args[index])
            else:
                if type(args[index]) in (float, int):
                    setattr(self, key, tuple([args[index]]))
                else:
                    setattr(self, key, tuple(args[index]))
    
    def _asdict(self):
        
        d=dict().fromkeys(self.__slots__)
        
        for attr in self.__slots__[:3]:
            try:
                d[attr]=",".join(list(getattr(self, attr)))
            except TypeError as exc:
                raise TypeError("{0}; {1}".format(exc, getattr(self, attr)))
        for attr in self.__slots__[3:5]:
            d[attr]=getattr(self, attr)
        for attr in self.__slots__[5:-1]:
            d[attr]=",".join("{0:,.2f}".format(x) for x in getattr(self,attr))
        d["distance"]=self.distance[0] #Last attribute
        return d
    
    def __str__(self):
        
        r = self._asdict()
        line=[]
        for key in self.__slots__:
            line.append(r[key])
        return "\t".join(line)
    
    def __repr__(self):
        
        t="result( "
        for key in self.__slots__:
            t+= "{0}={1} ".format(key,  getattr(self, key))
        t+=")"
        return t
    
 
    
    