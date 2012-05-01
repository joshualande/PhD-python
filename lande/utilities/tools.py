""" This file contains various function which I have found useful, but
    can't think of a better place to put. 
    
    Author: Joshua Lande <joshualande@gmail.com> """
import yaml
from collections import OrderedDict
import numpy as np
from uw.pulsar.phase_range import PhaseRange


def merge_dict(first,second):
    """ recursivly merge two dictionaries. If two conflicting non-dictionary 
        items are present, read in from first by default. """
    ret={}
    for k in set(first.keys()+second.keys()):
        if not first.has_key(k):
            ret[k]=second[k]
        elif not second.has_key(k):
            ret[k]=first[k]
        elif type(first[k])==dict and type(second[k]):
            ret[k]=merge_dict(first[k],second[k])
        else:
            ret[k]=first[k]
    return ret


def tolist(x):
    """ convenience function that takes in a 
        nested structure of lists and dictionaries
        and converts everything to its base objects.
        This is useful for dupming a file to yaml.

        (a) numpy arrays into python lists

            >>> type(tolist(np.asarray(123))) == int
            True
            >>> tolist(np.asarray([1,2,3])) == [1,2,3]
            True

        (b) numpy strings into python strings.

            >>> tolist([np.asarray('cat')])==['cat']
            True

        (c) an ordered dict to a dict

            >>> ordered=OrderedDict(a=1, b=2)
            >>> type(tolist(ordered)) == dict
            True

        (d) converts unicode to regular strings

            >>> type(u'a') == str
            False
            >>> type(tolist(u'a')) == str
            True

        (e) converts numbers & bools in strings to real represntation,
            (i.e. '123' -> 123)

            >>> type(tolist(np.asarray('123'))) == int
            True
            >>> type(tolist('123')) == int
            True
            >>> tolist('False') == False
            True
    """
    if isinstance(x,list):
        return map(tolist,x)
    elif isinstance(x,dict):
        return dict((tolist(k),tolist(v)) for k,v in x.items())
    elif isinstance(x,np.ndarray) or \
            isinstance(x,np.number):
        # note, call tolist again to convert strings of numbers to numbers
        return tolist(x.tolist())
    elif isinstance(x,PhaseRange):
        return x.tolist(dense=True)
    elif isinstance(x,OrderedDict):
        return dict(x)
    elif isinstance(x,basestring) or isinstance(x,np.str):
        x=str(x) # convert unicode & numpy strings 
        try:
            return int(x)
        except:
            try:
                return float(x)
            except:
                if x == 'True': return True
                elif x == 'False': return False
                else: return x
    else:
        return x


class OrderedDefaultDict(OrderedDict):
    """ Ordered default dictionary:

        Code orignally taken from:
            http://stackoverflow.com/questions/4126348/how-do-i-rewrite-this-function-to-implement-ordereddict
        But has been retaken from:
            http://www.techques.com/question/1-6190331/Can-I-do-an-ordered,-default-dict-in-Python
        to have a better __copy__ support.

        Some simple tests

            >>> x = OrderedDefaultDict(list)
            >>> x[1].append(2)
            >>> x['a'].append('b')
            >>> print x.keys()
            [1, 'a']

        For a while, deepcopying these objects crashes. This has been fixed:
            
            >>> import copy
            >>> y=copy.deepcopy(x)
            >>> print y.keys()
            [1, 'a']
            
    """
    def __init__(self, default_factory=None, *a, **kw):
        if (default_factory is not None and
            not hasattr(default_factory, '__call__')):
            raise TypeError('first argument must be callable')
        OrderedDict.__init__(self, *a, **kw)
        self.default_factory = default_factory

    def __getitem__(self, key):
        try:
            return OrderedDict.__getitem__(self, key)
        except KeyError:
            return self.__missing__(key)

    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = value = self.default_factory()
        return value

    def __reduce__(self):
        if self.default_factory is None:
            args = tuple()
        else:
            args = self.default_factory,
        return type(self), args, None, None, self.items()

    def copy(self):
        return self.__copy__()

    def __copy__(self):
        return type(self)(self.default_factory, self)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
