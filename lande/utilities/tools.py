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
        and converst all

        (a) numpy arrays into python lists
        (b) numpy strings into python scrings.
        (c) an ordered dict to a dict

        which is conveneint for duming a file to yaml.
    """
    if isinstance(x,list):
        return map(tolist,x)
    elif isinstance(x,dict):
        return dict((tolist(k),tolist(v)) for k,v in x.items())
    elif isinstance(x,np.ndarray) or \
            isinstance(x,np.number):
        return x.tolist()
    elif isinstance(x,np.str):
        return str(x)
    elif isinstance(x,PhaseRange):
        return x.tolist(dense=True)
    elif isinstance(x,OrderedDict):
        return dict(x)
    else:
        return x


# Taken form http://stackoverflow.com/questions/4126348/how-do-i-rewrite-this-function-to-implement-ordereddict
class OrderedDefaultdict(OrderedDict):
    def __init__(self, *args, **kwargs):
        newdefault = None
        newargs = ()
        if len(args):
            newdefault = args[0]
            if not callable(newdefault) and newdefault != None:
                raise TypeError('first argument must be callable or None')
            newargs = args[1:]
        self.default_factory = newdefault
        super(OrderedDefaultdict, self).__init__(*newargs, **kwargs)

    def __missing__ (self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = value = self.default_factory()
        return value

