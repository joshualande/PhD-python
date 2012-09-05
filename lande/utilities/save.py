import os
from os.path import expandvars
import argparse

import yaml
from numpy.core.records import fromarrays

from uw.utilities.makerec import makefits, fitsrec

from . tools import tolist


def loaddict(filename):
    filename = expandvars(filename)

    extension = os.path.splitext(filename)[-1]

    if extension == '.yaml':
        return yaml.load(open(filename,'r'))
    elif extension == '.hdf5':
        import h5py
        return h5py.File(filename,'r')
    elif extension == '.fits':
        return fitsrec(filename)
    elif extension == '.xml':
        from pyxml2obj import XMLin
        input=XMLin(open(filename).read())
        # pass through tolist for automatic type conversion
        return tolist(input)
    else:
        raise Exception("Unrecognized extension %s" % extension)


def savedict(filename, results, yaml_kwargs=dict()):
    """ Save a dictionary to a file. Choose file format based
        upon extension to filename. """
    is_results = lambda x: isinstance(x,dict) or isinstance(x,list)

    if isinstance(filename,str) and is_results(results):
        pass
    elif is_results(filename) and isinstance(results, str):
        filename, results = results, filename
    else:
        raise Exception("Unrecoginized types for filename and results")

    filename = expandvars(filename)

    extension = os.path.splitext(filename)[-1]

    if extension == '.yaml':
        open(filename, 'w').write(yaml.dump(tolist(results), **yaml_kwargs))
    elif extension == '.hdf5':
        if not isinstance(results, dict): raise Exception("Can only save dicts to hdf5 format.")
        import h5py
        f=h5py.File(filename,'w')
        for k,v in results.items(): f[k] = v
        f.close()
    elif extension == '.fits':
        if not isinstance(results, dict): raise Exception("Can only save dicts to fits format.")
        rec = fromarrays(results.values(), names=results.keys())
        makefits(rec, filename, clobber=True)
    elif extension == '.xml':
        from pyxml2obj import XMLout
        open(filename, 'w').write(XMLout(results))
    else:
        raise Exception("Unrecognized extension %s" % extension)
