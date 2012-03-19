import os

import yaml

import h5py

from . tools import tolist

def savedict(filename, results):
    """ Save a dictionary to a file. """
    is_results = lambda x: isinstance(x,dict) or isinstance(x,list)

    if isinstance(filename,str) and is_results(results):
        pass
    elif is_results(filename) and isinstance(results, str):
        filename, results = results, filename
    else:
        raise Exception("Unrecoginized types for filename and results")

    extension = os.path.splitext(filename)[-1]

    if extension == '.yaml':
        open(filename, 'w').write(yaml.dump(tolist(results)))
    elif extension == '.hdf5':
        if not isinstance(results, dict): raise Exception("Can only save dicts to hdf5 format.")
        f=h5py.File(filename,'w')
        for k,v in results.items(): f[k] = v
        f.close()
    else:
        raise Exception("Unrecognized extension %s" % extension)


