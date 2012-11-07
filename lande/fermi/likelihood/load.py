import numpy as np

import pyLikelihood

from uw.like.Models import FileFunction
import uw.like.Models
from uw.darkmatter.spectral import DMFitFunction

_funcFactory = pyLikelihood.SourceFactory_funcFactory()


def pointlike_dict_to_spectrum(d):
    if d['name'] == 'FileFunction':
        model = FileFunction(file=d['file'])
    elif d['name'] == 'DMFitFunction':
        model = DMFitFunction()
    else:
        model = uw.like.Models.__dict__[d['name']]()
    for k,v in d.items(): 
        if k not in ['name','method','file','covariance_matrix']:
            if len(k) > 10 and k[-10:] in ['_upper_err','_lower_err']:
                pass
            elif len(k) > 4 and k[-4:] == '_err': 
                model.set_error(k[:-4],v)
            else:
                model[k]=v
    return model
    
def gtlike_dict_to_spectrum(d):
    """ Load back as a pyLikelihood spectrum object
        a spectrum that has been saved by the spectrum_to_string
        object. This undoes the conversion of spectrum_to_dict """
    if d['name'] == 'FileFunction':
        spectrum = pyLikelihood.FileFunction()
        spectrum.readFunction(d['file'])
    else:
        spectrum=_funcFactory.create(d['name'])
    for k,v in d.items(): 
        if k not in ['name','method','file','covariance_matrix']:
            if len(k) > 10 and k[-10:] in ['_upper_err','_lower_err']:
                pass
            elif len(k) > 4 and k[-4:] == '_err': 
                param =  spectrum.getParam(k[:-4])
                param.setError(v/param.getScale())
            else:
                spectrum.getParam(k).setTrueValue(v)
    return spectrum

def dict_to_spectrum(d):
    """ Test out FileFunction:

            >>> from uw.like.Models import PowerLaw
            >>> pl= PowerLaw()
            >>> from tempfile import NamedTemporaryFile
            >>> tempfile = NamedTemporaryFile()
            >>> filename = tempfile.name
            >>> pl.save_profile(filename, 1, 100)

            >>> d = dict(file=filename, method='pointlike', Normalization=3, name='FileFunction')
            >>> model = dict_to_spectrum(d)
            >>> isinstance(model,FileFunction)
            True
            >>> model.file == filename
            True
            >>> np.allclose(model['Normalization'],3)
            True

            >>> d = dict(file=filename, method='gtlike', Normalization=3, name='FileFunction')
            >>> spectrum = dict_to_spectrum(d)
            >>> isinstance(spectrum,pyLikelihood.FileFunction)
            True
            >>> spectrum.filename() == filename
            True
            >>> np.allclose(spectrum.getParam('Normalization').getTrueValue(),3)
            True
    """
    assert d['method'] in ['gtlike','pointlike']
    if d['method'] == 'gtlike':
        return gtlike_dict_to_spectrum(d)
    if d['method'] == 'pointlike':
        return pointlike_dict_to_spectrum(d)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
