""" Module to convert from gtlike to pointlike spectral models (and vice-verca).
"""
import numpy as np

import pyLikelihood

from uw.like.Models import FileFunction
from uw.utilities.parmap import LimitMapper,LinearMapper
from uw.utilities.xml_parsers import XML_to_Model

from . save import spectrum_to_dict
from . load import dict_to_spectrum

_funcFactory = pyLikelihood.SourceFactory_funcFactory()

def gtlike_unscale_all_parameters(spectrum):
    """ kind of a kluge. """
    return dict_to_spectrum(spectrum_to_dict(spectrum))


def build_gtlike_spectrum(model):
    """ Convert a pointlike uw.like.Models.Model to a pyLikelihood spectral model object. 

        Convert a pointlike model to a gtlike model:

            >>> import pyLikelihood
            >>> from uw.like.Models import PowerLaw
            >>> pointlike_model = PowerLaw()

        You cannot convert models unless they have limits on all parameters:

            >>> spectrum = build_gtlike_spectrum(pointlike_model)
            Traceback (most recent call last):
                ...
            Exception: Unable to build gtlike model. Parameter Norm must have limits (current mapper = <class 'uw.utilities.parmap.LogMapper'>).

        Set Norm params:

            >>> pointlike_model.setp('Norm', 1e-8)
            >>> pointlike_model.set_limits('Norm', 1e-10, 1e-6, scale=1e-9)
            >>> pointlike_model.set_free('Norm', False)
            >>> pointlike_model.set_error('Norm', 2e-9)

        Set index params:

            >>> pointlike_model.setp('index', 1.1)
            >>> pointlike_model.set_limits('index', -2, 3)
            >>> pointlike_model.set_error('index',0.25)

        After setting limits, it should work

            >>> spectrum = build_gtlike_spectrum(pointlike_model)


        check prefactor:

            >>> param = spectrum.getParam('Prefactor')
            >>> param.isFree()
            False
            >>> np.allclose(param.getTrueValue(),1e-8)
            True
            >>> param.getScale()
            1e-09
            >>> param.getBounds()[0]*param.getScale()
            1e-10
            >>> param.getBounds()[1]*param.getScale()
            1e-06
            >>> np.allclose(param.error()*param.getScale(),2e-9)
            True

        Check index:

            >>> param=spectrum.getParam('Index')
            >>> gvalue=param.getTrueValue()
            >>> np.allclose(gvalue,-1.1)
            True

            >>> lower,upper=param.getBounds()
            >>> lower,upper=sorted([lower*param.getScale(),upper*param.getScale()])
            >>> np.allclose([lower,upper], [-3, 2])
            True

            >>> spectrum.getParam('Index').isFree() == pointlike_model.get_free('Index')
            True
            >>> np.allclose(spectrum.getParam('Index').error(), 0.25)
            True
        
        Test spectral values:

            >>> energies = np.logspace(1, 6, 10000)
            >>> from uw.darkmatter.spectral import DMFitFunction
            >>> np.allclose(DMFitFunction.call_pylike_spectrum(spectrum, energies),
            ...     pointlike_model(energies), rtol=1e-20) 
            True

        FileFunction was previously buggy, but now works:

            >>> from tempfile import NamedTemporaryFile
            >>> temp = NamedTemporaryFile()
            >>> filename = temp.name
            >>> pointlike_model.save_profile(filename, emin=1, emax=1e6)

            >>> ff_pointlike = FileFunction(normalization=9.5,file=filename, set_default_oomp_limits=True)
            >>> ff_gtlike = build_gtlike_spectrum(ff_pointlike)
            >>> np.allclose(DMFitFunction.call_pylike_spectrum(ff_gtlike, energies),
            ...     ff_pointlike(energies), rtol=1e-20) 
            True
    """
    for p in model.param_names:
        if not isinstance(model.get_mapper(p),LimitMapper):
            raise Exception("Unable to build gtlike model. Parameter %s must have limits (mapper=%s)." % (p,model.get_mapper(p)))


    gtlike_name = model.gtlike['name']

    if gtlike_name != 'FileFunction':
        spectrum = _funcFactory.create(gtlike_name)
    else:
        spectrum = pyLikelihood.FileFunction()
        spectrum.readFunction(model.file)

    for p,g in zip(model.param_names,model.gtlike['param_names']):
        param=spectrum.getParam(g)

        scale=model.get_scale_gtlike(p)
        param.setScale(scale)

        lower,upper=model.get_limits_gtlike(p)
        lower,upper=sorted([lower/scale,upper/scale])

        # NB, most robust to set scale, then value, then bounds
        param.setTrueValue(model.getp_gtlike(p))
        param.setBounds(lower,upper)

        param.setFree(model.get_free(p))
        param.setError(abs(model.error(p)/scale))

    for p,g in model.gtlike['extra_param_names'].items():
        param=spectrum.getParam(g)
        param.setScale(1)
        param.setTrueValue(model[p])
        param.setBounds(model[p],model[p])

    return spectrum


def build_pointlike_model(spectrum):
    """ Convert a gtlike model object to a pointlike
        model object.
        
            >>> spectrum = _funcFactory.create('PowerLaw')

            >>> param=spectrum.getParam('Prefactor')
            >>> param.setScale(10)
            >>> param.setTrueValue(1e-9)
            >>> param.setBounds(1e-11,1e-9)
            >>> param.setError(3e-11)
            >>> param.setFree(True)

            >>> param=spectrum.getParam('Index')
            >>> param.setScale(2)
            >>> param.setBounds(-10,5)
            >>> param.setTrueValue(-3)
            >>> param.setError(0.125)
            >>> param.setFree(False)

        Check spectral values:

            >>> model = build_pointlike_model(spectrum)
            >>> energies = np.logspace(1, 6, 10000)
            >>> from uw.darkmatter.spectral import DMFitFunction
            >>> np.allclose(DMFitFunction.call_pylike_spectrum(spectrum, energies),
            ...     model(energies), rtol=1e-20, atol=1e-20) 
            True

        Check prefactor:

            >>> model.get_scale('norm')
            10.0
            >>> np.allclose(model.get_limits('norm'),[1e-10, 1e-08])
            True
            >>> np.allclose(model.getp('norm'),1e-9)
            True
            >>> np.allclose(model.error('norm'),1e-10)
            True
            >>> model.get_free('norm')
            True

        Check index params:

            >>> model.get_scale('index')
            -2.0
            >>> model.get_limits('index')
            [-10.0, 20.0]
            >>> model.getp('index')
            3.0
            >>> model.error('index')
            0.25
            >>> model.get_free('index')
            False

        Example creating a FileFunction object:
        
        First, create file out of old model:

            >>> from tempfile import NamedTemporaryFile
            >>> temp = NamedTemporaryFile()
            >>> filename = temp.name
            >>> model.save_profile(filename, emin=1, emax=1e6)
        
        Now, make FileFunction:

            >>> spectrum = pyLikelihood.FileFunction()
            >>> spectrum.readFunction(filename)

        Set param values:

            >>> param=spectrum.getParam('Normalization')
            >>> param.setScale(2)
            >>> param.setTrueValue(4)
            >>> param.setBounds(.1,10)

            >>> model = build_pointlike_model(spectrum)

        Test spectral points:

            >>> np.allclose(DMFitFunction.call_pylike_spectrum(spectrum, energies),
            ...     model(energies), rtol=1e-20, atol=1e-20) 
            True

        Test param values:

            >>> model.get_scale('Normalization')
            2.0
            >>> model.getp('Normalization')
            4.0
            >>> model.get_limits('Normalization')
            [0.2, 20.0]


    """
    gtlike_name = spectrum.genericName()

    if gtlike_name == 'FileFunction':
        ff=pyLikelihood.FileFunction_cast(spectrum)
        filename=ff.filename()
        model = FileFunction(file=filename)
    else:
        model = XML_to_Model.modict[gtlike_name]()
    
    param_names = pyLikelihood.StringVector()
    spectrum.getParamNames(param_names)
    for gtlike_name in param_names:
        pointlike_name = model.get_pointlike_name(gtlike_name)
        
        param=spectrum.getParam(gtlike_name)

        if pointlike_name in model.default_extra_params.keys():
            # no mapping for extra params
            model.setp(pointlike_name,param.getTrueValue())
        else:
            model.setp_gtlike(pointlike_name,param.getTrueValue())

            if pointlike_name in model.param_names:
                model.set_mapper(pointlike_name, LinearMapper)
                model.set_limits_gtlike(
                    pointlike_name,
                    lower=param.getBounds()[0]*param.getScale(),
                    upper=param.getBounds()[1]*param.getScale(),
                    scale=param.getScale())
                model.set_error(
                    pointlike_name,
                    abs(param.error()*param.getScale()))
                model.set_free(pointlike_name, param.isFree())
    return model


if __name__ == "__main__":
    import doctest
    doctest.testmod()
