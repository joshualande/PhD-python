from uw.utilities.parmap import LimitMapper

import numpy as np

import pyLikelihood

_funcFactory = pyLikelihood.SourceFactory_funcFactory()

def build_gtlike_model(model):
    """ Convert a pointlike uw.like.Models.Model to a pyLikelihood spectral model object. 

        Convert a pointlike model to a gtlike model:

            >>> import pyLikelihood
            >>> from uw.like.Models import PowerLaw
            >>> pointlike_model = PowerLaw()

        Add in some weird stuff to the model

            >>> pointlike_model.set_limits('index', -2, 3)
            >>> pointlike_model.set_free('Norm', False)
        
        You cannot convert models unless they have limits on all parameters:

            >>> gtlike_model = build_gtlike_model(pointlike_model)
            Traceback (most recent call last):
                ...
            Exception: Unable to build gtlike model. Parameter Norm must have limits.

        After setting limits, it should work

            >>> pointlike_model.set_default_limits()
            >>> gtlike_model = build_gtlike_model(pointlike_model)

        Check to see if values & limits in model are correct:
            >>> param=gtlike_model.getParam('Index')
            >>> gvalue=param.getTrueValue()
            >>> pvalue=pointlike_model.gtlike['togtlike'][pointlike_model.name_mapper('Index')](pointlike_model['Index'])
            >>> np.allclose(gvalue,pvalue)
            True

            >>> lower,upper=param.getBounds()
            >>> lower,upper=lower*param.getScale(),upper*param.getScale()
            >>> np.allclose([lower,upper], pointlike_model.get_limits('Index'))
            True

        Test free:
            >>> gtlike_model.getParam('Prefactor').isFree() == pointlike_model.get_free('Norm')
            True
            >>> gtlike_model.getParam('Index').isFree() == pointlike_model.get_free('Index')


            >>> energies = np.logspace(1, 6, 10000)
            >>> from uw.darkmatter.spectral import DMFitFunction
            >>> np.allclose(DMFitFunction.call_pylike_spectrum(gtlike_model, energies),
            ...     pointlike_model(energies), rtol=1e-20, atol=1e-20) 
            True

    """
    for p in model.param_names:
        if not isinstance(model.get_mapper(p),LimitMapper):
            raise Exception("Unable to build gtlike model. Parameter %s must have limits." % p)


    gtlike_name = model.gtlike['name']

    assert gtlike_name != 'FileFunction'

    gtlike_model = _funcFactory.create(gtlike_name)


    for p,g in zip(model.param_names,model.gtlike['param_names']):

        param=gtlike_model.getParam(g)

        scale=model.get_scale_gtlike(p)
        param.setScale(scale)

        lower,upper=model.get_limits_gtlike(p)
        lower,upper=lower/scale,upper/scale
        param.setBounds(lower,upper)

        param.setTrueValue(model.getp_gtlike(p))

        param.setFree(model.get_free(p))

    for p,g in model.gtlike['extra_param_names'].items():
        gtlike_model.setParam(g,model[p])

    return gtlike_model

if __name__ == "__main__":
    import doctest
    doctest.testmod()
