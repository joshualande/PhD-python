import numpy as np

import pyLikelihood

_funcFactory = pyLikelihood.SourceFactory_funcFactory()

def build_gtlike_model(model):
    """ Convert a pointlike uw.like.Models.Model to a pyLikelihood spectral model object. 

        Convert a pointlike model to a gtlike model:
            >>> import pyLikelihood
            >>> from uw.like.Models import PowerLaw
            >>> pointlike_model = PowerLaw()
            >>> gtlike_model = build_gtlike_model(pointlike_model)

            >>> energies = np.logspace(1, 6, 10000)
            >>> from uw.darkmatter.spectral import DMFitFunction
            >>> np.allclose(DMFitFunction.call_pylike_spectrum(gtlike_model, energies),
            ...     pointlike_model(energies), rtol=1e-20, atol=1e-20) 
            True

    """
    gtlike_name = model.gtlike['name']

    assert gtlike_name != 'FileFunction'

    gtlike_model = _funcFactory.create(gtlike_name)


    for p,g in zip(model.param_names,model.gtlike['param_names']):
        gtlike_model.setParam(g,model.getp_gtlike(p))

    for p,g in model.gtlike['extra_param_names'].items():
        gtlike_model.setParam(g,model[p])

    return gtlike_model

if __name__ == "__main__":
    import doctest
    doctest.testmod()
