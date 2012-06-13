from uw.like.Models import Model,PowerLaw
from uw.like.SpatialModels import SpatialModel,Gaussian

def model_to_string(model):
    """ Convert a spectral model to a python string which could be used
        to recate the object.

            >>> m=PowerLaw(norm=8.81e-14, index=2.63, e0=3162)
            >>> print tostring(m)
            PowerLaw(norm=8.81e-14, index=2.63, e0=3162)
    """
    assert isinstance(model,Model)
    return model.name + '(' + \
            ', '.join('%s=%s' % (m.lower(),model[m]) for m in model.param_names + model.default_extra_params.keys()) + \
            ')'

def spatial_model_to_string(model):
    """ Same as model_to_string, but for SpatialModel objects

            >>> m=Gaussian(sigma=1, l=2, b=3)
            >>> print spatial_model_to_string(m)
            Gaussian(sigma=1, l=2, b=3)
            >>> 
    """
    assert isinstance(model,SpatialModel)
    return model.name + '(' + \
            ', '.join('%s=%s' % (m.lower(),model[m]) for m in model.param_names) + \
            ')'


def point_source_to_model(ps):
    pass

def extended_source_to_model(ps):
    pass

def soruce_to_model(source):
    pass

def tostring(model):
    """ Intellegently converts to a string the object.
    if isinstance(model,Model):
        return model_to_string(model)
    elif isinstance(model,SpatialModel):
        return spatial_model_to_string(model)
    elif isinstance(model,PointSource) or isinstance(model,ExtendedSource):
        return soruce_to_model(model)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
