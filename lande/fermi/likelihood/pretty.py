import sys

from uw.like.Models import Model,PowerLaw
from uw.like.SpatialModels import SpatialModel,Gaussian,SpatialMap
from uw.like.pointspec_helpers import PointSource 
from uw.like.roi_extended import ExtendedSource

def f(x):
    return pprint.pformat(x)

def pformat_model(model):
    """ Convert a spectral model to a python string which could be used
        to recate the object.
    """
    assert isinstance(model,Model)
    return model.name + '(' + \
            ', '.join('%s=%s' % (m.lower(),f(model[m])) for m in model.param_names + model.default_extra_params.keys()) + \
            ')'

def pformat_spatial_model(model):
    """ Same as model_to_string, but for SpatialModel objects
    """
    assert isinstance(model,SpatialModel)
    # For now, special case for SpatialMap:
    if isinstance(model,SpatialMap):
        return "SpatialMap(file=%s)" % f(model.file)
    else:
        return model.name + '(' + \
                ', '.join('%s=%s' % (m.lower(),f(model[m])) for m in model.param_names) + \
                ')'

def pformat_point_source(ps):
    pass

def pformat_extended_source(ps):
    pass

def pformat_source(source):
    pass

def pformat(model):
    """ Intellegently converts to a string any of the allowed
        objects.

        Example, pformat a spectral model

            >>> m=PowerLaw(norm=8.81e-14, index=2.63, e0=3162)
            >>> print pformat(m)
            PowerLaw(norm=8.81e-14, index=2.63, e0=3162)

        formatting a FileFunction is a little trickier:

            >>> from tempfile import NamedTemporaryFile
            >>> from uw.like.Models import FileFunction
            >>> tempfile = NamedTemporaryFile()
            >>> m.save_profile(tempfile.name, 10, 1e6)
            >>> m=FileFunction(normalization=1.0, file=tempfile.name)
            >>> print pformat(m).replace(tempfile.name,'[FILENAME]')
            FileFunction(normalization=1.0, file='[FILENAME]')

        pformat a spatial model

            >>> m=Gaussian(sigma=1, l=2, b=3)
            >>> print pformat(m)
            Gaussian(l=2.0, b=3.0, sigma=1.0)

        pformat a SpatialMap (tricky because a filename)

            >>> tempfile = NamedTemporaryFile()
            >>> m.save_template(tempfile.name)
            >>> m=SpatialMap(file=tempfile.name)
            >>> print pformat(m).replace(tempfile.name,'[FILENAME]')
            SpatialMap(file='[FILENAME]')

    """
    if isinstance(model,Model):
        return pformat_model(model)
    elif isinstance(model,SpatialModel):
        return pformat_spatial_model(model)
    elif isinstance(model,PointSource) or isinstance(model,ExtendedSource):
        return pformat_source(model)
    else:
        raise Exception("Unrecognized type %s" % type(model))

def pprint(model,out=sys.stdout):
    out.write(pformat(model))

if __name__ == "__main__":
    import doctest
    doctest.testmod()
