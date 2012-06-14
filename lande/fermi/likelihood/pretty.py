import sys

from skymaps import SkyDir
from uw.like.Models import Model,PowerLaw
from uw.like.SpatialModels import SpatialModel,Gaussian,SpatialMap
from uw.like.pointspec_helpers import PointSource 
from uw.like.roi_extended import ExtendedSource

def f(x):
    if isinstance(x,str):
        return "'%s'" % x
    else:
        return str(x)

def pformat_skydir(skydir,galactic=True):
    if galactic:
        return 'SkyDir(%s,%s,SkyDir.GALACTIC)' % (skydir.l(),skydir.b())
    else:
        return 'SkyDir(%s,%s)' % (skydir.ra(),skydir.dec())

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
    return 'PointSource(name=%s, skydir=%s, model=%s)' % (f(ps.name),
                                                          pformat_skydir(ps.skydir),
                                                          pformat_model(ps.model))

def pformat_extended_source(ps):
    return 'ExtendedSource(name=%s, spatial_model=%s, model=%s)' % (f(ps.name),
                                                                    pformat_spatial_model(ps.spatial_model),
                                                                    pformat_model(ps.model))

def pformat_source(source):
    if isinstance(source,PointSource):
        return pformat_point_source(source)
    elif isinstance(source,ExtendedSource):
        return pformat_extended_source(source)
    else:
        raise Exception("Unrecognized type %s" % type(source))

def pformat(model, *args, **kwargs):
    """ Intellegently converts to a string any of the allowed
        objects.

        Example, pformat a skydir:

            >>> sd=SkyDir(0,0)
            >>> print pformat(sd,galactic=False)
            SkyDir(0.0,0.0)

            >>> sd=SkyDir(0,0,SkyDir.GALACTIC)
            >>> print pformat(sd)
            SkyDir(6.36110936293e-15,-1.59027734073e-15,SkyDir.GALACTIC)
        
        pformat a spectral model:

            >>> m=PowerLaw(norm=8.81e-14, index=2.63, e0=3162)
            >>> print pformat(m)
            PowerLaw(norm=8.81e-14, index=2.63, e0=3162)

        formatting a FileFunction is a little trickier:

            >>> from tempfile import NamedTemporaryFile
            >>> from uw.like.Models import FileFunction
            >>> tempfile = NamedTemporaryFile()
            >>> m.save_profile(tempfile.name, 10, 1e6)
            >>> ff=FileFunction(normalization=1.0, file=tempfile.name)
            >>> print pformat(ff).replace(tempfile.name,'[FILENAME]')
            FileFunction(normalization=1.0, file='[FILENAME]')

        pformat a spatial model

            >>> sm=Gaussian(sigma=1, l=2, b=3)
            >>> print pformat(sm)
            Gaussian(l=2.0, b=3.0, sigma=1.0)

        pformat a SpatialMap (tricky because a filename)

            >>> tempfile = NamedTemporaryFile()
            >>> sm.save_template(tempfile.name)
            >>> map=SpatialMap(file=tempfile.name)
            >>> print pformat(map).replace(tempfile.name,'[FILENAME]')
            SpatialMap(file='[FILENAME]')

        pformat a PointSource:

            >>> ps = PointSource(name='source',skydir=sd, model=m)
            >>> print pformat(ps)
            PointSource(name='source', skydir=SkyDir(6.36110936293e-15,-1.59027734073e-15,SkyDir.GALACTIC), model=PowerLaw(norm=8.81e-14, index=2.63, e0=3162))

            >>> es = ExtendedSource(name='source',spatial_model=sm, model=m)
            >>> print pformat(es)
            ExtendedSource(name='source', spatial_model=Gaussian(l=2.0, b=3.0, sigma=1.0), model=PowerLaw(norm=8.81e-14, index=2.63, e0=3162))
    """
    if isinstance(model,SkyDir):
        func=pformat_skydir
    elif isinstance(model,Model):
        func=pformat_model
    elif isinstance(model,SpatialModel):
        func=pformat_spatial_model
    elif isinstance(model,PointSource) or isinstance(model,ExtendedSource):
        func=pformat_source
    else:
        raise Exception("Unrecognized type %s" % type(model))
    return func(model, *args, **kwargs)

def pprint(model,out=sys.stdout):
    out.write(pformat(model))

if __name__ == "__main__":
    import doctest
    doctest.testmod()
