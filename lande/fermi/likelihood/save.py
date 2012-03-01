import traceback
import sys

import pylab as P
import numpy as np

from skymaps import DiffuseFunction,IsotropicSpectrum,IsotropicPowerLaw

from uw.like.pointspec_helpers import PointSource
from uw.like.roi_extended import ExtendedSource
from uw.like.roi_diffuse import DiffuseSource
from uw.like.Models import DefaultModelValues

from SED import SED

from lande.pysed import units
from lande.utilities.tools import tolist
from . tools import gtlike_or_pointlike


def gtlike_get_full_energy_range(like): return like.energies[[0,-1]]
def pointlike_get_full_energy_range(roi): return roi.bin_edges[[0,-1]]
def get_full_energy_range(*args, **kwargs):
    return gtlike_or_pointlike(gtlike_get_full_energy_range, pointlike_get_full_energy_range, *args, **kwargs)


def pointlike_spectrum_to_dict(model, errors=False):
    """ Package of a spectral model into a handy
        python dictionary.

            >>> from uw.like.Models import PowerLaw
            >>> m=PowerLaw(norm=1, index=.5, index_offset=1)
            >>> d=pointlike_spectrum_to_dict(m)
            >>> print d['Norm']
            1.0
            >>> print d['Index']
            -0.5
            >>> d=pointlike_spectrum_to_dict(m, errors=False)
            >>> d.has_key('Index_err')
            False
    """
    d = dict(name = model.name)
    default = DefaultModelValues.simple_models[model.name]
    for k,v in default.items():
        if k == '_p': continue
        elif k == 'param_names': 
            for p in v: 
                d[p]=model[p]
                if errors:
                    d['%s_err' % p]=model.error(p)
        else: 
            d[k]=getattr(model,k)

    # stupid kluge:
    if d.has_key('index_offset'):
        index_offset = d.pop('index_offset')
        d['Index'] -= index_offset
        if errors:
            d['Index_err'] -= index_offset

    return tolist(d)

gtlike_spectrum_to_dict = SED.spectrum_to_dict


def gtlike_name_to_dict(like, name, *args, **kwargs):
    source = like.logLike.getSource(name)
    spectrum = source.spectrum()
    return gtlike_spectrum_to_dict(spectrum, *args, **kwargs)

def pointlike_name_to_dict(roi, name, *args, **kwargs):
    model = roi.get_model(name)
    return pointlike_spectrum_to_dict(model, *args, **kwargs)


def gtlike_fluxdict(like,name,emin=None,emax=None,flux_units='erg', error=True):

    if emin is None and emax is None: 
        emin, emax = get_full_energy_range(like)

    ce=lambda e: units.convert(e,'MeV',flux_units)
    f=dict(flux=like.flux(name,emin=emin,emax=emax),
           eflux=ce(like.energyFlux(name,emin=emin,emax=emax)),
           flux_units='ph/cm^2/s',
           eflux_units='%s/cm^2/s' % flux_units,
           emin=emin,
           emax=emax)

    if error:
        try:
            # incase the errors were not calculated
            f['flux_err']=like.fluxError(name,emin=emin,emax=emax)
            f['eflux_err']=ce(like.energyFluxError(name,emin=emin,emax=emax))
        except Exception, ex:
            print 'ERROR calculating flux error: ', ex
            traceback.print_exc(file=sys.stdout)
            f['flux_err']=-1
            f['eflux_err']=-1
    return tolist(f)


def gtlike_get_spatial_model_name(like, name):
    """ This code is adapted from
        the Likelihood file SourceModelBuilder.cxx's
        function
        SourceModelBuilder::addSpatialPart
    """
    source = like.logLike.getSource(name)

    fns = source.getSrcFuncs()

    assert fns.count("Position") or fns.count("SpatialDist")

    if fns.count("Position"):
        return "SkyDirFunction"

    elif fns.count("SpatialDist"):
        type = fns["SpatialDist"].genericName()
        return type

def pointlike_get_spatial_model_name(roi, name):
    source = roi.get_source(name)
    if isinstance(source,PointSource):
        return 'SkyDirFunction'
    elif isinstance(source,ExtendedSource):
        # this is only approximately true
        return 'SpatialMap'
    elif isinstance(source,DiffuseSource):
        dm = source.dmodel
        if hasattr(dm,'__len__') and len(dm)==1: dm = dm[0]
        if isinstance(dm,DiffuseFunction):
            return 'MapCubeFunction'
        elif isinstance(dm,IsotropicSpectrum) or \
                isinstance(dm,IsotropicPowerLaw):
            return 'ConstantValue'

def get_spatial_model_name(*args, **kwargs):
    return gtlike_or_pointlike(gtlike_get_spatial_model_name, pointlike_get_spatial_model_name, *args, **kwargs)

def pointlike_get_all_names(roi):
    """ Get a list of the names of all sources in the pointlike ROI. """
    return np.append(roi.psm.names,roi.dsm.names)

def gtlike_get_all_names(like):
    """ Get a list of the names of all sources in the gtlike ROI. """

def get_all_names(*args, **kwargs):
    return gtlike_or_pointlike(gtlike_get_all_names, pointlike_get_all_names, *args, **kwargs)
    return like.sourceNames()


def get_sources(like_or_roi):
    """ Get a list of point-like and extended sources
        in the ROI. """
    all_names=get_all_names(like_or_roi)
    all_ps = [i for i in all_names \
              if get_spatial_model_name(like_or_roi,i) in \
              ['SkyDirFunction','SpatialMap']]
    return all_ps
        

def get_background(like_or_roi):
    """ Get a list of the names of all background
        sources in the ROI. """
    all_names=get_all_names(like_or_roi)
    all_bg = [i for i in all_names \
              if get_spatial_model_name(like_or_roi,i) in \
              ['ConstantValue','MapCubeFunction']]
    return all_bg

def diffusedict(like_or_roi):
    """ Save out all diffuse sources. """

    f = dict()

    bgs = get_background(like_or_roi)

    for name in bgs:
        f[name] = name_to_dict(like_or_roi, name, errors=True)
    return tolist(f)

def gtlike_sourcedict(like, name, emin=None, emax=None, flux_units='erg'):
    from pyLikelihood import ParameterVector

    if emin is None and emax is None:
        emin, emax = get_full_energy_range(like)

    d=dict(
        TS=like.Ts(name,reoptimize=True),
        logLikelihood=like.logLike.value()
    )

    d['flux']=fluxdict(like,name,emin,emax,flux_units=flux_units)

    source = like.logLike.getSource(name)
    spectrum = source.spectrum()

    d['model']=spectrum_to_dict(spectrum, errors=True)

    d['diffuse'] = diffusedict(like)

    return tolist(d)

def sourcedict(*args, **kwargs):
    return gtlike_or_pointlike(gtlike_sourcedict, pointlike_sourcedict, *args, **kwargs)

def pointlike_model_to_flux(model, emin, emax, flux_units='erg', error=True):

    ce=lambda e: units.convert(e,'MeV',flux_units)
    f=dict()
    if error:
        f['flux'],f['flux_err']=model.i_flux(emin=emin,emax=emax,error=True)
        ef,ef_err=model.i_flux(emin=emin,emax=emax,e_weight=1,error=True)
        f['eflux'],f['eflux_err']=ce(ef),ce(ef_err)
    else:
        f['flux']=model.i_flux(emin=emin,emax=emax,error=False)
        ef=model.i_flux(emin=emin,emax=emax,e_weight=1,error=False)
        f['eflux']=ce(ef)

    f['flux_units']='ph/cm^2/s'
    f['eflux_units']='%s/cm^2/s' % flux_units
    f['emin'],f['emax']=emin,emax
    return tolist(f)

def pointlike_fluxdict(roi, which, emin=None, emax=None, *args, **kwargs):

    if emin is None and emax is None:
        emin, emax = get_full_energy_range(roi)

    model=roi.get_model(which)
    return tolist(pointlike_model_to_flux(model, emin, emax, *args, **kwargs))



def pointlike_sourcedict(roi, name, emin=None, emax=None, flux_units='erg'):
    d={}

    if emin is None and emax is None:
        emin, emax = get_full_energy_range(roi)

    source=roi.get_source(name)
    model=roi.get_model(name)

    old_quiet = roi.quiet; roi.quiet=True
    d['TS'] = roi.TS(name,quick=False)
    roi.quiet = old_quiet

    d['logLikelihood']=-roi.logLikelihood(roi.parameters())

    d['flux']=fluxdict(roi,name,emin,emax,flux_units)

    d['model']=spectrum_to_dict(model, errors=True)

    # Source position
    d['gal'] = [source.skydir.l(),source.skydir.b()]
    d['equ'] = [source.skydir.ra(),source.skydir.dec()]

    d['diffuse'] = diffusedict(roi)

    f = d['spatial_model'] = dict()
    if isinstance(source,ExtendedSource):
        # Extended Source parameters
        spatial_model = source.spatial_model
        for param in spatial_model.param_names:
            f[param]=spatial_model[param]
            f[param + '_err']=spatial_model.error(param)

    # add elliptical error, if they exist.
    # N.B. If no localization performed, this will return
    # an empty dictionary.
    # N.B. This method will do the wrong thing if you have recently relocalized
    # another source. This is rarely the case.
    f.update(roi.get_ellipse())

    return tolist(d)

def spectrum_to_dict(*args, **kwargs):
    return gtlike_or_pointlike(gtlike_spectrum_to_dict, pointlike_spectrum_to_dict, *args, **kwargs)

def name_to_dict(*args, **kwargs):
    return gtlike_or_pointlike(gtlike_name_to_dict, pointlike_name_to_dict, *args, **kwargs)


def fluxdict(*args, **kwargs):
    return gtlike_or_pointlike(gtlike_fluxdict, pointlike_fluxdict, *args, **kwargs)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
