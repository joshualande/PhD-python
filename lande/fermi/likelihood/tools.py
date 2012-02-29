
""" This file contains various function which I have found useful. """
import pprint
from math import ceil
import traceback
import sys

import pylab as P
import numpy as np
import pyfits as pf

from mpl_toolkits.axes_grid1 import AxesGrid

from skymaps import DiffuseFunction,IsotropicSpectrum,IsotropicPowerLaw

from uw.like.roi_analysis import ROIAnalysis
from uw.like.pointspec_helpers import PointSource
from uw.like.roi_extended import ExtendedSource
from uw.like.roi_diffuse import DiffuseSource
from uw.like.Models import Model,PowerLaw,ExpCutoff,DefaultModelValues
from uw.like.roi_state import PointlikeState

from lande.utilities.toolbag import tolist
from SED import SED

from lande.pysed import units

from . superstate import SuperState


def gtlike_or_pointlike(f_gtlike, f_pointlike, like_or_roi, *args, **kwargs):
    """ Note, like_or_roi can be either an ROI or a spectral model. """

    from pyLikelihood import Function
    from BinnedAnalysis import BinnedAnalysis

    if isinstance(like_or_roi, BinnedAnalysis) or \
       isinstance(like_or_roi,Function):
        f=f_gtlike
    elif isinstance(like_or_roi, ROIAnalysis) or \
       isinstance(like_or_roi,Model):
        f=f_pointlike
    else:
        raise Exception("like_or_roi must be of type BinnedAnalysis or ROIAnalysis")
    return f(like_or_roi, *args, **kwargs)


def gtlike_get_full_energy_range(like): return like.energies[[0,-1]]
def pointlike_get_full_energy_range(roi): return roi.bin_edges[[0,-1]]


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

def pointlike_get_all_names(roi):
    """ Get a list of the names of all sources in the pointlike ROI. """
    return np.append(roi.psm.names,roi.dsm.names)

def gtlike_get_all_names(like):
    """ Get a list of the names of all sources in the gtlike ROI. """
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

def gtlike_modify(like, name, free=True):
    """ Freeze a source in a gtlike ROI. 
    
        The method for modifying the ROI
        follows the code in SuperState.py 
        
        I am not sure why the modificaiton
        has to be done in this particular way. """
    from pyLikelihood import StringVector

    source = like.logLike.getSource(name)
    spectrum = like[name].src.spectrum()

    parNames = StringVector()
    spectrum.getParamNames(parNames)
    for parName in parNames:
        index = like.par_index(name, parName)
        par = like.params()[index]
        par.setFree(free)

    like.syncSrcParams(name)


def gtlike_upper_limit(like, name, cl, emin=None, emax=None, 
                       flux_units='erg', **kwargs):
    """
        N.B. spectral fit in this function instead
        of in upper limits code since my
        paranoid_gtlike_fit function is more robust. """

    print 'Calculating gtlike upper limit'

    saved_state = SuperState(like)
    source = like.logLike.getSource(name)

    try:
        import IntegralUpperLimit

        if emin is None and emax is None: 
            emin, emax = get_full_energy_range(like)

        # First, freeze spectrum (except for normalization)
        # of our soruce
        gtlike_fit_only_prefactor(like, name)

        # Spectral fit whole ROI

        paranoid_gtlike_fit(like)

        # Freeze everything but our source of interest
        for i in range(len(like.model.params)):
            like.model[i].setFree(False)
            like.syncSrcParams(like[i].srcName)

        # Note, I think freeze_all is redundant, but flag it just 
        # to be paranoid
        flux_ul, results = IntegralUpperLimit.calc_int(like, name, 
                                                       freeze_all=True,
                                                       skip_global_opt=True,
                                                       cl=cl,
                                                       emin=emin, 
                                                       emax=emax, 
                                                       **kwargs)

        prefactor=like.normPar(name)
        pref_ul = results['ul_value']*prefactor.getScale()
        prefactor.setTrueValue(pref_ul)

        flux_ul = like.flux(name,emin,emax)
        flux_units_string = 'ph/cm^2/s'

        eflux_ul = units.convert(like.energyFlux(name,emin,emax), 'MeV', flux_units)
        eflux_units_string = '%s/cm^2/s' % flux_units

        ul = dict(
            emin=emin, emax=emax,
            flux_units=flux_units_string, flux=flux_ul, 
            eflux_units=eflux_units_string, eflux=eflux_ul)

    except Exception, ex:
        print 'ERROR gtlike upper limit: ', ex
        traceback.print_exc(file=sys.stdout)
        ul = None
    finally:
        saved_state.restore()

    return tolist(ul)


def gtlike_powerlaw_upper_limit(like, name, powerlaw_index=2 , cl=0.95, emin=None, emax=None, 
                                flux_units='erg',
                                **kwargs):
    """ Wrap up calculating the flux upper limit for a powerlaw
        source.  This function employes the pyLikelihood function
        IntegralUpperLimit to calculate a Bayesian upper limit.

        The primary benefit of this function is that it replaces the
        spectral model automatically with a PowerLaw spectral model
        and fixes the index to -2. It then picks a better scale for the
        powerlaw and gives the upper limit calculation a more reasonable
        starting value, which helps the convergence.
    """
    print 'Calculating gtlike power-law upper limit'

    saved_state = SuperState(like)

    if emin is None and emax is None: 
        emin, emax = get_full_energy_range(like)

    e = np.sqrt(emin*emax)

    # assume a canonical dnde=1e-11 at 1GeV index 2 starting value
    dnde = PowerLaw(norm=1e-11, index=2,e_scale=1e3)

    like.setSpectrum(name,'PowerLaw')

    # fix index to 0
    index=like[like.par_index(name, 'Index')]
    index.setTrueValue(-1*powerlaw_index)

    # good starting guess for source
    prefactor=like[like.par_index(name, 'Prefactor')]
    prefactor.setScale(dnde(e))
    prefactor.setValue(1)
    # unbound the prefactor since the default range 1e-2 to 1e2 may not be big enough
    # in small phase ranges.
    prefactor.setBounds(1e-10,1e10)

    scale=like[like.par_index(name, 'Scale')]
    scale.setScale(1)
    scale.setValue(e)

    like.syncSrcParams(name)

    results = gtlike_upper_limit(like, name, cl, emin, emax, flux_units, **kwargs)
    if results is not None:
        results['powerlaw_index']=powerlaw_index

    saved_state.restore()

    return tolist(results)

def pointlike_upper_limit(roi, name, cl, emin=None, emax=None, flux_units='erg', **kwargs):

    if emin is None and emax is None:
        emin, emax = get_full_energy_range(roi)

    params = roi.parameters().copy()
    try:
        flux_ul = roi.upper_limit(which=name, confidence=cl, emin=emin, emax=emax, **kwargs)

        flux_units_string = 'ph/cm^2/s'

        ul = dict(
            emin=emin, emax=emax,
            flux_units=flux_units_string, 
            flux=flux_ul)

    except Exception, ex:
        print 'ERROR pointlike upper limit: ', ex
        traceback.print_exc(file=sys.stdout)
        ul = None
    finally:
        roi.set_parameters(params)
        roi.__update_state__()

    return tolist(ul)


def pointlike_powerlaw_upper_limit(roi, name, powerlaw_index=2, cl=0.95, emin=None, emax=None, 
                                   flux_units='erg', **kwargs):
    print 'Calculating pointlike upper limit'

    saved_state = PointlikeState(roi)

    """ Note keep old flux, because it is important to have
        the spectral model pushed into the upper_limit
        code reasonably close to the best fit flux. This
        is because initial likelihood (ll_0) is used to scale
        the likelihood so it has to be reasonably close to 
        the best value. """
    roi.modify(which=name, model=PowerLaw(index=powerlaw_index), keep_old_flux=True)

    ul = pointlike_upper_limit(roi, name, cl, emin, emax, flux_units, **kwargs)
    ul['powerlaw_index']=powerlaw_index

    saved_state.restore()

    return tolist(ul)




def freeze_insignificant_to_catalog(roi,catalog,exclude_names=[], min_ts=25):
    """ Replace all insigificant 2FGL catalog sources
        with the predictions of 2FGL and 
        the spectral shape of the source frozen. """
    any_changed=False
    for source in roi.get_sources():
        name = source.name

        if name in exclude_names: continue

        # Note only check sources with MORE than their
        # prefactor frozen!
        if np.any(source.model.free[1:]) and roi.TS(which=source)< min_ts:
            try:
                catalog_source = catalog.get_source(name)
            except StopIteration:
                pass
            else:
                print 'Freezing spectra of %s to 2FGL prediction b/c it is insignificant' % name
                roi.modify(which=name, model=catalog_source.model, keep_old_flux=False)
                fit_only_prefactor(roi, name)
                any_changed=True
    return any_changed

def freeze_bad_index_to_catalog(roi,catalog,exclude_names=[], min_ts=25):
    """ Fix the spectrum of all power-law catalog sources with a bad spectral
        index to the predictions from the catalog
        with the predictions of 2FGL and 
        the spectral shape of the source frozen. """
    any_changed=False
    for source in roi.get_sources():
        name = source.name

        if name in exclude_names: continue

        if isinstance(source.model,PowerLaw):
            index =source.model['index']
            if index < -5 or index > 5:
                try:
                    catalog_source = catalog.get_source(name)
                except StopIteration:
                    pass
                else:
                    print 'Freezing spectra of %s to 2FGL prediction b/c fit index is bad' % name
                    roi.modify(which=name, model = catalog_source.model, keep_old_flux=False)
                    fit_only_prefactor(roi, name)
                    any_changed=True
    return any_changed



def force_gradient(use_gradient):
    """ A kludge to force use_gradient everywhere! """
    from uw.like.roi_analysis import ROIAnalysis
    from lande.utilities.decorators import modify_defaults
    ROIAnalysis.fit=modify_defaults(use_gradient=use_gradient)(ROIAnalysis.fit)

def galstr(skydir):
    return 'SkyDir(%.3f,%.3f,SkyDir.GALACTIC)' % (skydir.l(),skydir.b())


def spectrum_to_dict(*args, **kwargs):
    return gtlike_or_pointlike(gtlike_spectrum_to_dict, pointlike_spectrum_to_dict, *args, **kwargs)

def name_to_dict(*args, **kwargs):
    return gtlike_or_pointlike(gtlike_name_to_dict, pointlike_name_to_dict, *args, **kwargs)

def get_full_energy_range(*args, **kwargs):
    return gtlike_or_pointlike(gtlike_get_full_energy_range, pointlike_get_full_energy_range, *args, **kwargs)

def fluxdict(*args, **kwargs):
    return gtlike_or_pointlike(gtlike_fluxdict, pointlike_fluxdict, *args, **kwargs)

def sourcedict(*args, **kwargs):
    return gtlike_or_pointlike(gtlike_sourcedict, pointlike_sourcedict, *args, **kwargs)

def powerlaw_upper_limit(*args, **kwargs):
    return gtlike_or_pointlike(gtlike_powerlaw_upper_limit, pointlike_powerlaw_upper_limit, *args, **kwargs)

def upper_limit(*args, **kwargs):
    return gtlike_or_pointlike(gtlike_upper_limit, pointlike_upper_limit, *args, **kwargs)

def get_all_names(*args, **kwargs):
    return gtlike_or_pointlike(gtlike_get_all_names, pointlike_get_all_names, *args, **kwargs)

def get_spatial_model_name(*args, **kwargs):
    return gtlike_or_pointlike(gtlike_get_spatial_model_name, pointlike_get_spatial_model_name, *args, **kwargs)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
