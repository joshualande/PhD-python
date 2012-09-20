import traceback
import sys

import pylab as P
import numpy as np
import pyfits

from skymaps import DiffuseFunction,IsotropicSpectrum,IsotropicPowerLaw,IsotropicConstant
from pyLikelihood import ParameterVector, SpatialMap_cast, PointSource_cast
import pyLikelihood

from uw.like.Models import PowerLaw, FileFunction, CompositeModel
from uw.like.pointspec_helpers import PointSource
from uw.like.roi_extended import ExtendedSource
from uw.like.roi_diffuse import DiffuseSource

from FluxDensity import FluxDensity

from lande.pysed import units
from lande.utilities.tools import tolist
from . tools import gtlike_or_pointlike
from . specplot import SpectrumPlotter

def gtlike_get_full_energy_range(like): return like.energies[[0,-1]]
def pointlike_get_full_energy_range(roi): return roi.bin_edges[[0,-1]]
def get_full_energy_range(*args, **kwargs):
    return gtlike_or_pointlike(gtlike_get_full_energy_range, pointlike_get_full_energy_range, *args, **kwargs)


def pointlike_spectrum_to_dict(model, errors=False):
    """ Package of a spectral model into a handy
        python dictionary.

            >>> m=PowerLaw(norm=1, index=-.5)
            >>> d=spectrum_to_dict(m)
            >>> print d['Norm']
            1.0
            >>> print d['Index']
            -0.5
            >>> d=spectrum_to_dict(m, errors=False)
            >>> d.has_key('Index_err')
            False

        Note, the way to save a ComositeModel is a little different.

            >>> from uw.like.Models import SumModel,LogParabola
            >>> pl=PowerLaw()
            >>> lp=LogParabola()
            >>> c = SumModel(pl,lp)
            >>> s=spectrum_to_dict(c)
            >>> s.keys()
            ['spectrum', 'method', 'name']
            >>> print s['method']
            pointlike
            >>> len(s['spectrum'])
            2
            >>> s['spectrum'][0] == spectrum_to_dict(pl)
            True
            >>> s['spectrum'][1] == spectrum_to_dict(lp)
            True
    """
    d = dict(name = model.name, method='pointlike')
    if isinstance(model,CompositeModel):
        d['spectrum'] = map(pointlike_spectrum_to_dict,model.models)
        return tolist(d)
    else:
        for p in model.param_names:
            d[p] = model[p]
            if errors:
                d['%s_err' % p] = model.error(p)
        for p in model.default_extra_params.keys():
            d[p] = getattr(model,p)
        if d['name'] == 'FileFunction': d['file'] = model.file

        return tolist(d)

def gtlike_spectrum_to_dict(spectrum, errors=False):
    """ Convert a pyLikelihood object to a python 
        dictionary which can be easily saved to a file. """
    parameters=ParameterVector()
    spectrum.getParams(parameters)
    d = dict(name = spectrum.genericName(), method='gtlike')
    for p in parameters: 
        d[p.getName()]= p.getTrueValue()
        if errors: 
            d['%s_err' % p.getName()]= p.error()*p.getScale() if p.isFree() else np.nan
        if d['name'] == 'FileFunction': 
            ff=pyLikelihood.FileFunction_cast(spectrum)
            d['file']=ff.filename()
    return d


def gtlike_name_to_spectral_dict(like, name, errors=False, minos_errors=False, covariance_matrix=False):
    source = like.logLike.getSource(name)
    spectrum = source.spectrum()
    d=gtlike_spectrum_to_dict(spectrum, errors)
    if minos_errors:
        parameters=ParameterVector()
        spectrum.getParams(parameters)
        for p in parameters: 
            pname = p.getName()
            if p.isFree():
                lower,upper=like.minosError(name, pname)
                try:
                    d['%s_lower_err' % pname] = -1*lower*p.getScale()
                    d['%s_upper_err' % pname] = upper*p.getScale()
                except Exception, ex:
                    print 'ERROR computing Minos errors on parameter %s for source %s:' % (pname,name), ex
                    traceback.print_exc(file=sys.stdout)
                    d['%s_lower_err' % pname] = np.nan
                    d['%s_upper_err' % pname] = np.nan
            else:
                d['%s_lower_err' % pname] = np.nan
                d['%s_upper_err' % pname] = np.nan
    if covariance_matrix:
        d['covariance_matrix'] = get_covariance_matrix(like, name)
    return d

def pointlike_name_to_spectral_dict(roi, name, errors=False, covariance_matrix=False):
    model = roi.get_model(name)
    d = pointlike_spectrum_to_dict(model, errors=errors)
    if covariance_matrix:
        d['covariance_matrix'] = get_covariance_matrix(roi, name)
    return d


def gtlike_flux_dict(like,name, emin=None,emax=None,flux_units='erg', energy_units='MeV',
                     errors=True, include_prefactor=False, prefactor_energy=None):
    """ Note, emin, emax, and prefactor_energy must be in MeV """

    if emin is None and emax is None: 
        emin, emax = get_full_energy_range(like)

    cef=lambda e: units.convert(e,'MeV',flux_units)
    ce=lambda e: units.convert(e,'MeV',energy_units)
    f=dict(flux=like.flux(name,emin=emin,emax=emax),
           flux_units='ph/cm^2/s',
           eflux=cef(like.energyFlux(name,emin=emin,emax=emax)),
           eflux_units='%s/cm^2/s' % flux_units,
           emin=ce(emin),
           emax=ce(emax),
           energy_units=energy_units)

    if errors:
        try:
            # incase the errors were not calculated
            f['flux_err']=like.fluxError(name,emin=emin,emax=emax)
            f['eflux_err']=cef(like.energyFluxError(name,emin=emin,emax=emax))
        except Exception, ex:
            print 'ERROR calculating flux error: ', ex
            traceback.print_exc(file=sys.stdout)
            f['flux_err']=-1
            f['eflux_err']=-1

    if include_prefactor:
        assert prefactor_energy is not None
        source = like.logLike.getSource(name)
        spectrum = source.spectrum()
        cp = lambda e: units.convert(e,'1/MeV','1/%s' % flux_units)
        f['prefactor'] = cp(SpectrumPlotter.get_dnde_mev(spectrum,prefactor_energy))
        f['prefactor_units'] = 'ph/cm^2/s/%s' % flux_units
        f['prefactor_energy'] = ce(prefactor_energy)
    return tolist(f)

def gtlike_powerlaw_prefactor_dict(like, name, flux_units='erg', errors=True, minos_errors=False):
    cp = lambda e: units.convert(e,'1/MeV','1/%s' % flux_units)

    source = like.logLike.getSource(name)
    spectrum = source.spectrum()
    assert spectrum.genericName() == 'PowerLaw'
    pref = spectrum.getParam('Prefactor')
    scale = spectrum.getParam('Scale')

    d=dict()
    d['prefactor'] = cp(pref.getTrueValue())
    if errors:
        d['prefactor_err'] = cp(pref.error()*pref.getScale())
    if minos_errors:
        try:
            lower,upper=like.minosError(name, 'Prefactor')
            d['prefactor_lower_err'] = cp(-1*lower*pref.getScale())
            d['prefactor_upper_err'] = cp(upper*pref.getScale())
        except Exception, ex:
            print 'ERROR computing Minos errors on parameter Prefactor for source %s:' % (name), ex
            d['prefactor_lower_err'] = np.nan
            d['prefactor_upper_err'] = np.nan

    d['prefactor_units'] = 'ph/cm^2/s/%s' % flux_units
    d['prefactor_energy'] = scale.getTrueValue()
    d['prefactor_energy_units'] = 'MeV'
    return d


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
                isinstance(dm,IsotropicPowerLaw) or \
                isinstance(dm,IsotropicConstant):
            return 'ConstantValue'
    else:
        raise Exception("Unknown spatial model for source %s" % source.name)

def get_spatial_model_name(*args, **kwargs):
    return gtlike_or_pointlike(gtlike_get_spatial_model_name, pointlike_get_spatial_model_name, *args, **kwargs)

def pointlike_get_all_names(roi):
    """ Get a list of the names of all sources in the pointlike ROI. """
    return np.append(roi.psm.names,roi.dsm.names)

def gtlike_get_all_names(like):
    """ Get a list of the names of all sources in the gtlike ROI. """
    return like.sourceNames()

def get_all_names(*args, **kwargs):
    return gtlike_or_pointlike(gtlike_get_all_names, pointlike_get_all_names, *args, **kwargs)


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

def diffuse_dict(like_or_roi):
    """ Save out all diffuse sources. """

    f = dict()
    bgs = get_background(like_or_roi)
    for name in bgs:
        f[name] = name_to_spectral_dict(like_or_roi, name, errors=True)
    return tolist(f)

def gtlike_ts_dict(like, name, verbosity=True):
    return dict(
        reoptimize=like.Ts(name,reoptimize=True, verbosity=verbosity),
        noreoptimize=like.Ts(name,reoptimize=False, verbosity=verbosity)
        )

def gtlike_source_dict(like, name, emin=None, emax=None, 
                       flux_units='erg', energy_units='MeV', 
                       errors=True, minos_errors=False, covariance_matrix=True,
                       save_TS=True, add_diffuse_dict=True,
                       verbosity=True):

    if emin is None and emax is None:
        emin, emax = get_full_energy_range(like)

    d=dict(
        logLikelihood=logLikelihood(like),
    )

    d['energy'] = energy_dict(emin=emin, emax=emax, energy_units=energy_units)
    
    d['spectrum']= name_to_spectral_dict(like, name, errors=errors, 
                                         minos_errors=minos_errors, covariance_matrix=covariance_matrix)

    if save_TS:
        d['TS']=gtlike_ts_dict(like, name, verbosity=verbosity)

    d['flux']=flux_dict(like,name,
                        emin=emin, emax=emax,
                        flux_units=flux_units, energy_units=energy_units, errors=errors)


    if add_diffuse_dict:
        d['diffuse'] = diffuse_dict(like)

    return tolist(d)

def source_dict(*args, **kwargs):
    return gtlike_or_pointlike(gtlike_source_dict, pointlike_source_dict, *args, **kwargs)

def pointlike_model_to_flux(model, emin, emax, flux_units='erg', energy_units='MeV', 
                            errors=True, include_prefactor=False, prefactor_energy=None):

    cef=lambda e: units.convert(e,'MeV',flux_units)
    ce=lambda e: units.convert(e,'MeV',energy_units)
    f=dict()
    if errors:
        f['flux'],f['flux_err']=model.i_flux(emin=emin,emax=emax,error=True)
        ef,ef_err=model.i_flux(emin=emin,emax=emax,e_weight=1,error=True)
        f['eflux'],f['eflux_err']=cef(ef),cef(ef_err)
    else:
        f['flux']=model.i_flux(emin=emin,emax=emax,error=False)
        ef=model.i_flux(emin=emin,emax=emax,e_weight=1,error=False)
        f['eflux']=cef(ef)

    f['flux_units']='ph/cm^2/s'
    f['eflux_units']='%s/cm^2/s' % flux_units
    f['energy_units']=energy_units
    f['emin'],f['emax']=ce(emin),ce(emax)

    if include_prefactor:
        assert prefactor_energy is not None
        cp = lambda e: units.convert(e,'1/MeV','1/%s' % flux_units)
        f['prefactor'] = cp(model(prefactor_energy))
        f['prefactor_units'] = 'ph/cm^2/s/%s' % flux_units
        f['prefactor_energy'] = ce(prefactor_energy)

    return tolist(f)

def pointlike_flux_dict(roi, which, emin=None, emax=None, *args, **kwargs):

    if emin is None and emax is None:
        emin, emax = get_full_energy_range(roi)

    model=roi.get_model(which)
    return tolist(pointlike_model_to_flux(model, emin, emax, *args, **kwargs))

def pointlike_powerlaw_prefactor_dict(roi, which, flux_units='erg', errors=True):
    model=roi.get_model(which)

    assert isinstance(model,PowerLaw)

    cp = lambda e: units.convert(e,'1/MeV','1/%s' % flux_units)
    d = dict()
    d['prefactor'] = cp(model['norm'])
    if errors:
        d['prefactor_err'] = cp(model.error('norm'))
    d['prefactor_units'] = 'ph/cm^2/s/%s' % flux_units
    d['prefactor_energy'] = model.e0
    d['prefactor_energy_units'] = 'MeV'
    return d


def energy_dict(emin, emax, energy_units='MeV'):
    ce=lambda e: units.convert(e,'MeV',energy_units)
    return dict(emin=ce(emin),
                emax=ce(emax),
                emiddle=ce(np.sqrt(emin*emax)),
                energy_units=energy_units)


def skydirdict(skydir):
    return tolist(dict(
        gal = [skydir.l(),skydir.b()],
        equ = [skydir.ra(),skydir.dec()]))

def pointlike_ts_dict(roi, name):
    return roi.TS(name,quick=False)

def pointlike_source_dict(roi, name, emin=None, emax=None, 
                          flux_units='erg', energy_units='MeV',
                          errors=True, covariance_matrix=True,
                          save_TS=True, 
                          add_diffuse_dict=True,
                          verbosity=True):
    d={}

    if emin is None and emax is None:
        emin, emax = get_full_energy_range(roi)

    old_quiet = roi.quiet; roi.quiet=True
    if save_TS:
        d['TS']=pointlike_ts_dict(roi,name)

    roi.quiet = old_quiet

    d['logLikelihood']=logLikelihood(roi)

    d['energy'] = energy_dict(emin=emin, emax=emax, energy_units=energy_units)

    d['flux']=flux_dict(roi, name, 
                        emin=emin, emax=emax, 
                        flux_units=flux_units, energy_units=energy_units, errors=errors)

    d['spectrum']= name_to_spectral_dict(roi, name, errors=errors, covariance_matrix=covariance_matrix)

    # Source position
    source = roi.get_source(name)
    d['position'] = skydirdict(source.skydir)

    if add_diffuse_dict:
        d['diffuse'] = diffuse_dict(roi)

    d['spatial_model'] = spatial_dict(source, roi, errors=errors)

    return tolist(d)

def spatial_dict(source, roi, errors=True):
    f = dict()
    if isinstance(source,ExtendedSource):
        # Extended Source parameters
        spatial_model = source.spatial_model
        for param in spatial_model.param_names:
            f[param]=spatial_model[param]
            if errors:
                f[param + '_err']=spatial_model.error(param)
        f['r68'] = spatial_model.r68()
        f['r99'] = spatial_model.r99()

    if hasattr(source,'localization'):
        f['ellipse'] = source.localization

    else:
        # add elliptical error, if they exist.
        # N.B. If no localization performed, this will return
        # an empty dictionary.
        # N.B. This method will do the wrong thing if you have recently relocalized
        # another source. This is rarely the case.
        f['ellipse'] = roi.get_ellipse()
    return tolist(f)

def ts_dict(*args, **kwargs):
    return gtlike_or_pointlike(gtlike_ts_dict, pointlike_ts_dict, *args, **kwargs)

def spectrum_to_dict(*args, **kwargs):
    """ Test out FileFunction

            >>> pl= PowerLaw()
            >>> from tempfile import NamedTemporaryFile
            >>> tempfile = NamedTemporaryFile()
            >>> filename = tempfile.name
            >>> pl.save_profile(filename, 1, 100)
            >>> model = FileFunction(file=filename, set_default_limits=True)
            >>> d = spectrum_to_dict(model)
            >>> d.pop('file') == filename
            True
            >>> np.allclose(d.pop('Normalization'),1)
            True
            >>> print d
            {'name': 'FileFunction', 'method': 'pointlike'}
            >>> from . models import build_gtlike_spectrum
            >>> spectrum = build_gtlike_spectrum(model)
            >>> d = spectrum_to_dict(spectrum)
            >>> d.pop('file') == filename
            True
            >>> np.allclose(d.pop('Normalization'),1)
            True
            >>> print d
            {'method': 'gtlike', 'name': 'FileFunction'}
    """
    return gtlike_or_pointlike(gtlike_spectrum_to_dict, pointlike_spectrum_to_dict, *args, **kwargs)

def name_to_spectral_dict(*args, **kwargs):
    return gtlike_or_pointlike(gtlike_name_to_spectral_dict, pointlike_name_to_spectral_dict, *args, **kwargs)

def flux_dict(*args, **kwargs):
    return gtlike_or_pointlike(gtlike_flux_dict, pointlike_flux_dict, *args, **kwargs)

def powerlaw_prefactor_dict(*args, **kwargs):
    return gtlike_or_pointlike(gtlike_powerlaw_prefactor_dict, pointlike_powerlaw_prefactor_dict, *args, **kwargs)

def pointlike_logLikelihood(roi): return -roi.logLikelihood(roi.parameters())

def gtlike_logLikelihood(like): return like.logLike.value()

def logLikelihood(*args, **kwargs):
    return gtlike_or_pointlike(gtlike_logLikelihood, pointlike_logLikelihood, *args, **kwargs)

def gtlike_get_roi_dir(like):
    dir=like.binnedData.countsMap.refDir()
    return dir
def pointlike_get_roi_dir(roi):
    return roi.roi_dir

def get_roi_dir(*args, **kwargs):
    return gtlike_or_pointlike(gtlike_get_roi_dir, pointlike_get_roi_dir, *args, **kwargs)


def gtlike_get_skydir(like, name):
    """ Get the skydir for a gtlike point or extended source. """
    if not name in get_sources(like):
        raise Exception("Unable to get skydir because %s is not a point or extended source." % name)

    spatial_model = get_spatial_model_name(like, name)
    assert spatial_model in ['SkyDirFunction', 'SpatialMap']

    source = like.logLike.getSource(name)

    if spatial_model == 'SkyDirFunction':
        return PointSource_cast(source).getDir()

    elif spatial_model == 'SpatialMap':
        spatial_map=SpatialMap_cast(source.getSrcFuncs()['SpatialDist'])
        
        # This is kind of ugly, but I can't find out how to get m_refDir out of WcsMap2 object
        filename = spatial_map.fitsFile()
        pf=pyfits.open(filename)
        h=pf['PRIMARY'].header
        crpix1=h['CRPIX1']
        crpix2=h['CRPIX2']
        return spatial_map.wcsmap().skyDir(crpix1,crpix2)

def pointlike_get_skydir(roi, name):
    source = roi.get_source(name)
    if not isinstance(source,PointSource) and not isinstance(source,ExtendedSource):
        raise Exception("Unable to get skydir because %s is not a point or extended source." % name)
    return source.skydir

def get_skydir(*args, **kwargs):
    return gtlike_or_pointlike(gtlike_get_skydir, pointlike_get_skydir, *args, **kwargs)


def gtlike_get_covariance_matrix(like, name):
    """ Get the covarince matrix. 

        We can mostly get this from FluxDensity, but
        the covariance matrix returned by FluxDensity
        is only for the free paramters. Here, we
        transform it to have the covariance matrix
        for all parameters, and set the covariance to 0
        when the parameter is free.
    """

    source = like.logLike.getSource(name)
    spectrum = source.spectrum()

    parameters=ParameterVector()
    spectrum.getParams(parameters)
    free = np.asarray([p.isFree() for p in parameters])
    scales = np.asarray([p.getScale() for p in parameters])
    scales_transpose = scales.reshape((scales.shape[0],1))

    cov_matrix = np.zeros([len(parameters),len(parameters)])

    try:
        fd = FluxDensity(like,name)
        cov_matrix[np.ix_(free,free)] = fd.covar

        # create absolute covariance matrix:
        cov_matrix = scales_transpose * cov_matrix * scales
    except Exception, ex:
        print 'ERROR unable to obtain covariance matrix for source %s:' % name, ex
        traceback.print_exc(file=sys.stdout)

    return tolist(cov_matrix)

def pointlike_get_covariance_matrix(roi, name):

    model = roi.get_model(name)
    return model.get_cov_matrix()

def get_covariance_matrix(*args, **kwargs):
    return gtlike_or_pointlike(gtlike_get_covariance_matrix, pointlike_get_covariance_matrix, *args, **kwargs)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
