import traceback
import sys

import pylab as P
import numpy as np

from uw.like.Models import Model,PowerLaw,ExpCutoff
from uw.like.roi_state import PointlikeState

from lande.utilities.toolbag import tolist
from SED import SED

from lande.pysed import units

from . superstate import SuperState

from . tools import gtlike_or_pointlike

def plot_gtlike_cutoff_test(cutoff_results, sed_results, filename=None, title=None, 
                            model_0_kwargs=dict(color='red', zorder=0),
                            model_1_kwargs=dict(color='blue', zorder=0),
                            sed_kwargs=dict(),
                            plot_kwargs=dict(),
                           ):
    """ Plots the cutoff test performed of a spectrum using the function
        gtlike_test_cutoff.

        Input:
            cutoff_dict: created by gtlike_test_cutoff
            sed_dict: created by LandeSED.todict(). Can also be a yaml
              file created by LandeSED.save().

            model_0_kwargs: kwargs for model_0's plot 
            model_1_kwargs: kwargs for model_0's plot 
            sed_kwargs: kwargs to pass into LandeSED
              E.G. flux_units, flux_units, figsize, ...
    """
    sed=LandeSED(sed_results,**sed_kwargs)

    axes=sed.plot(plot_spectral_fit=False, **plot_kwargs)
    axes.autoscale(enable=False)

    model_0 = LandeSED.dict_to_spectrum(cutoff_results['model_0'])
    model_1 = LandeSED.dict_to_spectrum(cutoff_results['model_1'])
    sed.plot_spectrum(model_0, **model_0_kwargs)
    sed.plot_spectrum(model_1, **model_1_kwargs)

    if title is None:
        axes.set_title('Gtlike Cutoff test for %s' % sed.name)
    else:
        axes.set_title(title)


    if filename is not None: P.savefig(filename)
    return axes

def pointlike_test_cutoff(roi, which, flux_units='erg'):
    print 'Testing cutoff in pointlike'
    emin,emax=get_full_energy_range(roi)
    d = {}

    saved_state = PointlikeState(roi)

    print 'these are probably not good startin values!'
    old_flux = roi.get_model(which).i_flux(emin,emax)
    m=PowerLaw(norm=1e-11, index=2, e0=np.sqrt(emin*emax))
    m.set_flux(old_flux,emin,emax)
    roi.modify(which=which, model=m,keep_old_flux=False)

    fit = lambda: roi.fit(estimate_errors=False)
    ll = lambda: -1*roi.logLikelihood(roi.parameters())
    def ts():
        old_quiet = roi.quiet; roi.quiet=True
        ts = roi.TS(which,quick=False)
        roi.quiet = old_quiet
        return ts

    spectrum = lambda: spectrum_to_dict(roi.get_model(which), errors=True)

    fit()
    d['ll_0'] = ll_0 = ll()
    d['TS_0'] = ts()
    d['model_0']=spectrum()
    d['flux_0']=fluxdict(roi,which,emin,emax,flux_units)

    m=ExpCutoff(n0=1e-11, gamma=1, cutoff=1000, e0=1000)
    m.set_flux(old_flux,emin,emax)
    roi.modify(which=which, 
               model=m,keep_old_flux=False)

    fit()
    d['ll_1'] = ll_1 = ll()
    d['TS_1'] = ts()
    d['model_1']=spectrum()
    d['flux_1']=fluxdict(roi,which,emin,emax,flux_units)

    d['TS_cutoff']=2*(ll_1-ll_0)

    saved_state.restore()

    return tolist(d)

def gtlike_test_cutoff(like, name, flux_units='erg'):
    print 'Testing cutoff in gtlike'

    saved_state = SuperState(like)

    d = {}

    try:
        emin, emax = get_full_energy_range(like)

        def fix(parname,value):
            par=like[like.par_index(name, parname)]
            par.setScale(1)
            par.setBounds(-1e100,1e100) # kind of lame, but i think this is necessary
            par.setTrueValue(value)
            par.setBounds(value,value)
            par.setFree(0)
            like.syncSrcParams(name)

        def set(parname,value,scale,lower,upper):
            """ Note, lower + upper are fractional limits if free=True. """
            par=like[like.par_index(name, parname)]
            par.setBounds(-1e100,1e100) # kind of lame, but i think this is necessary
            par.setScale(scale)
            par.setTrueValue(value)
            par.setBounds(lower,upper)
            par.setFree(1)
            like.syncSrcParams(name)

        def get_flux():
            return like.flux(name, emin, emax)

        def set_flux(flux):
            current_flux = get_flux()
            prefactor=like[like.par_index(name, 'Prefactor')]
            prefactor.setTrueValue(
                (flux/current_flux)*prefactor.getTrueValue())
            like.syncSrcParams(name)

        ll = lambda: like.logLike.value()
        ts = lambda: like.Ts(name,reoptimize=True)

        def spectrum():
            source = like.logLike.getSource(name)
            s=source.spectrum()
            return spectrum_to_dict(s, errors=True)

        old_flux = get_flux()

        like.setSpectrum(name,'PowerLaw')
        fix('Scale', np.sqrt(emin*emax))
        set('Prefactor',1e-11,1e-11, 1e-5, 1e5)
        set('Index',       -2,    1,   -5,   5)
        set_flux(old_flux)

        paranoid_gtlike_fit(like)
        d['ll_0'] = ll_0 = ll()
        d['TS_0'] = ts()
        d['model_0']=spectrum()
        d['flux_0']=fluxdict(like,name,emin,emax,flux_units)
        
        like.setSpectrum(name,'PLSuperExpCutoff')
        set('Prefactor', 1e-9,  1e-11,   1e-5,1e5)
        set('Index1',      -1,      1,     -5,  5)
        fix('Scale',     1000)
        set('Cutoff',    1000,      1,    1e2,1e8)
        fix('Index2',       1)
        set_flux(old_flux)

        paranoid_gtlike_fit(like)

        if ll() < ll_0:
            # if fit is worse than PowerLaw fit, then
            # restart fit with parameters almost
            # equal to best fit powerlaw
            m = d['model_0']
            set('Prefactor', m['Prefactor'],  1e-11,   1e-5,1e5)
            set('Index1',        m['Index'],      1,     -5,  5)
            fix('Scale',         m['Scale'])
            set('Cutoff',               1e6,      1,    1e2,1e8)
            fix('Index2',                 1)
            paranoid_gtlike_fit(like)

        d['ll_1'] = ll_1 = ll()
        d['TS_1'] = ts()
        d['model_1']=spectrum()
        d['flux_1']=fluxdict(like,name,emin,emax,flux_units)

        d['TS_cutoff']=2*(ll_1-ll_0)
    except Exception, ex:
        print 'ERROR gtlike test cutoff: ', ex
        traceback.print_exc(file=sys.stdout)
        d=-1
    finally:
        saved_state.restore()

    return tolist(d)

def test_cutoff(*args, **kwargs):
    return gtlike_or_pointlike(gtlike_test_cutoff, pointlike_test_cutoff, *args, **kwargs)

def fix_bad_cutoffs(roi, exclude_names):
    """ Loop over all sources. When ExpCutoff souce has cutoff>10TeV, convert to powerlaw. """
    any_changed=False
    for source in roi.get_sources():
        if source.name in exclude_names: continue

        model = source.model

        if np.any(model.free) and isinstance(model,ExpCutoff) and model['cutoff'] > 1e7:
            print 'Converting cutoff source %s to powerlaw because cutoff too high' % source.name
            new_model = PowerLaw(norm=model['norm'], index=model['index'], e0=model.e0)

            any_changed = True
            roi.modify(which=source, model=new_model, keep_old_flux=False)
    return any_changed
