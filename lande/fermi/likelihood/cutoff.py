import traceback
import sys

import pylab as P
import numpy as np

from uw.like.Models import Model,PowerLaw,PLSuperExpCutoff
from uw.like.roi_state import PointlikeState

from lande.utilities.tools import tolist

from lande.pysed import units
from lande.fermi.spectra.sed import SED

from . superstate import SuperState

from . tools import gtlike_or_pointlike
from . save import get_full_energy_range, spectrum_to_dict, fluxdict
from . fit import paranoid_gtlike_fit
from . printing import summary
from . models import build_gtlike_model


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
            sed_dict: created by GtlikeSED.todict(). Can also be a yaml
              file created by GtlikeSED.save().

            model_0_kwargs: kwargs for model_0's plot 
            model_1_kwargs: kwargs for model_0's plot 
            sed_kwargs: kwargs to pass into SED
              E.G. flux_units, flux_units, figsize, ...
    """
    sed=SED(sed_results,**sed_kwargs)

    axes=sed.plot(plot_spectral_fit=False, **plot_kwargs)
    axes.autoscale(enable=False)

    model_0 = SED.dict_to_spectrum(cutoff_results['model_0'])
    model_1 = SED.dict_to_spectrum(cutoff_results['model_1'])
    sed.plot_spectrum(model_0, axes=axes, **model_0_kwargs)
    sed.plot_spectrum(model_1, axes=axes, **model_1_kwargs)

    if title is None:
        axes.set_title('Gtlike Cutoff test for %s' % sed.name)
    else:
        axes.set_title(title)


    if filename is not None: P.savefig(filename)
    return axes

def pointlike_test_cutoff(roi, which, model0=None, model1=None, flux_units='erg'):
    print 'Testing cutoff in pointlike'
    emin,emax=get_full_energy_range(roi)
    d = {}

    saved_state = PointlikeState(roi)

    old_flux = roi.get_model(which).i_flux(emin,emax)

    if model0 is not None:
        pass
    elif isinstance(roi.get_model(which),PowerLaw):
        model0 = roi.get_model(which)
    else:
        model0=PowerLaw(norm=1e-11, index=2, e0=np.sqrt(emin*emax))
        model0.set_mapper('Index', PowerLaw.default_limits['Index'])
        model0.set_flux(old_flux,emin=emin,emax=emax)

    print "model0 is ",model0

    roi.modify(which=which, model=model0, keep_old_flux=False)

    fit = lambda: roi.fit(estimate_errors=False)
    ll = lambda: -1*roi.logLikelihood(roi.parameters())
    def ts():
        old_quiet = roi.quiet; roi.quiet=True
        ts = roi.TS(which,quick=False)
        roi.quiet = old_quiet
        return ts

    spectrum = lambda: spectrum_to_dict(roi.get_model(which), errors=True)

    print 'About to fit Model0'
    roi.print_summary()

    fit()
    
    print 'Done fitting Model0'
    roi.print_summary()

    d['ll_0'] = ll_0 = ll()
    d['TS_0'] = ts()
    d['model_0']=spectrum()
    d['flux_0']=fluxdict(roi,which,emin,emax,flux_units)

    if model1 is not None:
        pass
    else:
        model1=PLSuperExpCutoff(norm=1e-9, index=1, cutoff=1000, e0=1000, b=1)
        for p in ['Index', 'Cutoff', 'b']:
            model1.set_mapper(p, PLSuperExpCutoff.default_limits[p])
        model1.set_free('b', False)
        model1.set_flux(old_flux,emin=emin,emax=emax)

    print "model1 is ",model1

    roi.modify(which=which, model=model1, keep_old_flux=False)

    print 'About to fit Model1'
    roi.print_summary()

    fit()

    print 'Done fitting Model1'
    roi.print_summary()

    d['ll_1'] = ll_1 = ll()
    d['TS_1'] = ts()
    d['model_1']=spectrum()
    d['flux_1']=fluxdict(roi,which,emin,emax,flux_units)

    d['TS_cutoff']=2*(ll_1-ll_0)

    saved_state.restore()

    return tolist(d)

def gtlike_test_cutoff(like, name, model0=None, model1=None, flux_units='erg'):
    """ model0 and model1 must be pointlike model objects. """
    print 'Testing cutoff in gtlike'

    saved_state = SuperState(like)

    d = {}

    try:
        emin, emax = get_full_energy_range(like)

        def get_flux():
            return like.flux(name, emin, emax)

        ll = lambda: like.logLike.value()
        ts = lambda: like.Ts(name,reoptimize=True, verbosity=4)

        def spectrum():
            source = like.logLike.getSource(name)
            s=source.spectrum()
            return spectrum_to_dict(s, errors=True)

        old_flux = get_flux()
        if model0 is None:
            model0=PowerLaw(norm=1e-11, index=2, e0=np.sqrt(emin*emax), set_default_limits=True)
            model0.set_flux(old_flux,emin=emin,emax=emax)

        print 'model0 is',model0

        spectrum0=build_gtlike_model(model0)
        like.setSpectrum(name,spectrum0)

        print 'About to fit spectrum0'
        print summary(like)

        paranoid_gtlike_fit(like)

        print 'Done fitting spectrum0'
        print summary(like)

        d['ll_0'] = ll_0 = ll()
        d['TS_0'] = ts()
        d['model_0']=spectrum()
        d['flux_0']=fluxdict(like,name,emin,emax,flux_units)

        if model1 is None:
            model1=PLSuperExpCutoff(norm=1e-9, index=1, cutoff=1000, e0=1000, b=1, set_default_limits=True)
            model1.set_free('b', False)
            model1.set_flux(old_flux,emin=emin,emax=emax)

        print 'model1 is',model1

        spectrum1=build_gtlike_model(model1)
        like.setSpectrum(name,spectrum1)

        print 'About to fit model1'
        print summary(like)

        paranoid_gtlike_fit(like)

        if ll() < ll_0:
            # if fit is worse than PowerLaw fit, then
            # restart fit with parameters almost
            # equal to best fit powerlaw
            cutoff_plaw=PLSuperExpCutoff(b=1, set_default_limits=True)
            cutoff_plaw.set_free('b', False)
            cutoff_plaw.setp_gtlike('norm', d['model_0']['Prefactor'])
            cutoff_plaw.setp_gtlike('index', d['model_0']['Index'])
            cutoff_plaw.setp_gtlike('e0', d['model_0']['Scale'])
            cutoff_plaw.setp_gtlike('cutoff', 1e6)

            temp=build_gtlike_model(cutoff_plaw)
            like.setSpectrum(name,temp)

            print 'Redoing fit with cutoff same as plaw'
            print summary(like)

            paranoid_gtlike_fit(like)


        print 'Done fitting spectrum1'
        print summary(like)

        d['ll_1'] = ll_1 = ll()
        d['TS_1'] = ts()
        d['model_1']=spectrum()
        d['flux_1']=fluxdict(like,name,emin,emax,flux_units)

        d['TS_cutoff']=2*(ll_1-ll_0)

        print 'For cutoff test, TS_cutoff = ', d['TS_cutoff']
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
