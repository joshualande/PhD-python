import traceback
import sys

import numpy as np

from uw.like.Models import PowerLaw, PLSuperExpCutoff
from uw.like.roi_state import PointlikeState
from uw.like.roi_upper_limits import FluxUpperLimit

from lande.utilities.tools import tolist

from lande.pysed import units

from . superstate import SuperState
from . tools import gtlike_or_pointlike
from . save import get_full_energy_range, spectrum_to_dict, pointlike_model_to_flux, gtlike_fluxdict
from . fit import gtlike_allow_fit_only_prefactor, paranoid_gtlike_fit
from . models import build_gtlike_spectrum

def gtlike_upper_limit(like, name, cl=.95, emin=None, emax=None, 
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
        gtlike_allow_fit_only_prefactor(like, name)

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

        source=like.logLike.getSource(name)
        spectrum=source.spectrum()
        prefactor=spectrum.normPar()
        pref_ul = results['ul_value']*prefactor.getScale()
        prefactor.setTrueValue(pref_ul)

        ul = gtlike_fluxdict(like, name, emin=emin,emax=emax,flux_units=flux_units, error=False)

        ul['spectrum'] = spectrum_to_dict(spectrum)

    except Exception, ex:
        print 'ERROR gtlike upper limit: ', ex
        traceback.print_exc(file=sys.stdout)
        ul = None
    finally:
        saved_state.restore()

    return tolist(ul)


def gtlike_powerlaw_upper_limit(like, name, powerlaw_index=2 , cl=0.95, emin=None, emax=None, 
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
    dnde = PowerLaw(norm=1e-11, index=2,e0=1e3)

    like.setSpectrum(name,'PowerLaw')
    like.syncSrcParams(name)

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

    results = gtlike_upper_limit(like, name, emin=emin, emax=emax, **kwargs)

    saved_state.restore()

    return tolist(results)

def gtlike_cutoff_upper_limit(like, name, Index, Cutoff, b, emin=None, emax=None, **kwargs):
    print 'Calculating gtlike cutoff upper limit'

    saved_state = SuperState(like)

    if emin is None and emax is None: 
        emin, emax = get_full_energy_range(like)

    e = np.sqrt(emin*emax)

    cutoff_model = PLSuperExpCutoff(Index=Index, Cutoff=Cutoff, b=b, set_default_limits=True)
    cutoff_spectrum = build_gtlike_spectrum(cutoff_model)

    like.setSpectrum(name,cutoff_spectrum)
    like.syncSrcParams(name)

    results = gtlike_upper_limit(like, name, emin=emin, emax=emax, **kwargs)

    saved_state.restore()

    return tolist(results)


def pointlike_upper_limit(roi, name, cl=.95, emin=None, emax=None, flux_units='erg', **kwargs):

    # Just to be safe, make sure integral is over VERY large range
    if 'integral_min' not in kwargs: kwargs['integral_min']=-20
    if 'integral_max' not in kwargs: kwargs['integral_max']=-5

    if emin is None and emax is None:
        emin, emax = get_full_energy_range(roi)

    saved_state = PointlikeState(roi)

    try:
        ful = FluxUpperLimit(roi=roi, which=name, confidence=cl, **kwargs)
        model = ful.upper_limit_model

        ul = pointlike_model_to_flux(model, emin=emin, emax=emax, flux_units=flux_units, error=False)

        ul['model'] = spectrum_to_dict(model)

    except Exception, ex:
        print 'ERROR pointlike upper limit: ', ex
        traceback.print_exc(file=sys.stdout)
        ul = None
    finally:
        saved_state.restore(just_spectra=True)

    return tolist(ul)


def pointlike_powerlaw_upper_limit(roi, name, powerlaw_index=2, emin=None, emax=None, **kwargs):
    print 'Calculating pointlike upper limit'

    saved_state = PointlikeState(roi)

    """ Note keep old flux, because it is important to have
        the spectral model pushed into the upper_limit
        code reasonably close to the best fit flux. This
        is because initial likelihood (ll_0) is used to scale
        the likelihood so it has to be reasonably close to 
        the best value. """
    roi.modify(which=name, model=PowerLaw(index=powerlaw_index), keep_old_flux=True)

    ul = pointlike_upper_limit(roi, name, emin=emin, emax=emax, **kwargs)

    saved_state.restore(just_spectra=True)

    return tolist(ul)

def pointlike_cutoff_upper_limit(roi, name, Index, Cutoff, b, emin=None, emax=None, **kwargs):
    print 'Calculating pointlike upper limit'

    saved_state = PointlikeState(roi)

    cutoff_model = PLSuperExpCutoff(Index=Index, Cutoff=Cutoff, b=b)
    roi.modify(which=name, model=cutoff_model, keep_old_flux=True)

    ul = pointlike_upper_limit(roi, name, emin=emin, emax=emax, **kwargs)

    saved_state.restore(just_spectra=True)

    return tolist(ul)


def powerlaw_upper_limit(*args, **kwargs):
    return gtlike_or_pointlike(gtlike_powerlaw_upper_limit, pointlike_powerlaw_upper_limit, *args, **kwargs)

def cutoff_upper_limit(*args, **kwargs):
    return gtlike_or_pointlike(gtlike_cutoff_upper_limit, pointlike_cutoff_upper_limit, *args, **kwargs)

def upper_limit(*args, **kwargs):
    return gtlike_or_pointlike(gtlike_upper_limit, pointlike_upper_limit, *args, **kwargs)
