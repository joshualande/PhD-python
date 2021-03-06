
""" This file contains various function which I have found useful. """
import traceback
import sys

import numpy as np

from uw.like.Models import PowerLaw

from lande.utilities.tools import tolist

from . superstate import SuperState
from . save import logLikelihood
from . tools import gtlike_or_pointlike
from . printing import summary
from . modify import modify

class FitterException(Exception): pass


def paranoid_gtlike_fit(like, covar=True, niter=1, verbosity=False):
    """ Perform a sepctral fit in gtlike in
        a paranoid manner. 
        
        See here for description of method:
            http://fermi.gsfc.nasa.gov/ssc/data/analysis/documentation/Cicerone/Cicerone_Likelihood/Fitting_Models.html
    """
    if niter > 1:
        if verbosity: print 'Fitting %s times' % niter
        for i in range(niter):
            if verbosity: print "Fitting iteration %s" % i
            paranoid_gtlike_fit(like, covar, niter=1, verbosity=verbosity)
        return

    optverbosity = max(verbosity-1, 0) # see IntegralUpperLimit.py

    saved_state = SuperState(like)
    try:
        if verbosity: print 'First, fitting with minuit'
        init_likelihood = logLikelihood(like)
        like.fit(optverbosity, optimizer="MINUIT",covar=covar)
        final_likelihood = logLikelihood(like)
        if init_likelihood - final_likelihood > 10:
            raise FitterException("Error, the final likelihood=%.1d is much worse than the initial likelihood=%.1d (dLL=%.1d)." % (final_likelihood,init_likelihood,final_likelihood-init_likelihood))

    except Exception, ex:
        if verbosity: 
            print 'Minuit fit failed with optimizer=MINUIT, Try again with DRMNFB + NEWMINUIT!', ex
            print 'The bad ROI is:'
            print summary(like)
        traceback.print_exc(file=sys.stdout)
        saved_state.restore()

        saved_state = SuperState(like)
        try:
            init_likelihood = logLikelihood(like)
            if verbosity: print 'Resetting the ROI and Refitting, first with DRMNFB (and after with NEWMINUIT)'
            like.fit(optverbosity, optimizer='DRMNFB', covar=False)
            if verbosity: print 'Refitting, second with NEWMINUIT'
            like.fit(optverbosity, optimizer='NEWMINUIT', covar=covar)
            final_likelihood = logLikelihood(like)
            if init_likelihood - final_likelihood > 10:
                raise FitterException("Error, the final likelihood=%.1d is much worse than the initial likelihood=%.1d (dLL=%.1d)." % (final_likelihood,init_likelihood,final_likelihood-init_likelihood))

        except Exception, ex:
            traceback.print_exc(file=sys.stdout)
            if verbosity: 
                print 'ERROR spectral fitting with DRMNFB + NEWMINUIT: ', ex
                print 'The bad ROI is:'
                print summary(like)

            saved_state.restore()
            try:
                saved_state = SuperState(like)
                if verbosity: print 'Restting the ROI and Refitting with LBFGS'
                like.fit(optverbosity, optimizer='LBFGS', covar=False)
            except Exception, ex:
                print 'ERROR spectral fitting with LBFGS', ex
                traceback.print_exc(file=sys.stdout)
                saved_state.restore()


def gtlike_allow_fit_only_prefactor(like, name):
    """ Freeze everything but norm of source with name
        in pyLikelihood object. """
    modify(like, name, free=False)
    par = like.normPar(name)
    par.setFree(True)
    like.syncSrcParams(name)

def pointlike_allow_fit_only_prefactor(roi, which):
    model = roi.get_model(which)
    old_free = model.free
    new_free = np.zeros_like(old_free).astype(bool)
    new_free[0] = True
    roi.modify(which=which, free=new_free)

def pointlike_fit_only_source(roi, which, fit_only_prefactor=False, **kwargs):
    """ Fit the prefactor of source 'which'
        without varying any other parmters.
        
        Can help if one source has a very bad 
        starting value. """
    source = roi.get_source(which)
    model = roi.get_model(which)
    name = source.name

    frozen_sources = dict()
    for other_source in roi.psm.point_sources.tolist() + roi.dsm.diffuse_sources.tolist():
        other_model = roi.get_model(other_source)
        if np.any(other_model.free) and other_source.name != name:
            frozen_sources[other_source.name]=other_model.free.copy()
            roi.modify(which=other_source,free=False)

    if fit_only_prefactor:
        old_free = roi.get_model(which).free.copy()
        allow_fit_only_prefactor(roi, which)

    roi.fit(**kwargs)

    if fit_only_prefactor:
        roi.modify(which=which, free=old_free)

    for other_name,other_free in frozen_sources.items():
        roi.modify(which=other_name,free=other_free)

# for now, no gtlike implementation
fit_prefactor = lambda *args, **kwargs: pointlike_fit_only_source(*args, fit_only_prefactor=True, **kwargs)
fit_only_source = pointlike_fit_only_source

def allow_fit_only_prefactor(*args, **kwargs):
    return gtlike_or_pointlike(gtlike_allow_fit_only_prefactor, pointlike_allow_fit_only_prefactor, *args, **kwargs)



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
                allow_fit_only_prefactor(roi, name)
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
                    allow_fit_only_prefactor(roi, name)
                    any_changed=True
    return any_changed


def gtlike_setp(like, name, parname, value, scale, lower, upper, free):
    """ Set a paramter in gtlike. """
    par = like[like.par_index(name, parname)]
    par.setBounds(-1e100,1e100)
    par.setScale(scale)
    par.setTrueValue(value)
    par.setBounds(lower,upper)
    par.setFree(free)
    like.syncSrcParams(name)


