
""" This file contains various function which I have found useful. """
import traceback
import sys

import numpy as np

from lande.utilities.toolbag import tolist

from . superstate import SuperState

from . tools import gtlike_or_pointlike

def paranoid_gtlike_fit(like, covar=True):
    """ Perform a sepctral fit in gtlike in
        a paranoid manner. 
        
        See here for description of method:
            http://fermi.gsfc.nasa.gov/ssc/data/analysis/documentation/Cicerone/Cicerone_Likelihood/Fitting_Models.html
    """
    saved_state = SuperState(like)
    try:
        print 'First, fitting with minuit'
        like.fit(optimizer="MINUIT",covar=covar)
    except Exception, ex:
        print 'Minuit fit failed with optimizer=MINUIT, Try again with DRMNFB + NEWMINUIT!', ex
        traceback.print_exc(file=sys.stdout)
        saved_state.restore()

        try:
            saved_state = SuperState(like)
            print 'Refitting, first with DRMNFB'
            like.fit(optimizer='DRMNFB', covar=False)
            print 'Refitting, second with NEWMINUIT'
            like.fit(optimizer='NEWMINUIT', covar=covar)
        except Exception, ex:
            print 'ERROR spectral fitting with DRMNFB + NEWMINUIT: ', ex
            traceback.print_exc(file=sys.stdout)
            saved_state.restore()
            try:
                saved_state = SuperState(like)
                print 'Refitting with LBFGS'
                like.fit(optimizer='LBFGS', covar=False)
            except Exception, ex:
                print 'ERROR spectral fitting with LBFGS', ex
                traceback.print_exc(file=sys.stdout)
                saved_state.restore()


def gtlike_fit_only_prefactor(like, name):
    """ Freeze everything but norm of source with name
        in pyLikelihood object. """
    gtlike_modify(like, name, free=False)
    par = like.normPar(name)
    par.setFree(True)
    like.syncSrcParams(name)

def pointlike_fit_only_prefactor(roi, which):
    model = roi.get_model(which)
    old_free = model.free
    new_free = np.zeros_like(old_free).astype(bool)
    new_free[0] = True
    roi.modify(which=which, free=new_free)

def fit_prefactor(roi, which, *args, **kwargs):
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

    old_free = model.free.copy()
    new_free = np.zeros_like(model.free).astype(bool)
    new_free[0] = True
    roi.modify(which=which, free=new_free)

    roi.fit(*args, **kwargs)

    roi.modify(which=which, free=old_free)

    for other_name,other_free in frozen_sources.items():
        roi.modify(which=other_name,free=other_free)

def fit_only_prefactor(*args, **kwargs):
    return gtlike_or_pointlike(gtlike_fit_only_prefactor, pointlike_fit_only_prefactor, *args, **kwargs)

