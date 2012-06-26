""" Contain helper functions for iterating over an ROI
    and doing otherwise repeative tasks. """
import os
import sys

import numpy as np

from skymaps import SkyDir

from uw.like.pointspec_helpers import PointSource
from uw.like.Models import PowerLaw

from . fit import fit_prefactor
from . pretty import pformat

def new_ps(roi, name, l, b, tsmap_kwargs=dict(size=10)):
    """ A 'throw away' convenience function to add a new source to
        the ROI, localize it, and then print out a string which can be
        used to modify an ROI to add a new source. """

    skydir = SkyDir(l,b,SkyDir.GALACTIC)

    roi.print_summary(galactic=True)

    emin,emax=roi.bin_edges[[0,-1]]

    model=PowerLaw(e0=np.sqrt(emin*emax))

    ps = PointSource(
        name=name,
        skydir=skydir,
        model=model)

    roi.add_source(ps)

    roi.print_summary(galactic=True)

    fit_prefactor(roi, name, use_gradient=False)
    roi.fit(use_gradient=False)


    roi.print_summary(galactic=True)

    roi.localize(which=name, update=True)
    roi.fit(use_gradient=False)

    roi.print_summary(galactic=True)

    print roi

    ts = roi.TS(which=name, quick=False)
    print 'TS for source %s is %.1f' % (name,ts)

    path=os.path.abspath(sys.argv[0])
    print """
Code to create point source:

    # Analysis came from %s
    ps=%s
    roi.add_source(ps)
""" % (path,pformat(ps))

    roi.save('roi.dat')

    roi.plot_tsmap(filename='residual_tsmap.pdf', 
                   fitsfile='residual_tsmap.fits',
                   **tsmap_kwargs)


def free_src(roi, name, tsmap_kwargs=dict(size=10)):
    """ Throw away function to free a source in the ROI, refit
        the ROI, save it to a file and make a TS map. """

    roi.print_summary(galactic=True)
    if isinstance(name,basestring):
        roi.modify(which=name, free=True)
    else:
        for n in name: roi.modify(which=n, free=True)

    roi.fit(use_gradient=False)

    roi.print_summary(galactic=True)

    path=os.path.abspath(sys.argv[0])
    print """
Code to free sources:"

    # Analysis came from %s""" % (path)

    def print_name(n):
        print """    roi.modify(which="%s", free=True)""" % n
    if isinstance(name,basestring):
        print_name(name)
    else:
        for n in name: print_name(n)


    roi.save('roi.dat')

    roi.plot_tsmap(filename='residual_tsmap.pdf', 
                   fitsfile='residual_tsmap.fits',
                   **tsmap_kwargs)


def fix_convergence(roi, which, model, tsmap_kwargs=dict(size=10)):
    """ Fix the convergence for a source. """

    roi.modify(which=which, model=model, keep_old_flux=False)

    roi.fit(use_gradient=False)

    roi.print_summary()
    print roi

    path=os.path.abspath(sys.argv[0])
    print """
Code to fix convergence

    # Analysis came from %s
    roi.modify(which='%s', model=%s, keep_old_flux=False)
    """ % (path,which,pformat(model))

    roi.save('roi.dat')

    roi.plot_tsmap(filename='residual_tsmap.pdf', 
                   fitsfile='residual_tsmap.fits',
                   **tsmap_kwargs)

