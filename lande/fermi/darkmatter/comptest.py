""" Code to perform a likelihood-ratio test using 
    a comprehensive model.

    Author: Joshua Lande, Sheridan Zalewski, Alex Drlica-Wagner 
"""
import pylab as P

from uw.like.roi_state import PointlikeState
from uw.darkmatter.spectral import ComprehensiveModel
from uw.utilities import keyword_options

from lande.fermi.likelihood.fit import fit_prefactor 
from lande.fermi.likelihood.save import get_full_energy_range, spectrum_to_dict

from lande.fermi.spectra.pointlike import pointlike_sed_to_dict
from lande.fermi.spectra.sed import SED

from lande.utilities.tools import tolist

class ComprehensiveTest(object):
    """ Perform the likelihood ratio test using a "comprehensive model"

        Usage:
            roi = ...
            test = ComprehensiveTest(roi, 'source_name', PowerLaw(), DMFitFunction(), 
                                     fit_kwargs=dict(use_gradient=False))
            print test.todict()
            test.plot('test.pdf')
    """

    defaults = (
        ('fit_kwargs', dict(), 'Parameters to pass into the roi.fit function'),
        ('localize', False, 'Localize source with both hypothesis'),
        ('quiet', False, 'Shut this code up'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, roi, which, model0, model1, **kwargs):
        keyword_options.process(self, kwargs)

        state = PointlikeState(roi)

        self.roi=roi
        self.which=which
        self.model0=model0.copy()
        self.model1=model1.copy()

        self.compute()
        state.restore()

    def fit(self):
        roi=self.roi
        which=self.which

        fit_prefactor(roi, which, **self.fit_kwargs)
        roi.fit(**self.fit_kwargs)
        if self.localize: 
            roi.localize(which=which, update=True)
            roi.fit(**self.fit_kwargs)
        ll=-roi.logLikelihood(roi.parameters())
        return ll

    def compute(self):
        roi=self.roi
        which=self.which

        roi.modify(which=which, model=self.model0, keep_old_flux=False)

        if not self.quiet: print 'fit source %s with model0 (%s)' % (which,self.model0.name)
        self.ll_0 = self.fit()

        if not self.quiet:
            print roi

        self.sed_points = pointlike_sed_to_dict(roi.plot_sed(which=which))

        # TODO, get SED points

        emin,emax=get_full_energy_range(roi)

        flux=self.model0.i_flux(emin=emin,emax=emax)
        self.model1.set_flux(flux,emin=emin,emax=emax)
        self.comprehensive_model = ComprehensiveModel(self.model0.copy(), self.model1.copy())
        self.comprehensive_model.theta=0.5

        self.comprehensive_model_start = self.comprehensive_model.copy()

        if not self.quiet: print 'fit source %s with model1 (%s)' % (which,self.model1.name)
        roi.modify(which=self.which, model=self.comprehensive_model, keep_old_flux=False)

        self.ll_1 = self.fit()

        if not self.quiet:
            print roi

        self.TS_comp = 2*(self.ll_1-self.ll_0)

    def todict(self):
        return tolist(dict(
            sed_points = self.sed_points,
            model0 = spectrum_to_dict(self.model0),
            model1 = spectrum_to_dict(self.comprehensive_model),
            TS_comp = self.TS_comp,
            ll_1 = self.ll_1,
            ll_0 = self.ll_0,
        ))

    def plot(self, filename):
        sed = SED(self.sed_points)
        # Plot SED points
        axes=sed.plot()
        # Plot model0, model1
        sed.plot_spectrum(self.model0, axes=axes, label='model0')
        sed.plot_spectrum(self.comprehensive_model, axes=axes, label='model1')
        sed.plot_spectrum(self.comprehensive_model_start, axes=axes, label='model1 (start)')
        axes.legend()
        P.savefig(filename)

