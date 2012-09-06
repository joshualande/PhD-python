import traceback
import sys

import pylab as P
import numpy as np

from uw.utilities import keyword_options
from uw.like.Models import Model,PowerLaw,PLSuperExpCutoff
from uw.like.roi_state import PointlikeState

from lande.utilities.tools import tolist

from lande.pysed import units
from lande.fermi.spectra.sed import SED

from . superstate import SuperState

from . tools import gtlike_or_pointlike
from . save import get_full_energy_range, spectrum_to_dict, energy_dict, source_dict
from . fit import paranoid_gtlike_fit
from . printing import summary
from . models import build_gtlike_spectrum
from . basefit import BaseFitter
from . specplot import SpectrumPlotter,SpectralAxes

class CutoffTester(BaseFitter):

    defaults = BaseFitter.defaults + (
        ('energy_units', 'MeV', 'default units to plot energy flux (y axis) in.'),
        ('flux_units',  'erg', 'default units to plot energy (x axis) in'),
    )

    def plot(self, filename=None, axes=None, title=None, 
             fignum=None, figsize=(4,4),
             model_0_kwargs=dict(color='red', zorder=0),
             model_1_kwargs=dict(color='blue', zorder=0),
             sed_kwargs=dict(),
             plot_kwargs=dict(),):
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
        if axes is None:
            fig = P.figure(fignum,figsize)
            axes = SpectralAxes(fig=fig, rect=(0.22,0.15,0.75,0.8),
                         flux_units=self.flux_units,
                         energy_units=self.energy_units)
            fig.add_axes(axes)
            energy = self.results['energy']
            axes.set_xlim_units(energy['emin']*units.fromstring(energy['energy_units']), 
                                energy['emax']*units.fromstring(energy['energy_units']))

        sp = SpectrumPlotter(axes=axes)
        sp.plot(self.results['hypothesis_0']['spectrum'], **model_0_kwargs)
        sp.plot(self.results['hypothesis_1']['spectrum'], **model_1_kwargs)

        if title is not None: axes.set_title(title)
        if filename is not None: P.savefig(filename)
        return axes

class PointlikeCutoffTester(CutoffTester):

    defaults = CutoffTester.defaults + (
        ('model1',None,'starting value of spectral model'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, roi, name, *args, **kwargs):
        keyword_options.process(self, kwargs)

        self.roi = roi
        self.name = name
        self._calculate()

    def _calculate(self):
        roi = self.roi
        name = self.name
        
        if self.verbosity: print 'Testing cutoff in pointlike'
        emin,emax=get_full_energy_range(roi)

        self.results = d = dict(
            energy = energy_dict(emin=emin, emax=emax, energy_units=self.energy_units)
        )

        saved_state = PointlikeState(roi)

        old_flux = roi.get_model(name).i_flux(emin,emax)

        if not isinstance(roi.get_model(name),PowerLaw):

            model0=PowerLaw(norm=1e-11, index=2, e0=np.sqrt(emin*emax))
            model0.set_mapper('Index', PowerLaw.default_limits['Index'])
            model0.set_flux(old_flux,emin=emin,emax=emax)

            if self.verbosity: print "model0 is ",model0

            roi.modify(which=name, model=model0, keep_old_flux=False)

        fit = lambda: roi.fit(estimate_errors=False)
        ll = lambda: -1*roi.logLikelihood(roi.parameters())
        def ts():
            old_quiet = roi.quiet; roi.quiet=True
            ts = roi.TS(name,quick=False)
            roi.quiet = old_quiet
            return ts

        spectrum = lambda: spectrum_to_dict(roi.get_model(name), errors=True)

        if self.verbosity: 
            print 'About to fit Model0'
            roi.print_summary()

        fit()
        
        if self.verbosity:
            print 'Done fitting Model0'
            roi.print_summary()

        d['hypothesis_0'] = source_dict(roi, name, emin=emin, emax=emax,
                                        flux_units=self.flux_units,
                                        energy_units=self.energy_units,
                                        verbosity=self.verbosity)

        if self.model1 is not None:
            pass
        else:
            self.model1=PLSuperExpCutoff(norm=1e-9, index=1, cutoff=1000, e0=1000, b=1)
            for p in ['Index', 'Cutoff', 'b']:
                self.model1.set_mapper(p, PLSuperExpCutoff.default_limits[p])
            self.model1.set_free('b', False)
            self.model1.set_flux(old_flux,emin=emin,emax=emax)

        if self.verbosity: print "model1 is ",self.model1

        roi.modify(which=name, model=self.model1, keep_old_flux=False)

        if self.verbosity: 
            print 'About to fit Model1'
            roi.print_summary()

        fit()

        if self.verbosity:
            print 'Done fitting Model1'
            roi.print_summary()

        d['hypothesis_1'] = source_dict(roi, name, emin=emin, emax=emax,
                                        flux_units=self.flux_units,
                                        energy_units=self.energy_units,
                                        verbosity=self.verbosity)


        d['TS_cutoff']=2*(d['hypothesis_1']['TS']-d['hypothesis_0']['TS'])

        saved_state.restore()

class GtlikeCutoffTester(CutoffTester):

    defaults = CutoffTester.defaults + (
        ('model1',None,'starting value of spectral model. Must be pointlike model objects.'),
    )


    @keyword_options.decorate(defaults)
    def __init__(self, like, name, *args, **kwargs):
        keyword_options.process(self, kwargs)

        self.like = like
        self.name = name
        self._calculate()

    def _calculate(self):
        like = self.like
        name = self.name

        if self.verbosity: print 'Testing cutoff in gtlike'

        saved_state = SuperState(like)

        emin, emax = get_full_energy_range(like)

        self.results = d = dict(
            energy = energy_dict(emin=emin, emax=emax, energy_units=self.energy_units)
        )

        try:

            def get_flux():
                return like.flux(name, emin, emax)

            def spectrum():
                source = like.logLike.getSource(name)
                s=source.spectrum()
                return spectrum_to_dict(s, errors=True)

            old_flux = get_flux()

            if spectrum()['name'] == 'PowerLaw':
                pass
            else:
                model0=PowerLaw(norm=1e-11, index=2, e0=np.sqrt(emin*emax), set_default_limits=True)
                model0.set_flux(old_flux,emin=emin,emax=emax)

                if self.verbosity: print 'model0 is',model0

                spectrum0=build_gtlike_spectrum(model0)
                like.setSpectrum(name,spectrum0)

            if self.verbosity: 
                print 'About to fit spectrum0'
                print summary(like)

            paranoid_gtlike_fit(like, verbosity=self.verbosity)

            if self.verbosity: 
                print 'Done fitting spectrum0'
                print summary(like)

            d['hypothesis_0'] = source_dict(like, name, emin=emin, emax=emax,
                                            flux_units=self.flux_units,
                                            energy_units=self.energy_units,
                                            verbosity=self.verbosity)

            if self.model1 is None:
                self.model1=PLSuperExpCutoff(norm=1e-9, index=1, cutoff=1000, e0=1000, b=1, set_default_limits=True)
                self.model1.set_free('b', False)
                self.model1.set_flux(old_flux,emin=emin,emax=emax)

            if self.verbosity: 
                print 'model1 is',self.model1

            spectrum1=build_gtlike_spectrum(self.model1)
            like.setSpectrum(name,spectrum1)

            if self.verbosity: 
                print 'About to fit model1'
                print summary(like)

            paranoid_gtlike_fit(like, verbosity=self.verbosity)

            ll = like.logLike.value()

            if ll < d['hypothesis_0']['logLikelihood']:
                # if fit is worse than PowerLaw fit, then
                # restart fit with parameters almost
                # equal to best fit powerlaw
                cutoff_plaw=PLSuperExpCutoff(b=1)
                cutoff_plaw.set_free('b', False)
                cutoff_plaw.setp_gtlike('norm', d['hypothesis_0']['spectrum']['Prefactor'])
                cutoff_plaw.setp_gtlike('index', d['hypothesis_0']['spectrum']['Index'])
                cutoff_plaw.setp_gtlike('e0', d['hypothesis_0']['spectrum']['Scale'])
                cutoff_plaw.setp_gtlike('cutoff', 1e6)
                cutoff_plaw.set_default_limits(oomp_limits=True)

                temp=build_gtlike_spectrum(cutoff_plaw)
                like.setSpectrum(name,temp)

                if self.verbosity: 
                    print 'Redoing fit with cutoff same as plaw'
                    print summary(like)

                paranoid_gtlike_fit(like, verbosity=self.verbosity)

            if self.verbosity: 
                print 'Done fitting spectrum1'
                print summary(like)

            d['hypothesis_1'] = source_dict(like, name, emin=emin, emax=emax,
                                            flux_units=self.flux_units,
                                            energy_units=self.energy_units,
                                            verbosity=self.verbosity)

            d['TS_cutoff']=2*(d['hypothesis_1']['TS']['reoptimize']-d['hypothesis_0']['TS']['reoptimize'])

            if self.verbosity: 
                print 'For cutoff test, TS_cutoff = ', d['TS_cutoff']
        except Exception, ex:
            print 'ERROR gtlike test cutoff: ', ex
            traceback.print_exc(file=sys.stdout)
            self.results = None
        finally:
            saved_state.restore()


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
