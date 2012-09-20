from os.path import expandvars

import traceback
import sys

import sympy
import pylab as P
import numpy as np

from uw.like.Models import PowerLaw, PLSuperExpCutoff
from uw.like.roi_state import PointlikeState
from uw.like.roi_upper_limits import FluxUpperLimit
from uw.utilities import keyword_options
from uw.utilities.parmap import LimitMapper

from lande.utilities.tools import tolist
from lande.utilities.plotting import plot_points

from lande.pysed import units

from . printing import summary
from . superstate import SuperState
from . tools import gtlike_or_pointlike
from . save import get_full_energy_range, spectrum_to_dict, pointlike_model_to_flux, flux_dict
from . fit import gtlike_allow_fit_only_prefactor, paranoid_gtlike_fit
from . models import build_gtlike_spectrum, build_pointlike_model
from . basefit import BaseFitter
from . specplot import SpectralAxes, SpectrumPlotter

class UpperLimit(BaseFitter):

    defaults = BaseFitter.defaults + (
        ('energy_units', 'MeV', 'default units to plot energy flux (y axis) in.'),
        ('flux_units',  'erg', 'default units to plot energy (x axis) in'),
    )

    def plot(self, filename=None, axes=None, title=None,
             fignum=None, figsize=(4,4),
             spectral_kwargs=dict(color='red',zorder=1.9),
            ):
        """ Plot the upper limit. """
        file_energy_units = units.fromstring(self.results['energy_units'])
        emin = self.results['emin']*file_energy_units
        emax = self.results['emax']*file_energy_units
        if axes is None:
            fig = P.figure(fignum,figsize)
            axes = SpectralAxes(fig=fig, 
                                rect=(0.22,0.15,0.75,0.8),
                                flux_units=self.flux_units,
                                energy_units=self.energy_units)
            fig.add_axes(axes)
            axes.set_xlim_units(emin,emax)

        # plot the spectral model
        spectrum = self.results['spectrum']
        sp=SpectrumPlotter(axes=axes)
        sp.plot(spectrum, emin=emin, emax=emax, **spectral_kwargs)

        # plot in the middle an arrow pointing down
        e_middle = units.tosympy([np.sqrt(self.results['emin']*self.results['emax'])],file_energy_units)
        dnde = sp.get_dnde(spectrum, e_middle)
        energies, e2_dnde = axes.convert_points(e_middle, dnde)
        if 'autoscale' in spectral_kwargs: spectral_kwargs.pop('autoscale')
        if 'label' in spectral_kwargs: spectral_kwargs.pop('label')
        plot_points(x=energies, y=e2_dnde,
                    xlo=None, xhi=None,
                    y_lower_err=None, y_upper_err=None,
                    y_ul=e2_dnde,
                    significant=False,
                    axes=axes, **spectral_kwargs)

        if title is not None: axes.set_title(title)
        if filename is not None: 
            P.savefig(expandvars(filename))
        return axes


class GtlikeUpperLimit(UpperLimit):
    """ Compute gtlike upper limit for whatever spectral model is
        currently in like object. """

    defaults = UpperLimit.defaults + (
        ('cl', 0.95, 'confidence level'),
        ('emin', None, 'minimum energy (for quoted flux). Default is full energy range'),
        ('emax', None, 'maximum energy (for quoted flux). Default is full energy range'),
        ('include_prefactor', False, 'Compute prefactor upper limit'),
        ('prefactor_energy', None, 'Energy to compute prefactor energy at'),
        ('upper_limit_kwargs', dict(), 'Kwargs passed into IntegralUpperLimit.calc_int'),
        ('rescale_parameters_before_limit', True, '...')
    )

    @keyword_options.decorate(defaults)
    def __init__(self, like, name, **kwargs):
        keyword_options.process(self, kwargs)

        self.like = like
        self.name = name

        if self.emin is None and self.emax is None: 
            self.emin, self.emax = get_full_energy_range(like)
            self.e = np.sqrt(self.emin*self.emax)

        self._compute()

    def _compute(self):
        if self.verbosity: print 'Calculating gtlike upper limit'

        like = self.like
        name = self.name

        saved_state = SuperState(like)
        source = like.logLike.getSource(name)

        try:
            import IntegralUpperLimit

            # First, freeze spectrum (except for normalization)
            # of our soruce
            gtlike_allow_fit_only_prefactor(like, name)

            # Spectral fit whole ROI

            if self.verbosity:
                print 'Before fitting before Upper Limit:'
                print summary(like)

            """ N.B. spectral fit in this function instead
                of in upper limits code since my
                paranoid_gtlike_fit function is more robust. """
            paranoid_gtlike_fit(like, verbosity=self.verbosity)

            if self.verbosity:
                print 'After fitting before Upper Limit:'
                print summary(like)

            # Freeze everything but our source of interest
            for i in range(len(like.model.params)):
                like.model[i].setFree(False)
                like.syncSrcParams(like[i].srcName)

            if self.rescale_parameters_before_limit:
                source=like.logLike.getSource(name)
                spectrum=source.spectrum()
                model = build_pointlike_model(spectrum)
                if self.verbosity:
                    print 'Rescaling %s model parameters before limits' % model.name
                    print ' * Initial mappers:',model.mappers
                norm_name = model.param_names[0]
                default_norm_limits = model.default_limits[norm_name]
                lower,upper = model.get_limits(norm_name)
                norm_min = min(lower, default_norm_limits.lower)
                norm_max = max(upper, default_norm_limits.upper)

                model.set_limits(norm_name, norm_min, norm_max, scale=1)
                if self.verbosity:
                    print ' * New mappers:',model.mappers

                spectrum = build_gtlike_spectrum(model)
                like.setSpectrum(name,spectrum)
                like.syncSrcParams(name)

            # Note, I think freeze_all is redundant, but flag it just 
            # to be paranoid
            flux_ul, results = IntegralUpperLimit.calc_int(like, name, 
                                                           freeze_all=True,
                                                           skip_global_opt=True,
                                                           cl=self.cl,
                                                           emin=self.emin, 
                                                           emax=self.emax, 
                                                           verbosity=self.verbosity,
                                                           **self.upper_limit_kwargs)

            if self.verbosity:
                print 'After computing Upper limit:'
                print summary(like)

            source=like.logLike.getSource(name)
            spectrum=source.spectrum()
            prefactor=spectrum.normPar()
            pref_ul = results['ul_value']*prefactor.getScale()
            prefactor.setTrueValue(pref_ul)

            self.results = flux_dict(like, name, 
                                     emin=self.emin,emax=self.emax,
                                     flux_units=self.flux_units, 
                                     energy_units=self.energy_units,
                                     errors=False,
                                     include_prefactor=self.include_prefactor,
                                     prefactor_energy=self.prefactor_energy)

            self.results['spectrum'] = spectrum_to_dict(spectrum)
            self.results['confidence'] = self.cl

        except Exception, ex:
            print 'ERROR gtlike upper limit: ', ex
            traceback.print_exc(file=sys.stdout)
            self.results = None
        finally:
            saved_state.restore()



class GtlikePowerLawUpperLimit(GtlikeUpperLimit):

    defaults = GtlikeUpperLimit.defaults + (
        ('powerlaw_index', 2, 'PowerLaw index for assuemd powerlaw'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, *args, **kwargs):
        super(GtlikePowerLawUpperLimit,self).__init__(*args, **kwargs)

    def _compute(self):
        """ Wrap up calculating the flux upper limit for a powerlaw
            source.  This function employes the pyLikelihood function
            IntegralUpperLimit to calculate a Bayesian upper limit.

            The primary benefit of this function is that it replaces the
            spectral model automatically with a PowerLaw spectral model
            and fixes the index to -2. It then picks a better scale for the
            powerlaw and gives the upper limit calculation a more reasonable
            starting value, which helps the convergence.
        """
        if self.verbosity: print 'Calculating gtlike power-law upper limit'
        like = self.like
        name = self.name

        saved_state = SuperState(like)

        e = np.sqrt(self.emin*self.emax)

        """ I had tons of trouble getting a robust fitting algorithm. 
            The problem with computing upper limits is
                (a) Getting an initial fit of the region (with the spectral index
                    fixed) to converge
                (b) Getting the upper limit to integrate over a good range.

            This is what I found to be most robust way to compute upper limits:
                (a) Create a generic powerlaw model with the spectral index fixed
                    at the desired value (typically set to -2).
                    Note, don't set e0, use default. This ensures that the
                    prefactor range really does convert to a physically reasonable
                    range of parameters.
                (b) Give the spectral model the pointlike default spectral limits.
                    This is important because it gives the source a big enough range 
                    such that the upper limit can find a proper integration range.
                (c) Set the flux of the current model to equal the flux of the input model.
                    This starts the fitter at a reasonable value. Do this setting with 
                    the set_flux flag strict=False beacuse, in case the initial fit totally
                    failed to converge (flux -> 0), you don't want to put the starting
                    value fo the flux too far away from the true value.
                (d) Keep the lower and upper limit on the prefactor as the pointlike default limits,
                    but set the scale of the source to be the new 'norm' found by preserving
                    the flux. This ensures the fitter doesn't have too much trouble finding
                    the true minimum.
            
            Using this procedure, you get a reasonable parameter limits which allows the
            preliminary fit to converge and the upper limits code to integrate over a reasonable
            parameter range.
        """
        source = like.logLike.getSource(name)
        spectrum=source.spectrum()

        old_flux = like.flux(name, self.emin, self.emax)
        model = PowerLaw(index=self.powerlaw_index)
        model.set_flux(old_flux, emin=self.emin, emax=self.emax, strict=False)
        model.set_default_limits(oomp_limits=True)

        spectrum = build_gtlike_spectrum(model)
        like.setSpectrum(name,spectrum)
        like.syncSrcParams(name)

        results = super(GtlikePowerLawUpperLimit,self)._compute()

        saved_state.restore()


class GtlikeCutoffUpperLimit(GtlikeUpperLimit):

    defaults = GtlikeUpperLimit.defaults + (
        ('Index', 1.7, "'Index' parameter of PLSuperExpCutoff"),
        ('Cutoff', 3e3, "'Cutoff' parameter of PLSuperExpCutoff"),
        ('b', 1, "'b'parameter of PLSuperExpCutoff"),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, *args, **kwargs):
        super(GtlikeCutoffUpperLimit,self).__init__(*args, **kwargs)

    def _compute(self):
        if self.verbosity: print 'calculating gtlike cutoff upper limit'
        like = self.like
        name = self.name

        saved_state = SuperState(like)
        
        old_flux = like.flux(name, self.emin, self.emax)
        cutoff_model = PLSuperExpCutoff(Index=self.Index, Cutoff=self.Cutoff, b=self.b)
        cutoff_model.set_flux(old_flux, emin=self.emin, emax=self.emax, strict=False)
        cutoff_model.set_default_limits(oomp_limits=True)

        cutoff_spectrum = build_gtlike_spectrum(cutoff_model)
        like.setSpectrum(name,cutoff_spectrum)
        like.syncSrcParams(name)

        super(GtlikeCutoffUpperLimit,self)._compute()

        saved_state.restore()

class PointlikeUpperLimit(UpperLimit):
    defaults = UpperLimit.defaults + (
        ('cl', 0.95, 'confidence level'),
        ('emin', None, 'minimum energy (for quoted flux). Default is full energy range'),
        ('emax', None, 'maximum energy (for quoted flux). Default is full energy range'),
        ('flux_units', 'erg', 'Units to quote flux in'),
        ('include_prefactor', False, 'Compute prefactor upper limit'),
        ('prefactor_energy', None, 'Energy to compute prefactor energy at'),
        ('simps_points',   100, 'Number of integration points (per decade of energy).'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, roi, name, **kwargs):
        keyword_options.process(self, kwargs)

        self.roi = roi 
        self.name = name

        if self.emin is None and self.emax is None: 
            self.emin, self.emax = get_full_energy_range(roi)


        self._compute()

    def _compute(self):

        roi = self.roi
        name = self.name
        model = roi.get_model(name)

        saved_state = PointlikeState(roi)

        try:
            ful = FluxUpperLimit(roi=roi, 
                                 which=name, 
                                 confidence=self.cl, 
                                 simps_points=self.simps_points,
                                 verbosity=self.verbosity)
            model = ful.upper_limit_model

            self.results  = pointlike_model_to_flux(model, 
                                                    emin=self.emin, emax=self.emax, 
                                                    flux_units=self.flux_units, 
                                                    energy_units=self.energy_units, 
                                                    errors=False,
                                                    include_prefactor=self.include_prefactor,
                                                    prefactor_energy=self.prefactor_energy,
                                                   )
            self.results['confidence'] = self.cl

            self.results['spectrum'] = spectrum_to_dict(model)

        except Exception, ex:
            print 'ERROR pointlike upper limit: ', ex
            traceback.print_exc(file=sys.stdout)
            self.results = None
        finally:
            saved_state.restore(just_spectra=True)



class PointlikePowerLawUpperLimit(PointlikeUpperLimit):

    defaults = PointlikeUpperLimit.defaults + (
        ('powerlaw_index', 2, 'PowerLaw index for assuemd powerlaw'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, *args, **kwargs):
        super(PointlikePowerLawUpperLimit,self).__init__(*args, **kwargs)

    def _compute(self):
        if self.verbosity: 
            print 'Calculating pointlike upper limit'

        roi = self.roi
        name = self.name

        saved_state = PointlikeState(roi)

        """ Note keep old flux, because it is important to have
            the spectral model pushed into the upper_limit
            code reasonably close to the best fit flux. This
            is because initial likelihood (ll_0) is used to scale
            the likelihood so it has to be reasonably close to 
            the best value. """
        model = PowerLaw(index=self.powerlaw_index)
        roi.modify(which=name, model=model, keep_old_flux=True)

        super(PointlikePowerLawUpperLimit,self)._compute()

        saved_state.restore(just_spectra=True)


class PointlikeCutoffUpperLimit(PointlikeUpperLimit):

    defaults = PointlikeUpperLimit.defaults + (
        ('Index', 1.7, "'Index' parameter of PLSuperExpCutoff"),
        ('Cutoff', 3e3, "'Cutoff' parameter of PLSuperExpCutoff"),
        ('b', 1, "'b'parameter of PLSuperExpCutoff"),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, *args, **kwargs):
        super(PointlikeCutoffUpperLimit,self).__init__(*args, **kwargs)

    def _compute(self):
        if self.verbosity: print 'calculating pointlike cutoff upper limit'

        roi = self.roi
        name = self.name

        saved_state = PointlikeState(roi)

        cutoff_model = PLSuperExpCutoff(Index=self.Index, Cutoff=self.Cutoff, b=self.b)
        roi.modify(which=name, model=cutoff_model, keep_old_flux=True)

        super(PointlikeCutoffUpperLimit,self)._compute()

        saved_state.restore(just_spectra=True)

