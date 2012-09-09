import numpy as np
from matplotlib.axes import Axes

import pyLikelihood
from SED import SED as BaseGtlikeSED

# overload SED namespace for our base object
from uw.utilities import keyword_options
from uw.like.Models import Model

from lande.pysed import units

from . load import dict_to_spectrum
from . models import gtlike_unscale_all_parameters

class SEDException(Exception): 
    pass

class SpectralAxes(Axes):
    def __init__(self, energy_units='MeV', flux_units='erg', *args, **kwargs):
        self.energy_units = energy_units
        self.flux_units = flux_units

        self.energy_units_obj = units.fromstring(self.energy_units)
        self.flux_units_obj = units.fromstring(self.flux_units)


        super(SpectralAxes,self).__init__(*args, **kwargs)
        
        self.set_xscale('log')
        self.set_yscale('log')

        self.set_xlabel('Energy (%s)' % self.energy_units)
        self.set_ylabel('E$^2$ dN/dE (%s cm$^{-2}$ s$^{-1}$)' % self.flux_units)

    def set_xlim_units(self, emin, emax):
        self.set_xlim(float(emin/self.energy_units_obj), float(emax/self.energy_units_obj))

    def set_ylim_units(self, fmin, fmax):
        f = self.flux_units_obj/units.cm**2/units.s
        self.set_ylim(float(fmin/f), float(fmax/f))



class SpectrumPlotter(object):
    """ Plot spectra. """
    def __init__(self, axes):
        self.axes = axes

    @staticmethod
    def get_dnde(spectrum,energies):
        if isinstance(spectrum,pyLikelihood.Function):
            return BaseGtlikeSED.get_dnde(spectrum,energies)
        elif isinstance(spectrum,Model):
            return spectrum(energies)
        elif isinstance(spectrum,dict):
            return SpectrumPlotter.get_dnde(dict_to_spectrum(spectrum),energies)
        else:
            raise SEDException("Unrecognized type %s for spectrum." % type(spectrum))

    @staticmethod
    def get_dnde_error_gtlike(spectrum,covariance_matrix,energies):
        spectrum = gtlike_unscale_all_parameters(spectrum)

        dnde_err = np.empty_like(energies)
        for i,energy in enumerate(energies):

            # method taken from pyLikelihood.FluxDensity
            srcpars = pyLikelihood.StringVector()
            spectrum.getParamNames(srcpars)
            arg = pyLikelihood.dArg(energy)
            partials = np.array([spectrum.derivByParam(arg, x) for x in srcpars])
            dnde_err[i] = np.sqrt(np.dot(partials, np.dot(covariance_matrix, partials)))
        return dnde_err

    @staticmethod
    def get_dnde_error_pointlike(model,covariance_matrix,energies):
        dnde_err = np.empty_like(energies)
        for i,energy in enumerate(energies):
            partials = model.external_gradient(energy)
            dnde_err[i] = np.sqrt(np.dot(partials, np.dot(covariance_matrix, partials)))
        return dnde_err

    @staticmethod
    def get_dnde_error(spectrum,*args, **kwargs):
        if isinstance(spectrum,pyLikelihood.Function):
            return SpectrumPlotter.get_dnde_error_gtlike(spectrum, *args, **kwargs)
        elif isinstance(spectrum,Model):
            return SpectrumPlotter.get_dnde_error_pointlike(spectrum, *args, **kwargs)
        elif isinstance(spectrum,dict):
            return SpectrumPlotter.get_dnde_error(dict_to_spectrum(spectrum),*args, **kwargs)
        else:
            raise SEDException("Unrecognized type %s for spectrum." % type(spectrum))


    def convert_spectrum(self, function, emin, emax, npts):
        """ Plot a pyLikelihood spectrum onto a matplotlib axes
            whose x axis has units energy_units_obj and y_axes has units
            flux_units_obj/cm^2/s. 
            
            N.B. This function corrects for the fact that gtlike always
            internally represnts the spectrum in units of ph/cm^2/s/MeV. 
            
            
            # kind of ugly, but spectrum is ph/cm^2/s/MeV
            # and it gets mutlplied by energy_units_obj**2,
            # so we need to multiple overall spectrum by
            # flux_units_obj*MeV/energy_units_obj**2
        """

        if emin is None and emax is None:
            emin,emax= self.axes.get_xlim()
        else:
            emin = float(emin/self.axes.energy_units_obj)
            emax = float(emax/self.axes.energy_units_obj)

        energies = np.logspace(np.log10(emin), np.log10(emax), npts)

        e_units = units.tosympy(energies,self.axes.energy_units_obj)

        energies_mev = units.tonumpy(e_units, units.MeV)
        dnde_mev = function(energies_mev)

        # (a) convert to acutal units. gtlike spectra take energy in MeV, return flux in ph/cm^2/s/MeV
        dnde = units.tosympy(dnde_mev,units.ph/units.cm**2/units.s/units.MeV)
        # (b) create E^2 dN/dE in acutal units
        e2_dnde = dnde.multiply_elementwise(e_units).multiply_elementwise(e_units)
        # (c) convert to desired units
        e2_dnde = units.tonumpy(e2_dnde,self.axes.flux_units_obj/units.cm**2/units.s)

        return energies, e2_dnde


    def plot(self, spectrum, emin=None, emax=None, npts=100, autoscale=None, **kwargs):
        energies, e2_dnde = self.convert_spectrum(lambda e: self.get_dnde(spectrum, e), emin, emax, npts)

        if autoscale is not None:
            old_autoscale=self.axes.get_autoscale_on()
            self.axes.autoscale(autoscale)

        self.axes.plot(energies, e2_dnde, **kwargs)

        if autoscale is not None:
            self.axes.autoscale(old_autoscale)

    def plot_error(self, spectrum, covariance_matrix, emin=None, emax=None, npts=100, autoscale=None, **kwargs):
        energies, e2_dnde = self.convert_spectrum(lambda e: self.get_dnde(spectrum, e), emin, emax, npts)
        energies, e2_dnde_error = self.convert_spectrum(lambda e: self.get_dnde_error(spectrum, covariance_matrix, e), emin, emax, npts)
        
        # clip very small values, problems with log scale otherwise
        low=e2_dnde-e2_dnde_error
        low[low<1e-100]=1e-100

        high=e2_dnde+e2_dnde_error

        if autoscale is not None:
            old_autoscale=self.axes.get_autoscale_on()
            self.axes.autoscale(autoscale)

        self.axes.fill_between(energies, low, high, **kwargs)

        if autoscale is not None:
            self.axes.autoscale(old_autoscale)
