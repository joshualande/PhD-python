import collections

import numpy as np
# not sure why, but works better than matplotlib.axes.Axes
# when using mpl_toolkits.axes_grid.axes_grid.,Grid
from mpl_toolkits.axisartist import HostAxes


import pyLikelihood

# overload SED namespace for our base object
from uw.utilities import keyword_options
from uw.like.Models import Model

from lande.pysed import units
from lande.pysed.helper import logspace_units

from . load import dict_to_spectrum

class SEDException(Exception): 
    pass

class SpectralAxes(HostAxes):
    def __init__(self, *args, **kwargs):
        self.energy_units = kwargs.pop('energy_units', 'MeV')
        self.flux_units = kwargs.pop('flux_units', 'erg')

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

    def plot_points_mev(self, energies_mev, dnde_mev, **kwargs):
        """ Take in energies in MeV and spectal points (in ph/cm^2/s/MeV)
            and plot them on the axes. """

        energies, e2_dnde = self.convert_points_mev(energies_mev, dnde_mev)

        self.plot(energies, e2_dnde, **kwargs)

    def convert_points_mev(self, energies_mev, dnde_mev):

        # (a) convert to acutal units. gtlike spectra take energy in MeV, return flux in ph/cm^2/s/MeV
        return self.convert_points(units.tosympy(energies_mev,units.MeV),
                                   units.tosympy(dnde_mev,units.ph/units.cm**2/units.s/units.MeV))

    def convert_energies(self, energies):
        return units.tonumpy(energies,self.energy_units_obj)

    def convert_points(self, energies, dnde):
        # (b) create E^2 dN/dE in acutal units
        e2_dnde = dnde.multiply_elementwise(energies).multiply_elementwise(energies)
        # (c) convert to desired units
        e2_dnde = units.tonumpy(e2_dnde,self.flux_units_obj/units.cm**2/units.s)

        energies = self.convert_energies(energies)

        return energies, e2_dnde



class SpectrumPlotter(object):
    """ Plot spectra. """
    def __init__(self, axes):
        self.axes = axes

    @staticmethod
    def get_dnde(spectrum,energies):
        """ assume energies has energy units and return flux in flux units. """
        energies=units.tonumpy(energies,units.MeV)
        dnde=SpectrumPlotter.get_dnde_mev(spectrum,energies)
        return units.tosympy(dnde,units.ph/units.cm**2/units.s/units.MeV)

    @staticmethod
    def get_dnde_error(spectrum,covariance_matrix,energies):
        """ assume energies has energy units and return flux in flux units. """
        energies=units.tonumpy(energies,units.MeV)
        dnde_error=SpectrumPlotter.get_dnde_error_mev(spectrum,covariance_matrix,energies)
        return units.tosympy(dnde_error,units.ph/units.cm**2/units.s/units.MeV)

    @staticmethod
    def get_dnde_mev(spectrum,energies):
        """ asume energy in mev and return flux in units of ph/cm**2/s/MeV. """
        if isinstance(spectrum,dict):
            return SpectrumPlotter.get_dnde_mev(dict_to_spectrum(spectrum),energies)
        elif isinstance(spectrum,pyLikelihood.Function):
            return SpectrumPlotter.get_dnde_mev_gtlike(spectrum,energies)
        elif isinstance(spectrum,Model):
            return SpectrumPlotter.get_dnde_mev_pointlike(spectrum,energies)
        else:
            raise SEDException("Unrecognized type %s for spectrum." % type(spectrum))

    @staticmethod
    def get_dnde_mev_gtlike(spectrum,energies):
        """ Returns the spectrum in units of ph/cm^2/s/MeV. """
        if isinstance(energies, collections.Iterable):
            return np.asarray([SpectrumPlotter.get_dnde_mev_gtlike(spectrum,i) for i in energies])
        return spectrum(pyLikelihood.dArg(energies))

    @staticmethod
    def get_dnde_mev_pointlike(spectrum,energies):
        return spectrum(energies)

    @staticmethod
    def get_dnde_error_mev_gtlike(spectrum,covariance_matrix,energies):
        """ asume energy in mev and return flux in units of ph/cm**2/s/MeV. """
        from . models import gtlike_unscale_all_parameters
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
    def get_dnde_error_mev_pointlike(model,covariance_matrix,energies):
        dnde_err = np.empty_like(energies)
        for i,energy in enumerate(energies):
            partials = model.external_gradient(energy)
            dnde_err[i] = np.sqrt(np.dot(partials, np.dot(covariance_matrix, partials)))
        return dnde_err

    @staticmethod
    def get_dnde_error_mev(spectrum, covariance_matrix, energies):
        if isinstance(spectrum,pyLikelihood.Function):
            return SpectrumPlotter.get_dnde_error_mev_gtlike(spectrum, covariance_matrix, energies)
        elif isinstance(spectrum,Model):
            return SpectrumPlotter.get_dnde_error_mev_pointlike(spectrum, covariance_matrix, energies)
        elif isinstance(spectrum,dict):
            return SpectrumPlotter.get_dnde_error_mev(dict_to_spectrum(spectrum), covariance_matrix, energies)
        else:
            raise SEDException("Unrecognized type %s for spectrum." % type(spectrum))

    def plot(self, spectrum, emin=None, emax=None, npts=100, autoscale=None, **kwargs):

        if emin is None and emax is None:
            emin, emax=self.axes.get_xlim()
            emin = emin*self.axes.energy_units_obj
            emax = emax*self.axes.energy_units_obj

        energies = logspace_units(emin, emax, npts)
        dnde = self.get_dnde(spectrum, energies)
        energies, e2_dnde = self.axes.convert_points(energies, dnde)

        if autoscale is not None:
            old_autoscale=self.axes.get_autoscale_on()
            self.axes.autoscale(autoscale)

        self.axes.plot(energies, e2_dnde, **kwargs)

        if autoscale is not None:
            self.axes.autoscale(old_autoscale)

    def plot_error(self, spectrum, covariance_matrix=None, emin=None, emax=None, npts=100, autoscale=None, **kwargs):

        if covariance_matrix is None: covariance_matrix = spectrum['covariance_matrix']

        if emin is None and emax is None:
            emin, emax=self.axes.get_xlim()
            emin = emin*self.axes.energy_units_obj
            emax = emax*self.axes.energy_units_obj

        energies_units = logspace_units(emin, emax, npts)
        dnde = self.get_dnde(spectrum, energies_units)
        dnde_error = self.get_dnde_error(spectrum, covariance_matrix,energies_units)
        energies, e2_dnde = self.axes.convert_points(energies_units, dnde)
        energies, e2_dnde_error = self.axes.convert_points(energies_units, dnde_error)
        
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
