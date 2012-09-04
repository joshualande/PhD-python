import numpy as np
from matplotlib.axes import Axes

import pyLikelihood
from SED import SED as BaseGtlikeSED

# overload SED namespace for our base object
from uw.utilities import keyword_options
from uw.like.Models import Model

from lande.pysed import units

from . save import dict_to_spectrum


class SpectralAxes(Axes):
    def __init__(self, energy_units, flux_units, *args, **kwargs):
        self.energy_units = energy_units
        self.flux_units = flux_units

        super(SpectralAxes,self).__init__(*args, **kwargs)
        
        self.set_xscale('log')
        self.set_yscale('log')

        self.set_xlabel('Energy (%s)' % self.energy_units)
        self.set_ylabel('E$^2$ dN/dE (%s cm$^{-2}$ s$^{-1}$)' % self.flux_units)

class SpectrumPlotter(object):
    """ Plot spectra. """
    defaults= (
        ('energy_units', 'MeV', 'default units to plot energy flux (y axis) in.'),
        ('flux_units',  'erg', 'default units to plot energy (x axis) in'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, **kwargs):
        keyword_options.process(self, kwargs)

        self.energy_units_obj     = units.fromstring(self.energy_units)
        self.flux_units_obj       = units.fromstring(self.flux_units)

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
    def _plot_spectrum(spectrum, axes, energy_units_obj, flux_units_obj, emin=None, emax=None, npts=100, **kwargs):
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
            emin,emax= axes.get_xlim()
        energies = np.logspace(np.log10(emin), np.log10(emax), npts)

        e_units = units.tosympy(energies,energy_units_obj)

        energies_mev = units.tonumpy(e_units, units.MeV)

        # (a) convert to acutal units. gtlike spectra take energy in MeV, return flux in ph/cm^2/s/MeV
        dnde = units.tosympy(SpectrumPlotter.get_dnde(spectrum,energies_mev),units.ph/units.cm**2/units.s/units.MeV)
        # (b) create E^2 dN/dE in acutal units
        e2_dnde = dnde.multiply_elementwise(e_units).multiply_elementwise(e_units)
        # (c) convert to desired units
        e2_dnde = units.tonumpy(e2_dnde,flux_units_obj/units.cm**2/units.s)
        axes.plot(energies, e2_dnde, **kwargs)

    def plot(self, spectrum, axes, emin=None, emax=None, **kwargs):
        if axes is None: 
            axes=self.axes
        SpectrumPlotter._plot_spectrum(spectrum, axes, self.energy_units_obj, self.flux_units_obj, emin, emax, **kwargs)

