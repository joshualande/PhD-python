""" Implements a subclass of SED which has nicer features. 

    Author: J. Lande
"""
from os.path import expandvars

import pylab as P
import numpy as np

# overload SED namespace for our base object
from lande.pysed import units

from uw.utilities import keyword_options

from lande.utilities.plotting import plot_points

from lande.fermi.likelihood.basefit import BaseFitter
from lande.fermi.likelihood.specplot import SpectralAxes, SpectrumPlotter

class SED(BaseFitter):
    """ Base object for plotting SEDs.

        The input must be a YAML-formatted text file (or python dictionary).

        TODO: DOCUMENT FORMAT

        Something about optional significant flag, what to do with asymetrical errors.

    """
    defaults = BaseFitter.defaults + (
        ('energy_units', 'MeV', 'default units to plot energy flux (y axis) in.'),
        ('flux_units',  'erg', 'default units to plot energy (x axis) in'),
    )

    def plot_spectral_error(self, axes, **kwargs):
        spectral_error_kwargs=dict(color='red',alpha=0.5,zorder=1.8)
        spectral_error_kwargs.update(kwargs)

        sp=SpectrumPlotter(axes=axes)
        sp.plot_error(self.results['spectrum'], self.results['spectrum']['covariance_matrix'],
                      autoscale=False, **spectral_error_kwargs)

    def plot_spectral_fit(self,axes, **kwargs):

        spectral_kwargs=dict(color='red',zorder=1.9)
        spectral_kwargs.update(kwargs)

        sp=SpectrumPlotter(axes=axes)
        sp.plot(self.results['spectrum'], 
                autoscale=False, **spectral_kwargs)

    def plot_points(self, axes, **kwargs):
        """ Plot the SED using matpotlib. """
        data_kwargs=dict(color='black')
        data_kwargs.update(kwargs)


        edict = self.results['Energy']
        fdict = self.results['dNdE']

        file_energy_units = units.fromstring(edict['Units'])
        file_flux_units = units.fromstring(fdict['Units'])

        ce = lambda x: units.convert(np.asarray(x),file_energy_units, axes.energy_units_obj)

        # get energy part
        energy = ce(edict['Value'])
        if 'Lower' in edict and 'Upper' in edict:
            lower_energy = ce(edict['Lower'])
            upper_energy = ce(edict['Upper'])
            has_energy_errors = True
        else:
            has_energy_errors = False

        # get spectral part

        cf = lambda y: units.convert(energy**2*np.asarray(y),
                                     axes.energy_units_obj**2*file_flux_units,
                                     axes.flux_units_obj/units.cm**2/units.s)

        dnde = cf(fdict['Value'])

        if 'Lower_Error' in fdict and 'Upper_Error' in fdict:
            # assymetric errors
            dnde_lower_err = cf(fdict['Lower_Error'])
            dnde_upper_err = cf(fdict['Upper_Error'])
            has_assymetric_errors = True
        else:
            has_assymetric_errors = False
            dnde_err = cf(fdict['Average_Error'])

        # get limits, otherwise assume all significant
        if 'Upper_Limit' in fdict and 'Significant' in self.results:
            dnde_ul = cf(fdict['Upper_Limit'])
            significant = np.asarray(self.results['Significant'])
            has_upper_limits=True
        else:
            has_upper_limits=False

        plot_points(
            x=energy,
            xlo=lower_energy if has_energy_errors else None,
            xhi=upper_energy if has_energy_errors else None,
            y=dnde,
            y_lower_err=dnde_lower_err if has_assymetric_errors else dnde_err,
            y_upper_err=dnde_upper_err if has_assymetric_errors else dnde_err,
            y_ul=dnde_ul if has_upper_limits else None,
            significant=significant if has_upper_limits else np.ones(len(energy),dtype=bool),
            axes=axes, **data_kwargs)

    def plot(self, filename=None, axes=None, title=None,
             fignum=None, figsize=(4,4),
             plot_spectral_fit=True,
             plot_spectral_error=True,
             data_kwargs=dict(),
             spectral_kwargs=dict(),
             spectral_error_kwargs=dict(),
            ):

        if axes is None:
            fig = P.figure(fignum,figsize)
            axes = SpectralAxes(fig=fig, 
                                rect=(0.22,0.15,0.75,0.8),
                                flux_units=self.flux_units,
                                energy_units=self.energy_units)
            fig.add_axes(axes)

            edict = self.results['Energy']
            file_energy_units = units.fromstring(edict['Units'])

            if 'Lower' in edict and 'Upper' in edict:
                axes.set_xlim_units(edict['Lower'][0]*file_energy_units, edict['Upper'][-1]*file_energy_units)
            else:
                axes.set_xlim_units(edict['Energy'][0]*file_energy_units, edict['Energy'][-1]*file_energy_units)

        self.plot_points(axes=axes, **data_kwargs)

        if plot_spectral_fit:
            self.plot_spectral_fit(axes=axes, **spectral_kwargs)
        if plot_spectral_error:
            self.plot_spectral_error(axes=axes, **spectral_error_kwargs)

        if title is not None: axes.set_title(title)
        if filename is not None: 
            P.savefig(expandvars(filename))
        return axes
