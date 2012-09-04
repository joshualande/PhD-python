""" Implements a subclass of SED which has nicer features. 

    Author: J. Lande
"""
import pylab as P
import numpy as np

# overload SED namespace for our base object
from SED import SED as BaseGtlikeSED

from lande.pysed import units

from uw.utilities import keyword_options

from lande.fermi.likelihood.base import BaseFitter
from lande.fermi.likelihood.specplot import SpectralAxes, SpectrumPlotter, set_xlim_mev

class SEDException(Exception): 
    pass


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

    @keyword_options.decorate(defaults)
    def __init__(self, results, *args, **kwargs):
        keyword_options.process(self, kwargs)

        self.energy_units_obj = units.fromstring(self.energy_units)
        self.flux_units_obj = units.fromstring(self.flux_units)

        super(SED,self).__init__(results)


    def plot(self, filename=None, axes=None, title=None,
             fignum=None, figsize=(4,4),
             plot_spectral_fit=True,
             data_kwargs=dict(),
             spectral_kwargs=dict(color='red',zorder=1.9)):
        """ Plot the SED using matpotlib. """
        edict = self.results['Energy']
        ce = lambda x: units.convert(np.asarray(x),units.fromstring(edict['Units']), self.energy_units_obj)

        # get energy part
        energy = ce(edict['Value'])
        if 'Lower' in edict and 'Upper' in edict:
            lower_energy = ce(edict['Lower'])
            upper_energy = ce(edict['Upper'])
            has_energy_errors = True
        else:
            has_energy_errors = False

        # get spectral part

        fdict = self.results['dNdE']

        cf = lambda y: units.convert(energy**2*np.asarray(y),
                                     self.energy_units_obj**2*units.fromstring(fdict['Units']),
                                     self.flux_units_obj/units.cm**2/units.s)

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


        if axes is None:
            fig = P.figure(fignum,figsize)
            axes = SpectralAxes(fig=fig, 
                                rect=(0.22,0.15,0.75,0.8),
                                flux_units=self.flux_units,
                                energy_units=self.energy_units)
            fig.add_axes(axes)

            if has_energy_errors:
                # use BaseGtlikeSED.set_xlim to add 10% on either side.
                set_xlim_mev(axes, lower_energy[0], upper_energy[-1], self.energy_units)
            else:
                set_xlim_mev(axes, energy[0], energy[-1], self.energy_units)

        self.axes = axes

        if plot_spectral_fit and 'spectrum' in self.results:
            sp=SpectrumPlotter(energy_units=self.energy_units, flux_units=self.flux_units)
            sp.plot(self.results['spectrum'], axes=axes, **spectral_kwargs)


        BaseGtlikeSED._plot_points(
            x=energy,
            xlo=lower_energy if has_energy_errors else None,
            xhi=upper_energy if has_energy_errors else None,
            y=dnde,
            y_lower_err=dnde_lower_err if has_assymetric_errors else dnde_err,
            y_upper_err=dnde_upper_err if has_assymetric_errors else dnde_err,
            y_ul=dnde_ul if has_upper_limits else None,
            significant=significant if has_upper_limits else np.ones(len(energy),dtype=bool),
            axes=axes, **data_kwargs)

        if title is not None: axes.set_title(title)
        if filename is not None: P.savefig(filename)
        return axes



