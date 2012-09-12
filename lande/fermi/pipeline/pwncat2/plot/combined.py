from os.path import expandvars
import pylab as P

from uw.utilities import keyword_options

from lande.pysed import units

from lande.fermi.spectra.sed import SED
from lande.fermi.likelihood.specplot import SpectralAxes, SpectrumPlotter
from lande.fermi.likelihood.limits import UpperLimit
from lande.fermi.likelihood.bandfitter import BandFitter


class CombinedGtlikeSpectralPlotter(object):

    defaults = (
        ('fignum', None, 'matplotlib figure number'),
        ('axes', None, 'matplotlib axes'),
        ('energy_units', 'MeV', 'default units to plot energy flux (y axis) in.'),
        ('flux_units',  'erg', 'default units to plot energy (x axis) in'),
    )


    @keyword_options.decorate(defaults)
    def __init__(self, results, hypothesis, **kwargs):
        keyword_options.process(self, kwargs)

        self.results = results
        self.hypothesis = hypothesis

        self.spectral_model = self.results[self.hypothesis]['gtlike']['spectrum']

        self.cutoff_model = self.results[self.hypothesis]['gtlike']['test_cutoff']['hypothesis_1']['spectrum']

        self.sed_4bpd = self.results[self.hypothesis]['gtlike']['seds']['4bpd']

        self.powerlaw_limit = self.results[self.hypothesis]['gtlike']['powerlaw_upper_limit']
        self.cutoff_limit = self.results[self.hypothesis]['gtlike']['cutoff_upper_limit']

        self.bandfit = self.results[self.hypothesis]['gtlike']['bandfit']

    def plot(self, filename=None, axes=None,
            fignum=None, figsize=(5.5,4.5),
            spectral_kwargs=dict(color='red', zorder=1.9),
            cutoff_kwargs=dict(color='blue', zorder=1.9),
            ):

        if axes is None:
            fig = P.figure(fignum,figsize)
            axes = SpectralAxes(fig=fig,
                                rect=(0.22,0.15,0.75,0.8),
                                flux_units=self.flux_units,
                                energy_units=self.energy_units)
            fig.add_axes(axes)

            axes.set_xlim_units(100*units.MeV, 10**5.5*units.MeV)
    
        # plot sed
        sed = SED(self.sed_4bpd)
        sed.plot(axes=axes,
                 plot_spectral_fit=False, plot_spectral_error=False,
                )

        # plot spectral model
        sp=SpectrumPlotter(axes=axes)
        sp.plot(self.spectral_model, autoscale=False, **spectral_kwargs)
        sp.plot_error(self.spectral_model, autoscale=False, alpha=0.1, **spectral_kwargs)

        sp.plot(self.cutoff_model, autoscale=False, **cutoff_kwargs)
        sp.plot_error(self.cutoff_model, autoscale=False, alpha=0.1, **cutoff_kwargs)

        # plot limits
        ul=UpperLimit(self.powerlaw_limit)
        ul.plot(axes=axes, spectral_kwargs=dict(color='orange', zorder=1.9, autoscale=False))

        ul=UpperLimit(self.cutoff_limit)
        ul.plot(axes=axes, spectral_kwargs=dict(color='purple', zorder=1.9, autoscale=False))

        bf = BandFitter(self.bandfit)
        bf.plot(axes=axes, spectral_kwargs=dict(color='green',zorder=1.9), spectral_error_kwargs=dict(color='green', alpha=0.1))



        if filename is not None:
            P.savefig(expandvars(filename))
        return axes
