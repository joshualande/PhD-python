from os.path import expandvars
import pylab as P

from uw.utilities import keyword_options

from lande.pysed import units

from lande.fermi.spectra.sed import SED
from lande.fermi.likelihood.specplot import SpectralAxes, SpectrumPlotter
from lande.fermi.likelihood.limits import UpperLimit
from lande.fermi.likelihood.bandfitter import BandFitter


class CombinedSpectralPlotter(object):

    defaults = (
        ('fignum', None, 'matplotlib figure number'),
        ('axes', None, 'matplotlib axes'),
        ('energy_units', 'MeV', 'default units to plot energy flux (y axis) in.'),
        ('flux_units',  'erg', 'default units to plot energy (x axis) in'),
    )


    @keyword_options.decorate(defaults)
    def __init__(self, code, results, hypothesis, **kwargs):
        keyword_options.process(self, kwargs)

        self.code = code
        self.results = results
        self.hypothesis = hypothesis

        assert self.code in ['gtlike','pointlike']

        self._load_results()

    def _load_results(self):
        self.spectral_model = self.results[self.hypothesis][self.code]['spectrum']

        if self.hypothesis != 'extended':
            self.cutoff_model = self.results[self.hypothesis][self.code]['test_cutoff']['hypothesis_1']['spectrum']

        self.sed_4bpd = self.results[self.hypothesis][self.code]['seds']['4bpd']

        if self.hypothesis == 'at_pulsar':
            self.powerlaw_limit = self.results[self.hypothesis][self.code]['powerlaw_upper_limit']
            self.cutoff_limit = self.results[self.hypothesis][self.code]['cutoff_upper_limit']

        #if self.code == 'gtlike':
        #    self.bandfits = self.results[self.hypothesis][self.code]['bandfits']

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
        sed.plot_points(axes=axes)

        if self.hypothesis == 'at_pulsar':
            # plot power-law limits
            ul=UpperLimit(self.powerlaw_limit)
            ul.plot(axes=axes, color='orange', zorder=1.9)

        axes.autoscale(False)

        # plot spectral model
        sp=SpectrumPlotter(axes=axes)
        sp.plot(self.spectral_model, **spectral_kwargs)
        sp.plot_error(self.spectral_model, alpha=0.25, **spectral_kwargs)

        if self.hypothesis != 'extended':
            # plot cutoff model
            sp.plot(self.cutoff_model, **cutoff_kwargs)
            sp.plot_error(self.cutoff_model, alpha=0.25, **cutoff_kwargs)

        if self.hypothesis == 'at_pulsar':
            ul=UpperLimit(self.cutoff_limit)
            ul.plot(axes=axes, color='purple', zorder=1.9)

        #if self.code == 'gtlike':
        #    bf = BandFitter(self.bandfits)
        #    bf.plot(axes=axes, 
        #            spectral_kwargs=dict(color='green',zorder=1.9), 
        #            spectral_error_kwargs=dict(color='green', alpha=0.25))

        if filename is not None:
            P.savefig(expandvars(filename))
        return axes
