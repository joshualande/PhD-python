from os.path import expandvars
import pylab as P

from uw.utilities import keyword_options

from lande.fermi.spectra.sed import SED


class CombinedPlotter(object):

    defaults = (
        ('fignum', None, 'matplotlib figure number'),
        ('axes', None, 'matplotlib axes'),
        ('figsize',             (5.5,4.5), 'Size of figure in inches'),
        ('energy_units', 'MeV', 'default units to plot energy flux (y axis) in.'),
        ('flux_units',  'erg', 'default units to plot energy (x axis) in'),
    )


    @keyword_options.decorate(defaults)
    def __init__(self, results, hypothesis, **kwargs):
        keyword_options.process(self, kwargs)

        self.results = results
        self.hypothesis = hypothesis

        if self.axes is None:
            fig = P.figure(self.fignum,self.figsize)
            P.clf()
            self.axes = fig.add_subplot(111)

        self.plot_sed()
        self.plot_bandfit()
        self.plot_spectra()
        self.plot_ul()

    def plot_bandfit(self):
        pass

    def plot_spectra(self):
        pass

    def plot_sed(self):
        sed = SED(
            self.results[self.hypothesis]['gtlike']['4bpd'],
            energy_units=self.energy_units,
            flux_units=self.flux_units,
        )
        sed.plot(axes=self.axes)

        pass

    def plot_ul(self):
        pass

    def save(self, filename):
        P.savefig(expandvars(filename))
