import numpy as np
import pylab as P
import yaml

from SED import SED as BaseGtlikeSED

from uw.utilities import keyword_options
from uw.like.Models import PowerLaw

from lande.pysed import units
from lande.utilities.tools import tolist

from . models import build_gtlike_spectrum, build_pointlike_model
from . save import dict_to_spectrum, source_dict
from . limits import GtlikePowerLawUpperLimit
from . fit import paranoid_gtlike_fit
from . superstate import SuperState
from . specplot import SpectrumPlotter,SpectralAxes,set_xlim_mev
from . base import BaseFitter

class BandFitter(BaseFitter):

    defaults= BaseFitter.defaults + (
        ('energy_units', 'MeV', 'default units to plot energy flux (y axis) in.'),
        ('flux_units',  'erg', 'default units to plot energy (x axis) in'),
    )

    def plot(self, filename=None, axes=None, title=None, 
             fignum=None, figsize=(4,4), 
             spectral_kwargs=dict(color='red',zorder=1.9),
             ):

        if axes is None:
            fig = P.figure(fignum,figsize)

            axes = SpectralAxes(fig=fig, rect=(0.22,0.15,0.75,0.8),
                         flux_units=self.flux_units,
                         energy_units=self.energy_units)
            fig.add_axes(axes)
            axes.set_xlim_units(
                self.results['bands'][0]['energy']['emin']*units.fromstring(self.results['bands'][0]['energy']['energy_units']),
                self.results['bands'][-1]['energy']['emax']*units.fromstring(self.results['bands'][-1]['energy']['energy_units'])
            )
            
        sp = SpectrumPlotter(axes=axes)

        for i,r in enumerate(self.results['bands']):
            emin=r['energy']['emin']*units.fromstring(r['energy']['energy_units'])
            emax=r['energy']['emax']*units.fromstring(r['energy']['energy_units'])
            spectrum = dict_to_spectrum(r['spectrum'])

            if i > 0 and 'label' in spectral_kwargs: 
                # only one label
                spectral_kwargs.pop('label') 

            sp.plot(spectrum, emin=emin, emax=emax, **spectral_kwargs)

        if title is not None: axes.set_title(title)
        if filename is not None: P.savefig(filename)
        return axes



class GtlikeBandFitter(BandFitter):
    """ Performs a gtlike spectral analysis fitting
        the source as a power law in several independent
        energy bins. 

        Note, this code assumes that the intial source (named 'name')
        is a reasonable approximation to the best spectra and
        uses that spectra as a starting value for the fit (but
        allows the fit to vary by a factor of 10^4 in either direciton).
        """
    ul_choices = BaseGtlikeSED.ul_choices

    defaults = BandFitter.defaults + (
        ('ul_algorithm','bayesian',"choices = 'frequentist', 'bayesian' "),
        ('ul_confidence',.95,'confidence level for upper limit.'),
        ('upper_limit_index',2,'what index to assume when computing upper limits.'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, like, name, bin_edges, **kwargs):
        """ Parameters:
            * like - pyLikelihood object
            * name - source to make an SED for
            * bin_edges - if specified, calculate the SED in these bins.
        """
        keyword_options.process(self, kwargs)
        self.like               = like
        self.name               = name

        if not BaseGtlikeSED.good_binning(self.like, bin_edges):
            raise Exception("bin_edges is not commensurate with the underlying energy binning of pyLikelihood.")

        source=self.like.logLike.getSource(name) 
        self.init_spectrum=source.spectrum()
        self.init_model=build_pointlike_model(self.init_spectrum)
        self.init_energes = self.like.energies[[0,-1]]
            
        bin_edges = np.asarray(bin_edges)

        self.lower_energy=bin_edges[:-1]
        self.upper_energy=bin_edges[1:]
        self.middle_energy=np.sqrt(self.lower_energy*self.upper_energy)

        if self.ul_algorithm not in self.ul_choices:
            raise Exception("Upper Limit Algorithm %s not in %s" % (self.ul_algorithm,str(self.ul_choices)))

        empty = lambda: np.empty_like(self.middle_energy)

        self._calculate()

    def _calculate(self):
        """ Compute the flux data points for each energy. """

        like         = self.like
        name         = self.name

        # Freeze all sources except one to make sed of.
        all_sources = like.sourceNames()

        if name not in all_sources:
            raise Exception("Cannot find source %s in list of sources" % name)

        saved_state = SuperState(like)

        self.results = dict(
            name=name,
            bands=[],
        )

        for i,(emin,emax,e_middle) in enumerate(zip(self.lower_energy,self.upper_energy,self.middle_energy)):
            if self.verbosity: print 'Calculating spectrum from %.0dMeV to %.0dMeV' % (emin,emax)


            like.setEnergyRange(float(emin)+1, float(emax)-1)

            # build new powerlaw, assuming initial model has a reasonable spectrum
            init_dnde = self.init_model(e_middle)
            model = PowerLaw(norm=init_dnde, index=2, e0=e_middle)
            model.set_limits('norm',1e-4*init_dnde,1e4*init_dnde, scale=init_dnde)
            model.set_limits('index',-5,5)
            spectrum = build_gtlike_spectrum(model)

            like.setSpectrum(name,spectrum)
            like.syncSrcParams(name)

            paranoid_gtlike_fit(like, verbosity=self.verbosity)

            r = source_dict(like, name, emin=emin, emax=emax,
                            flux_units=self.flux_units,
                            energy_units=self.energy_units,
                            verbosity=self.verbosity)

            if self.verbosity: print 'Calculating upper limit from %.0dMeV to %.0dMeV' % (emin,emax)
            g = GtlikePowerLawUpperLimit(like, name,
                                         powerlaw_index=self.upper_limit_index,
                                         cl=self.ul_confidence,
                                         emin=emin,emax=emax,
                                         flux_units=self.flux_units,
                                         energy_units=self.energy_units)
            r['upper_limit'] = g.todict()
            self.results['bands'].append(r)

        # revert to old model
        like.setEnergyRange(*self.init_energes)
        saved_state.restore()
