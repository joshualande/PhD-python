import numpy as np

from SED import SED

from uw.like.Models import PowerLaw

from lande.utilities.tools import tolist

from . models import build_gtlike_spectrum, build_pointlike_model
from . save import fluxdict
from . limits import powerlaw_upper_limit
from . fit import paranoid_gtlike_fit
from . superstate import SuperState


class BandFitter(object):
    """ Performs a gtlike spectral analysis fitting
        the source as a power law in several independent
        energy bins. 

        Note, this code assumes that the intial source (named 'name')
        is a reasonable approximation to the best spectra and
        uses that spectra as a starting value for the fit (but
        allows the fit to vary by a factor of 10^4 in either direciton).
        """
    ul_choices = SED.ul_choices

    def __init__(self, like, name, bin_edges,
                 ul_algorithm='bayesian',
                 ul_confidence=.95,
                 upper_limit_index=2,
                 flux_units='erg'):
        """ Parameters:
            * like - pyLikelihood object
            * name - source to make an SED for
            * bin_edges - if specified, calculate the SED in these bins.
            * ul_algorithm - choices = 'frequentist', 'bayesian' 
            * ul_confidence - confidence level for upper limit. 
            * upper_limit_index - what index to assume when computing upper limits.
            * flux_units - what units to report flux in 
            """
        self.like               = like
        self.name               = name
        self.ul_algorithm       = ul_algorithm
        self.ul_confidence      = ul_confidence
        self.upper_limit_index     = upper_limit_index
        self.flux_units         = flux_units

        source=self.like.logLike.getSource(name) 
        self.init_spectrum=source.spectrum()
        self.init_model=build_pointlike_model(self.init_spectrum)
        self.init_energes = self.like.energies[[0,-1]]

        if not SED.good_binning(self.like, bin_edges):
            raise Exception("bin_edges is not commensurate with the underlying energy binning of pyLikelihood.")
            
        bin_edges = np.asarray(bin_edges)
        self.energy = np.sqrt(bin_edges[1:]*bin_edges[:-1])

        self.lower_energy=bin_edges[:-1]
        self.upper_energy=bin_edges[1:]

        if ul_algorithm not in self.ul_choices:
            raise Exception("Upper Limit Algorithm %s not in %s" % (ul_algorithm,str(self.ul_choices)))

        empty = lambda: np.empty_like(self.energy)

        self.index, self.index_err                = empty(), empty()
        self.flux,   self.flux_err,  self.flux_ul = empty(), empty(), empty()
        self.eflux, self.eflux_err, self.eflux_ul = empty(), empty(), empty()
        self.ts=empty()

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

        for i,(lower,upper) in enumerate(zip(self.lower_energy,self.upper_energy)):
            print 'Calculating spectrum from %.0dMeV to %.0dMeV' % (lower,upper)

            e = np.sqrt(lower*upper)

            like.setEnergyRange(float(lower)+1, float(upper)-1)

            # build new powerlaw


            init_dnde = self.init_model(e)
            model = PowerLaw(norm=init_dnde, index=2, e0=e)
            model.set_limits('norm',1e-4*init_dnde,1e4*init_dnde, scale=init_dnde)
            model.set_limits('index',-5,5)
            spectrum = build_gtlike_spectrum(model)

            like.setSpectrum(name,spectrum)
            like.syncSrcParams(name)

            paranoid_gtlike_fit(like)

            self.ts[i]=like.Ts(name,reoptimize=False, verbosity=4)

            index=spectrum.getParam('Index')
            self.index[i] = index.getTrueValue()
            self.index_err[i] = index.error()*index.getScale()

            fd = fluxdict(like,name,
                          emin=lower,emax=upper,
                          flux_units=self.flux_units)

            self.flux[i] = fd['flux']
            self.flux_err[i] = fd['flux_err']

            self.eflux[i] = fd['eflux']
            self.eflux_err[i] = fd['eflux_err']

            print 'Calculating upper limit from %.0dMeV to %.0dMeV' % (lower,upper)
            ul_dict = powerlaw_upper_limit(like, name, 
                                           powerlaw_index=self.upper_limit_index,
                                           cl=self.ul_confidence,
                                           emin=lower,emax=upper,
                                           flux_units=self.flux_units)
            if ul_dict != None:
                self.flux_ul[i] = ul_dict['flux']
                self.eflux_ul[i] = ul_dict['eflux']
            else:
                self.flux_ul[i] = np.nan
                self.eflux_ul[i] = np.nan

        # revert to old model
        like.setEnergyRange(*self.init_energes)
        saved_state.restore()

    def todict(self):
        """ Pacakge up the results of the SED fit into
            a nice dictionary. """
        return tolist(
            dict(
                name=self.name,
                energy=dict(
                    lower=self.lower_energy,
                    upper=self.upper_energy,
                    value=self.energy,
                    units='MeV'),
                flux=dict(
                    value=self.flux,
                    error=self.flux_err,
                    upper_limit=self.flux_ul,
                    units='ph/cm^2/s'),
                eflux=dict(
                    value=self.eflux,
                    error=self.eflux_err,
                    upper_limit=self.eflux_ul,
                    units='%s/cm^2/s' % self.flux_units),
                index=dict(
                    value=self.index,
                    error=self.index_err),
                TS=self.ts,
                upper_limit_index=self.upper_limit_index,
                ul_confidence=self.ul_confidence,
                )
            )
