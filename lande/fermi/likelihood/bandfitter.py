import numpy as np

from SED import SED
from lande.utilities.tools import tolist

from . save import fluxdict
from . limits import powerlaw_upper_limit
from . fit import paranoid_gtlike_fit
from . superstate import SuperState


class BandFitter(object):
    """ Performs a gtlike spectral analysis fitting
        the source as a power law in several independent
        energy bins. """
    ul_choices = SED.ul_choices

    def __init__(self, like, name, bin_edges,
                 ul_algorithm='bayesian',
                 ul_confidence=.95,
                 flux_units='erg'):
        """ Parameters:
            * like - pyLikelihood object
            * name - source to make an SED for
            * bin_edges - if specified, calculate the SED in these bins.
            * ul_algorithm - choices = 'frequentist', 'bayesian' 
            * ul_confidence - confidence level for upper limit. """
        self.name               = name
        self.ul_algorithm       = ul_algorithm
        self.ul_confidence      = ul_confidence
        self.flux_units         = flux_units


        if not SED.good_binning(like, bin_edges):
            raise Exception("bin_edges is not commensurate with the underlying energy binning of pyLikelihood.")
            
        bin_edges = np.asarray(bin_edges)
        self.energy = np.sqrt(bin_edges[1:]*bin_edges[:-1])

        self.lower_energy=bin_edges[:-1]
        self.upper_energy=bin_edges[1:]

        if ul_algorithm not in self.ul_choices:
            raise Exception("Upper Limit Algorithm %s not in %s" % (ul_algorithm,str(self.ul_choices)))

        empty = lambda: np.empty_like(self.energy)

        self.index, self.index_err=empty(), empty()
        self.flux, self.flux_err, self.flux_ul=empty(), empty(), empty()
        self.eflux, self.eflux_err, self.eflux_ul=empty(), empty(), empty()
        self.ts=empty()

        self._calculate(like)

    def _calculate(self,like):
        """ Compute the flux data points for each energy. """

        name    = self.name
        init_energes = like.energies[[0,-1]]

        # Freeze all sources except one to make sed of.
        all_sources = like.sourceNames()

        if name not in all_sources:
            raise Exception("Cannot find source %s in list of sources" % name)

        saved_state = SuperState(like)

        like.setSpectrum(name,'PowerLaw')

        # assume a canonical dnde=1e-11 at 1GeV index 2 starting value
        dnde = lambda e: 1e-11*(e/1e3)**-2

        def get(parameter):
            return like[like.par_index(name, parameter)]

        def set(parameter, value, scale, lower, upper):
            p=get(parameter)
            p.setBounds(-1e100, 1e100)
            p.setScale(scale)
            p.setTrueValue(value)
            p.setBounds(lower, upper)
            like.syncSrcParams(name)

        for i,(lower,upper) in enumerate(zip(self.lower_energy,self.upper_energy)):

            e = np.sqrt(lower*upper)

            print 'Calculating spectrum from %.0dMeV to %.0dMeV' % (lower,upper)
            
            set('Prefactor', dnde(e), dnde(e), 1e-10, 1e10)
            set('Index', -2, 1, -5, 5)
            set('Scale', e, 1, -1e-10, 1e10)

            like.setEnergyRange(float(lower)+1, float(upper)-1)

            paranoid_gtlike_fit(like)

            self.ts[i]=like.Ts(name,reoptimize=False, verbosity=4)

            index=get('Index')
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
                                           powerlaw_index=2,
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
        like.setEnergyRange(*init_energes)
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
                )
            )
