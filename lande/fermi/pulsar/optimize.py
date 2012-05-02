import numpy  as np

from uw.pulsar.stats import hm

from lande.fermi.pulsar.data import get_phases

class OptimizePhases(object):
    """ very simple object to load in an ft1 file and
        optimize the radius & energy to find the
        best pulsations. """

    def __init__(self, ft1, skydir, emax,
                 emins=np.linspace(100,1000,21),
                 rads=np.linspace(0.1,2,20),
                 verbose=False,
                ):

        self.ft1 = ft1
        self.skydir = skydir
        self.emax = emax

        stats = np.empty([len(emins),len(rads)])

        for iemin,emin in enumerate(emins):
            for irad,radius in enumerate(rads):
                phases = get_phases(ft1, skydir, emin, emax, radius)
                stat = hm(phases) if len(phases) > 0 else 0
                if verbose: print 'emin=%s, radius=%s, stat=%s, len=%s, n0=%s' % (emin,radius,stat,len(phases),np.sum(phases==0))
                stats[iemin,irad] = stat

        a = np.argmax(stats)
        coord_e, coord_r = np.unravel_index(a, stats.shape)

        self.optimal_emin = emins[coord_e]
        self.optimal_radius = rads[coord_r]
        self.optimal_h = np.max(stats)

