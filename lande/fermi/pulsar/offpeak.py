import numpy as np

import BayesianBlocks

from uw.pulsar.phase_range import PhaseRange

from . plotting import plot_phaseogram


class OffPeakBB(object):
    """ Algorithm to compute the off peak window
        of a pulsar using a Bayesian block analysis. 
        
        N.B. I ran into some preliminary problems using
        this algorithm running Bayesian blocks in
        unbinned mode becasue it would try to find
        very small structure in the light curves
        which I knew to be unphysical.
        The reasonable solution I decided upon was
        to run Bayesian Blocks in binned mode
        with a binsize of 1/100 in pulsar phase.
        This washes out any structure on smaller
        scales and appaers to help
        the method only find physically large enough
        structure in the pulsar light curve. """

    @staticmethod
    def find_phase_range(xx, yy):
        """ xx and yy are the bayesian block decomposition of
            a pulsar light curve. 
            
            The off pulse phase range is defined 
            as the lowest block with 10% removed from either
            side. """

        # first, get phase range by finding the lowest valey
        min_bin = np.argmin(yy)
        phase_min = xx[min_bin]
        phase_max = xx[min_bin+1]

        # worry about possibility of wraping around
        if (min_bin == 0) and np.allclose(yy[0],yy[-1]):
            phase_min = xx[-2] - 1

        phase_range = phase_max - phase_min
        phase_min, phase_max = phase_min + 0.1*phase_range, phase_max - 0.1*phase_range

        return PhaseRange(phase_min, phase_max)

    @staticmethod
    def get_binned_blocks(phases, ncpPrior):
        """ Apply binned baysian blocks to the pulsar phases
            Return the Bayesian block data. 
        
            According to Jeff Scargal, a easy way to accomplish
            Bayesian blocks on a periodic 
        """

        phases = np.concatenate((phases, phases-1, phases+1))

        nbins=200

        tstart=-1
        bins = np.linspace(-1,2,nbins*3+1)
        bin_content = np.histogram(phases, bins=bins)[0]
        bin_sizes = np.ones_like(phases)*(bins[1]-bins[0])
        bb = BayesianBlocks.BayesianBlocks(tstart,bin_content.tolist(),bin_sizes.tolist())
        xx, yy = bb.lightCurve(ncpPrior)

        xx = np.asarray(xx)
        yy = np.asarray(yy)

        assert not np.any(np.isinf(yy))

        first = np.where(xx>0)[0][0]-1
        last = np.where(xx<1)[0][-1]+1

        cut_xx = xx[first:last+1]
        cut_yy = yy[first:last+1]

        # sanity checks

        assert (cut_xx[0] <= 0) and (cut_xx[-1] >= 1)

        if not np.allclose(cut_yy[0],cut_yy[-1]):
            print 'Warning, Bayesian blocks do not appear to be periodic!'

        cut_xx[0] = 0
        cut_xx[-1] = 1

        return cut_xx, cut_yy

    def __init__(self,phases,ncpPrior=5):
        """ phases is the numpy array of pulsar phases. """

        self.phases = phases
        self.ncpPrior = ncpPrior

        assert np.all((phases < 1) & (phases >= 0))

        self.xx, self.yy = self.get_binned_blocks(phases, ncpPrior)

        self.off_peak = self.find_phase_range(self.xx, self.yy)

        self.blocks = dict(xx = self.xx, yy = self.yy)


def plot_phaseogram_blocks(ft1, blocks, blocks_kwargs=dict(), repeat_phase=False, **kwargs):

    # plot bins
    axes, bins = plot_phaseogram(ft1, repeat_phase=repeat_phase, **kwargs)
    binsz = bins[1]-bins[0]

    # plot blocks
    xx=blocks['xx']
    yy=blocks['yy']
    if repeat_phase: xx, yy = np.append(xx, xx+1), np.append(yy, yy)

    k=dict(color='blue')
    k.update(blocks_kwargs)
    axes.plot(np.asarray(xx),np.asarray(yy)*binsz, **k)


