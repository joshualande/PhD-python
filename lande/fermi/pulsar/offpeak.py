import numpy as np
from scipy.stats import poisson

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
    def find_phase_range(xx, yy, phases):
        """ xx and yy are the bayesian block decomposition of
            a pulsar light curve. 
            
            The off pulse phase range is defined 
            as the lowest block with 10% removed from either
            side. """
        ranges = [PhaseRange(a,b) for a,b in zip(xx[0::2],xx[1::2])]
        heights = yy[::2]

        if np.allclose(heights[0],heights[-1]):
            ranges[0] += ranges.pop(-1)
            heights = heights[:-1]

        sorted = np.argsort(heights)
        min_phase = ranges[sorted[0]]

        if len(sorted) < 3:
            # if only 2 blocks, no need to merge
            phase = min_phase
        else:

            second_min_phase = ranges[sorted[1]]

            ncounts = len([p for p in phases if p in min_phase])
            second_ncounts = len([p for p in phases if p in second_min_phase])

            predicted_second_counts = ncounts*second_min_phase.phase_fraction/min_phase.phase_fraction

            prob=0.01
            if (poisson.sf(second_ncounts, predicted_second_counts) < prob) or \
               second_min_phase.phase_fraction < 0.5*min_phase.phase_fraction:
                phase = min_phase
            else:
                phase  = min_phase + second_min_phase

        return phase.trim(fraction=0.1)

    @staticmethod
    def get_blocks(phases, ncpPrior, method='binned'):
        """ Apply binned baysian blocks to the pulsar phases
            Return the Bayesian block data. 
        
            According to Jeff Scargal, a easy way to accomplish
            Bayesian blocks on a periodic 
        """

        if method == 'unbinned':
            # use set(...) to remove duplicate phases which break baysian block algorithm.
            # in principle, there shouldn't be any duplicate phases because we
            # are measuring a continuous varaibale. But it is of negligible harm
            # to strip out one or two photons if needed.
            before=len(phases)
            phases=np.sort(list(set(phases)))
            after=len(phases)
            if before != after:
                print 'Warning, stripping %s/%s duplicate photons from the list' % (before-after,before)
        elif method == 'binned':
            if len(phases) < 25:
                nbins=10
            elif len(phases) < 50:
                nbins=25
            elif len(phases) < 100:
                nbins=50
            elif len(phases) < 1000:
                nbins=100
            else:
                nbins=200

        phases = np.concatenate((phases, phases-1, phases+1))

        if method == 'binned':

            tstart=-1
            bins = np.linspace(-1,2,nbins*3+1)
            bin_content = np.histogram(phases, bins=bins)[0]
            bin_sizes = np.ones_like(phases)*(bins[1]-bins[0])
            bb = BayesianBlocks.BayesianBlocks(tstart,bin_content.tolist(),bin_sizes.tolist())
            xx, yy = bb.lightCurve(ncpPrior)

        elif method == 'unbinned':

            bb = BayesianBlocks.BayesianBlocks(phases)
            xx, yy = bb.lightCurve(ncpPrior)

        xx = np.asarray(xx)
        yy = np.asarray(yy)

        assert not np.any(np.isinf(yy))

        bin_size = 1./nbins

        # I ran into floating point problems with blocks that
        # went just a hair above 0, so it would go from like [-0.75, 1e-10].
        # And this would create unphysically small blocks. The
        # easy solution to this is to just add this small tolerance.
        tolerance = 1e-5*bin_size

        first = np.where(xx>0 + tolerance)[0][0]-1
        last = np.where(xx<1 - tolerance)[0][-1]+1

        cut_xx = xx[first:last+1]
        cut_yy = yy[first:last+1]

        # sanity checks
        assert (cut_xx[0] <= tolerance) and (cut_xx[-1] >= 1 - tolerance)

        if not np.allclose(cut_yy[0],cut_yy[-1]):
            print 'Warning, Bayesian blocks do not appear to be periodic!'

        cut_xx[0] = 0
        cut_xx[-1] = 1

        return cut_xx, cut_yy

    def __init__(self,phases,ncpPrior=6):
        """ phases is the numpy array of pulsar phases. """

        self.phases = phases
        self.ncpPrior = ncpPrior

        assert np.all((phases < 1) & (phases >= 0))

        self.xx, self.yy = self.get_blocks(phases, ncpPrior)

        self.off_peak = self.find_phase_range(self.xx, self.yy, self.phases)

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


