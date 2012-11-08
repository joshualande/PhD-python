import traceback
import sys

import numpy as np
from scipy.stats import poisson

import BayesianBlocks

from uw.pulsar.phase_range import PhaseRange

from . plotting import plot_phaseogram


class BlockException(Exception):
    pass

class PeriodicBlocks(object):
    """ Applies Baysian blocks to
        data which is periodic between 0 and 1. """

    def __init__(self, phases, ncpPrior):
        if np.any((phases >= 1)|(phases < 0)):
            raise BlockException("Bad phases %s" % phases[(phases<0)|(phases>=1)])

        self.phases = phases
        self.ncpPrior = ncpPrior

        self.xx, self.yy = self.get_blocks()

        # Array of phase range for each block
        self.ranges = [PhaseRange(a,b) for a,b in zip(self.xx[0::2],self.xx[1::2])]

        # array of heights
        self.heights = self.yy[::2]

    def get_blocks(self):
        """ Apply binned baysian blocks to the pulsar phases
            Return the Bayesian block data. 
        
            According to Jeff Scargal, a easy way to accomplish
            Bayesian blocks on a periodic 
        """
        loop_phases = np.concatenate((self.phases, self.phases-1, self.phases+1))

        xx, yy = self._get_loop_blocks(loop_phases, self.phases)

        xx = np.asarray(xx)
        yy = np.asarray(yy)

        if np.any(np.isinf(yy)):
            raise BlockException('Error, some yy are inf. xx=%s, yy=%s' % (str(xx),str(yy)))


        bin_size = 0.0001

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


class PeriodicUnbinnedBlocks(PeriodicBlocks):

    def _get_loop_blocks(self, loop_phases, phases):
            # use set(...) to remove duplicate phases which break baysian block algorithm.
            # in principle, there shouldn't be any duplicate phases because we
            # are measuring a continuous varaibale. But it is of negligible harm
            # to strip out one or two photons if needed.

        before=len(loop_phases)
        loop_phases=np.sort(list(set(loop_phases)))
        after=len(loop_phases)
        if before != after:
            print 'Warning, stripping %s/%s duplicate photons from the list' % (before-after,before)

        bb = BayesianBlocks.BayesianBlocks(loop_phases)
        xx, yy = bb.lightCurve(self.ncpPrior)

        return xx, yy


class PeriodicBinnedBlocks(PeriodicBlocks):

    def _get_loop_blocks(self, loop_phases, phases):
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

        print 'phases',loop_phases

        tstart=-1
        bins = np.linspace(-1,2,nbins*3+1)
        bin_content = np.histogram(loop_phases, bins=bins)[0]
        bin_sizes = np.ones_like(loop_phases)*(bins[1]-bins[0])
        bb = BayesianBlocks.BayesianBlocks(tstart,bin_content.tolist(),bin_sizes.tolist())
        xx, yy = bb.lightCurve(self.ncpPrior)

        return xx, yy



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
    def find_phase_range(periodic_blocks, phases, prob_2peak=0.01):
        """ xx and yy are the bayesian block decomposition of
            a pulsar light curve. 
            
            The off pulse phase range is defined 
            as the lowest block with 10% removed from either
            side. """
        xx = periodic_blocks.xx
        yy = periodic_blocks.yy

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

            prob =poisson.sf(second_ncounts, predicted_second_counts)
            print 'Probability of there being a second peak from %s is %s' % (str(second_min_phase),prob)
            region_too_small = second_min_phase.phase_fraction < 0.5*min_phase.phase_fraction
            height_too_different = prob < prob_2peak

            if height_too_different or region_too_small:
                if height_too_different:
                    print 'Rejecting second peak - heights are inconsistent'
                if region_too_small:
                    print 'Rejecting second peak - region too small'
                phase = min_phase
            else:
                print 'Adding second peak!'
                phase  = min_phase + second_min_phase

        return phase.trim(fraction=0.1)


    def __init__(self,phases,ncpPrior=8):
        """ phases is the numpy array of pulsar phases. """

        self.phases = phases
        self.ncpPrior = ncpPrior
        self.actual_ncpPrior = ncpPrior
        
        print 'Trying Baysian blocks with ncpPrior=%s' % self.ncpPrior
        while True:
            try:
                self.periodic_blocks = PeriodicBinnedBlocks(phases=phases, ncpPrior=self.actual_ncpPrior)

                print '  xx, yy', self.periodic_blocks.xx, self.periodic_blocks.yy
            except BlockException, ex:
                print 'Baysian blocks failed with ncpPrior=%s because PeriodicBinnedBlocks raised an exception' % self.actual_ncpPrior
                traceback.print_exc(file=sys.stdout)
                good_blocks = False
            else:
                if len(self.periodic_blocks.xx) <= 2:
                    print 'Baysian blocks failed with ncpPrior=%s because the entire interval is constant' % self.actual_ncpPrior
                    good_blocks = False
                elif np.any(np.isinf(self.periodic_blocks.yy)):
                    print 'Baysian blocks failed with ncpPrior=%s because some of the blocks have infinite height' % self.actual_ncpPrior
                    good_blocks = False
                else:
                    good_blocks = True

            if good_blocks:
                break
            else:
                print 'Bayesian blocks failed with ncpPrior=%s, trying again with ncpPrior=%s' % (self.actual_ncpPrior,self.actual_ncpPrior-1)
                self.actual_ncpPrior -= 1 
                if self.actual_ncpPrior < 1:
                    raise Exception("No good ncpPriors!")


        self.off_peak = self.find_phase_range(self.periodic_blocks, phases)

        self.blocks = dict(xx = self.periodic_blocks.xx, yy = self.periodic_blocks.yy)


def plot_phaseogram_blocks(ft1, blocks=None, blocks_kwargs=dict(), repeat_phase=False, **kwargs):

    # plot bins
    axes, bins = plot_phaseogram(ft1, repeat_phase=repeat_phase, **kwargs)
    binsz = bins[1]-bins[0]

    if blocks is not None:
        # plot blocks
        xx=np.asarray(blocks['xx'])
        yy=np.asarray(blocks['yy'])
        if repeat_phase: xx, yy = np.append(xx, xx+1), np.append(yy, yy)

        k=dict(color='blue')
        k.update(blocks_kwargs)
        axes.plot(np.asarray(xx),np.asarray(yy)*binsz, **k)


