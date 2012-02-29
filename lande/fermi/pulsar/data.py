from os.path import expandvars

import pylab as P
import numpy as np

from uw.utilities.fitstools import rad_extract
from uw.pulsar.phase_range import PhaseRange

def plot_phaseogram(ft1, nbins=100, filename=None, title=None, off_pulse=None, axes=None, **kwargs):
    """ Simple code to plot a phaseogram. """
    phases = get_phases(ft1, **kwargs)

    if axes is None:

        fig = P.figure(None, figsize=(5,5))
        axes = fig.add_subplot(111)

    bins = np.linspace(0,1,nbins+1)
    axes.hist(phases,bins=bins,histtype='step',ec='red',normed=True,lw=1)
    axes.set_xlim(0,1)

    if title is not None: 
        axes.set_title(title)

    axes.set_xlabel('Phase')
    axes.set_ylabel('Counts')

    axes.grid=True

    if off_pulse is not None:
        PhaseRange(off_pulse).axvspan(axes=axes, alpha=0.5, color='blue')

    if filename is not None:
        P.savefig(filename)

    return axes, bins

def plot_phase_vs_time(ft1, filename, title=None, off_pulse=None, **kwargs):
    """ Simple code to plot phase vs time. """
    phases, times = get_phases_and_times(ft1, **kwargs)

    # here, put a 2d histogram
    fig = P.figure(None, figsize=(5,5))
    fig.subplots_adjust(left=0.2)
    axes = fig.add_subplot(111)

    # Note about 2D histograms: 
    #  http://www.physics.ucdavis.edu/~dwittman/Matplotlib-examples/
    hist, xedges, yedges = np.histogram2d(phases, times, bins=(50,50), range=[[0,1], [min(times), max(times)]])
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1] ]
    axes.imshow(hist.T,extent=extent,interpolation='nearest',origin='lower', aspect='auto')

    if off_pulse is not None:
        PhaseRange(off_pulse).axvspan(axes=axes, alpha=0.5, color='white')

    axes.set_xlabel('phase')
    axes.set_ylabel('MJD')

    if title is not None: 
        axes.set_title(title)

    P.savefig(filename)
    return axes

from uw.like.roi_image import memoize
@memoize
def get_all_phases(ft1, skydir):
    """ Cache photons = faster """
    ed = rad_extract(expandvars(ft1),skydir,radius_function=180,return_cols=['PULSE_PHASE', 'TIME'])
    return ed


def get_phases_and_times(ft1, skydir, emin, emax, radius):
    ed = get_all_phases(ft1, skydir)
    all_phases = ed['PULSE_PHASE']
    all_times = ed['TIME']

    cut = (ed['ENERGY'] >= emin) & (ed['ENERGY'] <= emax) & (ed['DIFFERENCES'] <= np.radians(radius))
    return all_phases[cut], all_times[cut]

def get_phases(*args, **kwargs):
    return get_phases_and_times(*args, **kwargs)[0]



