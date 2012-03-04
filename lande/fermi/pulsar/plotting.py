import pylab as P
import numpy as np

from . data import get_phases, get_phases_and_times
from uw.pulsar.phase_range import PhaseRange

def plot_phaseogram(ft1, nbins=100, filename=None, title=None, off_peak=None, 
                    data_kwargs=dict(),
                    off_peak_kwargs=dict(), repeat_phase=False, axes=None, **kwargs):
    """ Simple code to plot a phaseogram. """
    phases = get_phases(ft1, **kwargs)

    if axes is None:

        fig = P.figure(None, figsize=(5,5))
        axes = fig.add_subplot(111)

    bins = np.linspace(0,1,nbins+1)

    counts,bins=np.histogram(phases,bins=bins)
    x = np.asarray(zip(bins[:-1],bins[1:])).flatten()
    y = np.asarray(zip(counts,counts)).flatten()

    if repeat_phase:
        x,y = np.append(x, x+1), np.append(y, y)

    
    k = dict(color='black')
    k.update(data_kwargs)
    P.plot(x,y, **k)

    axes.set_ylim(0)
    axes.set_xlim(0,2 if repeat_phase else 1)

    if title is not None: 
        axes.set_title(title)

    axes.set_xlabel('Phase')
    axes.set_ylabel('Counts')

    axes.grid=True

    if off_peak is not None:
        kwargs=dict(alpha=0.25, color='blue')
        kwargs.update(off_peak_kwargs)
        PhaseRange(off_peak).axvspan(axes=axes, 
                                     phase_offsets=[0,1] if repeat_phase else 1,
                                     **kwargs)

    if filename is not None:
        P.savefig(filename)

    return axes, bins

def plot_phase_vs_time(ft1, filename, title=None, off_peak=None, 
                       off_peak_kwargs=dict(), **kwargs):
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

    if off_peak is not None:
        kwargs=dict(alpha=0.5, color='white')
        kwargs.update(off_peak_kwargs)
        PhaseRange(off_peak).axvspan(axes=axes, **kwargs)

    axes.set_xlabel('phase')
    axes.set_ylabel('MJD')

    if title is not None: 
        axes.set_title(title)

    P.savefig(filename)
    return axes
