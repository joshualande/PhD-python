import traceback
import sys
import os

from lande.fermi.likelihood.save import get_full_energy_range

from lande.fermi.data.plotting import ROITSMapBandPlotter, ROISourceBandPlotter, ROISourcesBandPlotter

def smooth_plots(roi, name, hypothesis, dirdict, size):
    """ smoothed counts maps """

    extra='%s_%s_%sdeg' % (hypothesis, name, size)

    smooth_kwargs = dict(which=name, 
                         override_center=roi.roi_dir,
                         size=size,
                         colorbar_radius=1)

    title = 'Source Map for %s (%s)' % (name,hypothesis);print title
    roi.plot_source(filename='%s/source_%s.png' % (dirdict['plots'], extra), 
                    title=title,
                    **smooth_kwargs)

    title = 'Sources Map for %s (%s)' % (name,hypothesis);print title
    roi.plot_sources(filename='%s/sources_%s.png' % (dirdict['plots'], extra), 
                     title=title,
                     **smooth_kwargs)

    title = 'Band Source Map for %s (%s)' % (name,hypothesis);print title
    ROISourceBandPlotter(roi, bin_edges=[1e3,1e4,10**5.5], 
                         title=title,
                         **smooth_kwargs).show(filename='%s/band_source_%s.png' % (dirdict['plots'],extra))

    title = 'Band Sources Map for %s (%s)' % (name,hypothesis);print title
    ROISourcesBandPlotter(roi, bin_edges=[1e3,1e4,10**5.5],
                         title=title,
                          **smooth_kwargs).show(filename='%s/band_sources_%s.png' % (dirdict['plots'],extra))
