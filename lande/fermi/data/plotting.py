import os

import pylab as P
import numpy as np
import pywcsgrid2

from mpl_toolkits.axes_grid1.axes_grid import ImageGrid

from uw.utilities import keyword_options
 
from uw.like.roi_plotting import ROITSMapPlotter, ROISmoothedSources, ROISmoothedSource, ROISignificance
from uw.like.roi_state import PointlikeState

class ROIBandPlotter(object):

    def __init__(self,roi,bin_edges,**kwargs):

        self.roi = roi
        keyword_options.process(self, kwargs)

        self.bin_edges = bin_edges
        self.nplots = len(bin_edges) - 1

        for e in bin_edges:
            if not np.any(np.abs(e-roi.bin_edges) < 0.5):
                raise Exception("Energy %.1f inconsistent with ROI energy binning." % e)

        self.lower_energies = bin_edges[:-1]
        self.upper_energies = bin_edges[1:]

        state = PointlikeState(roi)
 
        # step 1, test consistentcy of each energy with binning in pointlike

        self.maps = []
        for i,(lower,upper) in enumerate(zip(self.lower_energies, self.upper_energies)):
            roi.change_binning(fit_emin=lower,fit_emax=upper)
            self.maps.append(self.object(roi,title='',**kwargs))

        state.restore()

    def show(self,filename=None):
        print 'figsize = ',self.figsize

        self.fig = fig = P.figure(self.fignum,self.figsize)
        P.clf()

        header = self.maps[0].header

        self.grid = grid = ImageGrid(fig, (1, 1, 1), 
                                     nrows_ncols = (1, self.nplots),
                                     axes_pad=0.1, share_all=True,
                                     cbar_location="top",
                                     cbar_mode="each",
                                     cbar_size="7%",
                                     cbar_pad="2%",
                                     axes_class=(pywcsgrid2.Axes,
                                                 dict(header=header)))

        for i,(map,lower,upper) in enumerate(zip(self.maps,self.lower_energies,self.upper_energies)):
            map.show(axes=grid[i], cax=grid[i].cax)
            format_energy=lambda x: '%.1f' % (x/1000.) if x < 1000 else '%.0f' % (x/1000.)
            lower_string=format_energy(lower)
            upper_string=format_energy(upper)
            grid[i].add_inner_title("%s-%s GeV" % (lower_string,upper_string), loc=2)

        if filename is not None: P.savefig(filename)

class ROITSMapBandPlotter(ROIBandPlotter):
    object = ROITSMapPlotter
    defaults = object.defaults 
    defaults=keyword_options.change_defaults(defaults,'figsize',(9,4))

class ROISourcesBandPlotter(ROIBandPlotter):
    object = ROISmoothedSources
    defaults = object.defaults 
    defaults=keyword_options.change_defaults(defaults,'figsize',(9,4))

class ROISourceBandPlotter(ROIBandPlotter):
    object = ROISmoothedSource
    defaults = object.defaults 
    defaults=keyword_options.change_defaults(defaults,'figsize',(9,4))

class ROISignificanceBandPlotter(ROIBandPlotter):
    object = ROISignificance
    defaults = object.defaults 
    defaults=keyword_options.change_defaults(defaults,'figsize',(9,4))





