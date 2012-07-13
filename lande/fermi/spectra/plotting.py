""" This file contains various function which I have found useful. """
from math import ceil

import pylab as P

from mpl_toolkits.axes_grid1 import AxesGrid

def pointlike_plot_all_seds(roi, filename=None, ncols=4, **kwargs):
    """ Create an SED of all sources in the ROI as a big plot. """
    
    sources=roi.get_sources()
    nrows = int(ceil(float(len(sources))/ncols))
    
    fig = P.figure(figsize=(2.5*ncols,2*nrows),frameon=False)
    tit=P.title("All seds of the sources included in the region.\nRed : fitted sources\nBlue : non fitted sources inside the counts map\nBlack : sources outside of the counts map")

    tit.axes.get_xaxis().set_visible(False)
    tit.axes.get_yaxis().set_visible(False)
    grid = AxesGrid(fig, 111,
                    aspect=False,
                    nrows_ncols = (nrows, ncols),
                    axes_pad = 0.1,
                    add_all=True,
                    label_mode = "L")


    from celgal import dist
    from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
    use_ergs=True
    energy_flux_unit = None
    if energy_flux_unit is None:
        energy_flux_unit = 'erg' if use_ergs else 'MeV'
    assert energy_flux_unit in ('erg', 'MeV', 'GeV', 'eV') , 'unrecognized energy flux unit'
    energy_flux_factor = dict(erg=1.602e-6, MeV=1, eV=1e6, GeV=1e-3)[energy_flux_unit]
    dir=["bottom","top","left","right"]
    for i,which in enumerate(sources):
        source=roi.get_source(which=which)
        distance=dist([source.skydir.ra(),source.skydir.dec()],[roi.sa.roi_dir.ra(),roi.sa.roi_dir.dec()])
        axis = (80, 5e5, 1e-7*energy_flux_factor,3.0e-4*energy_flux_factor)
        axes=grid[i]
        if distance<float(roi.sa.maxROI):
            if len(source.model.get_parameters())!=0:
                for axem in dir:
                    axes.spines[axem].set_color('red')
            else :
                for axem in dir:
                    axes.spines[axem].set_color('blue')
        else :
            for axem in dir:
                axes.spines[axem].set_color('black')
        at = AnchoredText("%s"%which.name,
                          prop=dict(size=8), frameon=True,
                          loc=2,
                          )
        at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
        axes.add_artist(at)
                                
            
        roi.plot_sed(which,axes=grid[i],axis=axis,title=which.name,energy_flux_unit=energy_flux_unit,**kwargs)
        
    if filename is not None: P.savefig(filename)
    

plot_all_seds = pointlike_plot_all_seds # for now
