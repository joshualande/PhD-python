from os.path import expandvars
from collections import Iterable

import numpy as np

from matplotlib.figure import Figure
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from matplotlib.patheffects import withStroke
from mpl_toolkits.axes_grid.axes_grid import Grid, AxesGrid, ImageGrid

from . arrays import nzip

def histpoints(data, bins):
    """ Create a series of x,y, points which will draw a histogram.
        Useful because the points can be directly plotted with
        the plot function. This is similar to the matplotlib hist
        function, but allows more direct control. 
        
            >>> data = [0.25, 0.5, 1.25]
            >>> bins = [0,1,2]
            >>> x,y=histpoints(data, bins)
            >>> np.all(x == [0, 0, 1, 1, 2, 2])
            True
            >>> np.all(y == [0, 2, 2, 1, 1, 0])
            True
    """
    binned_data, bins = np.histogram(data, bins=bins)

    x = nzip(bins, bins)
    y = np.concatenate(([0], nzip(binned_data, binned_data), [0]))
    return x,y

def plot_ds9_contour(ax,contour,**kwargs):
    """ Parse a ds9 format contour file. Kwargs goes into the plot function. """
    lines=open(expandvars(contour)).readlines()
    ras,decs=[],[]
    for line in lines:
        if line.strip() is '':
            ax['fk5'].plot(ras,decs,'-',**kwargs)
            ras,decs=[],[]
        else:
            ra,dec=line.split()
            ras.append(float(ra)); decs.append(float(dec))

def fix_axesgrid(grid):
    """ Remove the ticks which overlap with nearby axes. """
    if grid._direction != 'row': 
        raise Exception("Not implemented")

    if grid._refax is not None:
        raise Exception("This function does not work for AxesGrids with share_all=True")

    nrows,ncols=grid._nrows,grid._ncols

    for row in range(nrows):
        for col in range(ncols):
            ax = grid[row*ncols + col]
            if row != 0 and col==0:
                ax.set_yticks(ax.get_yticks()[0:-1])
            if col != ncols-1 and row==nrows-1:
                ax.set_xticks(ax.get_xticks()[0:-1])


def label_axes(plots, stroke=True, **kwargs):
    """ Add "(a)" to first plot, "(b)" to second, ... """

    text_kwargs=dict(frameon=False, loc=2, prop=dict(size=14))
    text_kwargs.update(kwargs)

    if isinstance(plots, Iterable) or isinstance(plots, Grid) or \
       isinstance(plots, AxesGrid) or isinstance(plots, ImageGrid):
        plot_list=plots
    elif isinstance(plots,Figure):
        plot_list=plots.axes
    else:
        raise Exception("Unrecognized plot list.")

    for i,g in enumerate(plot_list):
        _at = AnchoredText('(%s)' % chr(i+97), **text_kwargs)

        if stroke:
            _at.txt._text.set_path_effects([withStroke(foreground="w", linewidth=3)])

        g.add_artist(_at)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
