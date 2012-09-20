import traceback
import sys
import os

from lande.fermi.likelihood.save import get_full_energy_range

from lande.fermi.data.plotting import ROITSMapBandPlotter, ROISourceBandPlotter, ROISourcesBandPlotter

from . binning import all_energy, one_bin_per_dec

def overlay_on_plot(axes, pulsar_position):
    """ Function to overlay on all plots
        * The pulsar position
        * New non-2FGL sources addded to the ROI. """
    axes['gal'].plot([pulsar_position.l()],[pulsar_position.b()], 
                     marker='*', color='green',
                     markeredgecolor='white', markersize=12, zorder=1)


def tsmap_plots(roi, name, hypothesis, datadir, plotdir, size, new_sources, tsmap_pixelsize=0.1, **common_kwargs):
    """ TS maps """
    emin, emax = get_full_energy_range(roi)

    extra='%s_%s_%sdeg' % (hypothesis,name,size)

    tsmap_kwargs = dict(size=size, pixelsize=tsmap_pixelsize, **common_kwargs)

    def _plot(filename,title):
        try:
            print 'Making Band %s TS Map (size=%s, pixelsize=%s)' % (title,size,tsmap_pixelsize)
            roi.plot_tsmap(filename='%s/tsmap_%s_%s.png' % (plotdir,filename,extra), 
                           fitsfile='%s/tsmap_%s_%s.fits' % (datadir,filename,extra),
                           title='%s TS Map for %s (%s)' % (title,name,hypothesis),
                           **tsmap_kwargs)
        except Exception, ex:
            print 'ERROR making TS map:', ex
            traceback.print_exc(file=sys.stdout) 

        if all_energy(emin,emax):
            try:
                print 'Making Band %s TS Map (size=%s, pixelsize=%s)' % (title,size,tsmap_pixelsize)
                ROITSMapBandPlotter(roi,  
                                    title='Band %s TS Map for %s (%s)' % (title,name,hypothesis),
                                    bin_edges=one_bin_per_dec(emin,emax),
                                    **tsmap_kwargs).show(filename='%s/band_tsmap_%s_%s.png' % (plotdir,filename,extra))
            except Exception, ex:
                print 'ERROR making band TS map:', ex
                traceback.print_exc(file=sys.stdout) 

    # reisidual ts map
    _plot('residual','Residual')

    # source ts map
    roi.zero_source(which=name)
    _plot('source','Source')
    roi.unzero_source(which=name)

    # ts map which shouls nearby 2FGL sources
    for source in new_sources:
        roi.zero_source(which=source.name)
    _plot('newsrc','newsrc')
    for source in new_sources:
        roi.unzero_source(which=source.name)

def counts_plots(roi, name, hypothesis, datadir, plotdir, size, pixelsize, **common_kwargs):
    """ Counts plots """
    emin, emax = get_full_energy_range(roi)
    
    extra='%s_%s_%sdeg_%sdeg' % (hypothesis,name,size,pixelsize)

    counts_kwargs = dict(size=size, **common_kwargs)
    title = 'Counts Residual for %s (%s)' % (name,hypothesis)
    print title
    roi.plot_counts_map(filename="%s/counts_residual_%s.png"%(plotdir,extra),
                        countsfile="%s/counts_residual_%s.fits"%(datadir,extra),
                        modelfile="%s/model_residual_%s.fits"%(datadir,extra),
                        pixelsize=pixelsize,
                        title=title,
                        **counts_kwargs)
    roi.zero_source(which=name)

    title = 'Counts Source for %s (%s)' % (name,hypothesis)
    print title
    roi.plot_counts_map(filename="%s/counts_source_%s.png"%(plotdir,extra),
                        countsfile="%s/counts_source_%s.fits"%(datadir,extra),
                        modelfile="%s/model_source_%s.fits"%(datadir,extra),
                        pixelsize=pixelsize,
                        title=title,
                        **counts_kwargs)
    roi.unzero_source(which=name)

    title = 'Slice for %s (%s)' % (name,hypothesis)
    print title
    roi.plot_slice(which=name,
                   pixelsize=pixelsize,
                   filename="%s/counts_slice_%s.png"%(plotdir,extra),
                   datafile='%s/counts_slice_%s.dat'%(datadir,extra),
                   title=title)

    title = 'Radial Integral for %s (%s)' % (name,hypothesis)
    print title
    roi.plot_radial_integral(which=name,
                             pixelsize=pixelsize,
                             filename="%s/radial_integral_%s.png"%(plotdir,extra),
                             datafile='%s/radial_integral_%s.dat'%(datadir,extra),
                             title=title)
    try:
        title = 'Spectra for %s (%s)' % (name,hypothesis)
        print title
        roi.plot_counts_spectra(filename="%s/spectra_%s_%s.png"%(plotdir,hypothesis, name),
                               title=title)
    except Exception, ex:
        print 'ERROR with plot_counts_spectra: ', ex
        traceback.print_exc(file=sys.stdout) 

def smooth_plots(roi, name, hypothesis, datadir, plotdir, size, kernel_rad, **common_kwargs):
    """ smoothed counts maps """
    emin, emax = get_full_energy_range(roi)

    extra='%s_%s_%sdeg_%sdeg' % (hypothesis, name, size,kernel_rad)

    smooth_kwargs = dict(which=name, 
                         override_center=roi.roi_dir,
                         size=size,
                         colorbar_radius=1, # most interesting within one degrees
                         kernel_rad=kernel_rad,
                         **common_kwargs)

    title = 'Source Map for %s (%s)' % (name,hypothesis)
    print title
    roi.plot_source(filename='%s/source_%s.png' % (plotdir, extra), 
                    title=title,
                    **smooth_kwargs)

    title = 'Sources Map for %s (%s)' % (name,hypothesis)
    print title
    roi.plot_sources(filename='%s/sources_%s.png' % (plotdir, extra), 
                     title=title,
                     **smooth_kwargs)

    if all_energy(emin,emax):
        title = 'Band Source Map for %s (%s)' % (name,hypothesis)
        print title
        ROISourceBandPlotter(roi, bin_edges=one_bin_per_dec(emin,emax), 
                             title=title,
                             **smooth_kwargs).show(filename='%s/band_source_%s.png' % (plotdir,extra))
        title = 'Band Sources Map for %s (%s)' % (name,hypothesis)
        print title
        ROISourcesBandPlotter(roi, bin_edges=one_bin_per_dec(emin,emax), 
                             title=title,
                              **smooth_kwargs).show(filename='%s/band_sources_%s.png' % (plotdir,extra))

def setup_plotting(roi, name, hypothesis, 
          pulsar_position, new_sources,
          datadir, plotdir):

    extra_overlay = lambda ax: overlay_on_plot(ax, pulsar_position=pulsar_position)

    # Override marker for new sources to be red stars
    override_kwargs = {source.name:dict(color='red',marker='*') for source in new_sources}

    common_kwargs = dict(extra_overlay=extra_overlay, 
                         overlay_kwargs=dict(override_kwargs=override_kwargs))

    args = (roi, name, hypothesis, datadir, plotdir)

    return args, common_kwargs

def plots(roi, name, hypothesis, 
          pulsar_position, new_sources,
          datadir, plotdir,
          size=5):

    print 'Making plots for hypothesis %s' % hypothesis

    args, common_kwargs = setup_plotting(roi, name, hypothesis, 
                                         pulsar_position, new_sources,
                                         datadir, plotdir)

    smooth_plots(*args, kernel_rad=0.1, size=size, **common_kwargs)
    counts_plots(*args, pixelsize=0.1, size=size, **common_kwargs)

    smooth_plots(*args, kernel_rad=0.25, size=size, **common_kwargs)
    counts_plots(*args, pixelsize=0.25, size=size, **common_kwargs)

def tsmaps(roi, name, hypothesis, 
          pulsar_position, new_sources,
          datadir, plotdir):

    print 'Making tsmaps for hypothesis %s' % hypothesis

    args, common_kwargs = setup_plotting(roi, name, hypothesis, 
                                         pulsar_position, new_sources,
                                         datadir, plotdir)

    for size in [5,10]:
        tsmap_plots(*args, tsmap_pixelsize=0.1, size=size, new_sources=new_sources, **common_kwargs)
