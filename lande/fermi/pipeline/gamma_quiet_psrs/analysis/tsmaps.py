from lande.fermi.data.plotting import ROITSMapBandPlotter

def tsmap_plots(roi, name, hypothesis, dirdict, size, tsmap_pixelsize=0.1):
    """ TS maps """

    extra='%s_%s_%sdeg' % (hypothesis,name,size)

    tsmap_kwargs = dict(size=size, pixelsize=tsmap_pixelsize)

    def _plot(filename,title):
        print 'Making Band %s TS Map (size=%s, pixelsize=%s)' % (title,size,tsmap_pixelsize)
        roi.plot_tsmap(filename='%s/tsmap_%s_%s.png' % (dirdict['plots'],filename,extra), 
                       fitsfile='%s/tsmap_%s_%s.fits' % (dirdict['plots'],filename,extra),
                       title='%s TS Map for %s (%s)' % (title,name,hypothesis),
                       **tsmap_kwargs)

        print 'Making Band %s TS Map (size=%s, pixelsize=%s)' % (title,size,tsmap_pixelsize)
        ROITSMapBandPlotter(roi,  
                            title='Band %s TS Map for %s (%s)' % (title,name,hypothesis),
                            bin_edges=[1e3,1e4,10**5.5],
                            **tsmap_kwargs).show(filename='%s/band_tsmap_%s_%s.png' % (dirdict['plots'],filename,extra))

    # reisidual ts map
    _plot('residual','Residual')

    # source ts map
    roi.zero_source(which=name)
    _plot('source','Source')
    roi.unzero_source(which=name)
