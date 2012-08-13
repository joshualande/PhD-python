from lande.fermi.likelihood.save import sourcedict
from lande.fermi.likelihood.printing import summary

def pointlike_counts(roi, name, plotdir,size):

    extra='%s_%sdeg' % (name,size)
    roi.plot_sources(filename='%s/sources_%s.png' % (plotdir,extra),
                     title='Sources Map for %s' % (name),
                     size=size)

def pointlike_tsmap(roi, name, plotdir,size,pixelsize):

    extra='%s_%sdeg' % (name,size)
    roi.plot_tsmap(filename='%s/tsmap_%s_%s.png' % (plotdir,'residual',extra), 
                   fitsfile='%s/tsmap_%s_%s.fits' % (plotdir,'residual',extra),
                   title='%s TS Map for %s' % ('Residual',name),
                   size=size,pixelsize=pixelsize)

def pointlike_plots(roi, name):
    pointlike_counts(roi, name, size=5)
    pointlike_tsmap(roi, name, size=5, pixelsize=0.1)

def pointlike_analysis(roi, name, plotdir):

    print_summary = lambda: roi.print_summary(galactic=True, maxdist=10)

    print_summary()
    roi.fit(use_gradient=False)
    print_summary()

    results  = sourcedict(roi, name)

    results['powerlaw_upper_limit'] = powerlaw_upper_limit(roi, name, emin=emin, emax=emax, cl=.95)

    roi.plot_sed(which=name, filename='%s/sed_pointlike_%s.png' % (plotdir,name), use_ergs=True)

    return results

def gtlike_analysis(roi, name, plotdir):
    gtlike = Gtlike(roi)
    like = gtlike

    print_summary = lambda: print summary(like, maxdist=10)

    print_summary()
    paranoid_gtlike_fit(like)
    print_summary()

    results=sourcedict(like, name)

    results['powerlaw_upper_limit'] = powerlaw_upper_limit(like, name, emin=emin, emax=emax, cl=.95, delta_log_like_limits=10)

    sed = GtlikeSED(like, name, always_upper_limit=True)
    sed.plot('%s/sed_gtlike_%s.png' % (plotdir,name))
    sed.save('%s/sed_gtlike_%s.yaml' % (plotdir,name))

    return results

