# Not entirely sure why, but pyLikelihood
# will get a swig::stop_iteraiton error
# unless it is imported first.
from lande.fermi.likelihood.roi_gtlike import Gtlike

import traceback
import sys
import os
import pylab as P
import numpy as np

import yaml

from uw.like.roi_state import PointlikeState
from uw.pulsar.phase_range import PhaseRange
from uw.like.SpatialModels import Gaussian

from lande.utilities.tools import tolist

from lande.fermi.likelihood.fit import paranoid_gtlike_fit, fit_prefactor, fit_only_source, freeze_insignificant_to_catalog, freeze_bad_index_to_catalog
from lande.fermi.likelihood.save import source_dict, get_full_energy_range
from lande.fermi.likelihood.limits import GtlikePowerLawUpperLimit, GtlikeCutoffUpperLimit, PointlikePowerLawUpperLimit, PointlikeCutoffUpperLimit


from lande.fermi.likelihood.localize import GridLocalize, paranoid_localize
from lande.fermi.likelihood.cutoff import PointlikeCutoffTester, GtlikeCutoffTester
from lande.fermi.likelihood.bandfitter import GtlikeBandFitter
from lande.fermi.likelihood.printing import summary
from lande.fermi.pulsar.plotting import plot_phaseogram, plot_phase_vs_time
from lande.fermi.spectra.gtlike import GtlikeSED
from lande.fermi.spectra.pointlike import PointlikeSED
from lande.fermi.data.plotting import ROITSMapBandPlotter, ROISourceBandPlotter, ROISourcesBandPlotter
from lande.fermi.likelihood.free import freeze_far_away, unfreeze_far_away

close=lambda x,y: np.allclose(x,y, rtol=0, atol=1)
all_energy=lambda emin,emax: close([emin,emax],[1e2,10**5.5]) or close([emin,emax],[1e2,1e5])
high_energy=lambda emin,emax: close([emin,emax],[10**4,10**5.5])
higher_energy=lambda emin,emax: close([emin,emax],[10**4.5,10**5.5])

def one_bin_per_dec(emin,emax):
    assert close(emin,1e2) and (close(emax,1e5) or close(emax,10**5.5))
    if close(emax,1e5):
        return np.logspace(2,5,4)
    elif close(emax,10**5.5):
        return [1e2,1e3,1e4,10**5.5]

def two_bin_per_dec(emin,emax):
    assert close(emin,1e2) and (close(emax,1e5) or close(emax,10**5.5))
    if close(emax,1e5):
        return np.logspace(2,5,7)
    elif close(emax,10**5.5):
        return np.logspace(2,5.5,8)

def four_bin_per_dec(emin,emax):
    assert close(emin,1e2) and (close(emax,1e5) or close(emax,10**5.5))
    if close(emax,1e5):
        return np.logspace(2,5,13)
    elif close(emax,10**5.5):
        return np.logspace(2,5.5,15)

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
        roi.plot_tsmap(filename='%s/tsmap_%s_%s.png' % (plotdir,filename,extra), 
                       fitsfile='%s/tsmap_%s_%s.fits' % (datadir,filename,extra),
                       title='%s TS Map for %s (%s)' % (title,name,hypothesis),
                       **tsmap_kwargs)

        if all_energy(emin,emax):
            ROITSMapBandPlotter(roi,  
                                title='Band %s TS Map for %s (%s)' % (title,name,hypothesis),
                                bin_edges=one_bin_per_dec(emin,emax),
                                **tsmap_kwargs).show(filename='%s/band_tsmap_%s_%s.png' % (plotdir,filename,extra))

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
    roi.plot_counts_map(filename="%s/counts_residual_%s.png"%(plotdir,extra),
                        countsfile="%s/counts_residual_%s.fits"%(datadir,extra),
                        modelfile="%s/model_residual_%s.fits"%(datadir,extra),
                        pixelsize=pixelsize,
                        title='Counts Residual for %s (%s)' % (name,hypothesis),
                        **counts_kwargs)
    roi.zero_source(which=name)

    roi.plot_counts_map(filename="%s/counts_source_%s.png"%(plotdir,extra),
                        countsfile="%s/counts_source_%s.fits"%(datadir,extra),
                        modelfile="%s/model_source_%s.fits"%(datadir,extra),
                        pixelsize=pixelsize,
                        title='Counts Source for %s (%s)' % (name,hypothesis),
                        **counts_kwargs)
    roi.unzero_source(which=name)

    roi.plot_slice(which=name,
                   pixelsize=pixelsize,
                   filename="%s/counts_slice_%s.png"%(plotdir,extra),
                   datafile='%s/counts_slice_%s.dat'%(datadir,extra),
                   title='Slice for %s (%s)' % (name,hypothesis))

    roi.plot_radial_integral(which=name,
                             pixelsize=pixelsize,
                             filename="%s/radial_integral_%s.png"%(plotdir,extra),
                             datafile='%s/radial_integral_%s.dat'%(datadir,extra),
                             title='Radial Integral for %s (%s)' % (name,hypothesis))
    try:
        roi.plot_counts_spectra(filename="%s/spectra_%s_%s.png"%(plotdir,hypothesis, name),
                               title='Spectra for %s (%s)' % (name,hypothesis))
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

    roi.plot_source(filename='%s/source_%s.png' % (plotdir, extra), 
                    title='Source Map for %s (%s)' % (name,hypothesis),
                    **smooth_kwargs)
    roi.plot_sources(filename='%s/sources_%s.png' % (plotdir, extra), 
                     title='Sources Map for %s (%s)' % (name,hypothesis),
                     **smooth_kwargs)

    if all_energy(emin,emax):
        ROISourceBandPlotter(roi, bin_edges=one_bin_per_dec(emin,emax), 
                             title='Band Source Map for %s (%s)' % (name,hypothesis),
                             **smooth_kwargs).show(filename='%s/band_source_%s.png' % (plotdir,extra))
        ROISourcesBandPlotter(roi, bin_edges=one_bin_per_dec(emin,emax), 
                             title='Band Sources Map for %s (%s)' % (name,hypothesis),
                              **smooth_kwargs).show(filename='%s/band_sources_%s.png' % (plotdir,extra))


def plots(roi, name, hypothesis, 
          pulsar_position, new_sources,
          do_plots, do_tsmap,
          datadir='data', plotdir='plots'):

    print 'Making plots for hypothesis %s' % hypothesis

    extra_overlay = lambda ax: overlay_on_plot(ax, pulsar_position=pulsar_position)

    # Override marker for new sources to be red stars
    override_kwargs = {source.name:dict(color='red',marker='*') for source in new_sources}

    common_kwargs = dict(extra_overlay=extra_overlay, 
                         overlay_kwargs=dict(override_kwargs=override_kwargs))

    for dir in [datadir, plotdir]: 
        if not os.path.exists(dir): os.makedirs(dir)

    args = (roi, name, hypothesis, datadir, plotdir)

    if do_plots:
        for size in [5]:
            smooth_plots(*args, kernel_rad=0.1, size=size, **common_kwargs)
            counts_plots(*args, pixelsize=0.1, size=size, **common_kwargs)

            smooth_plots(*args, kernel_rad=0.25, size=size, **common_kwargs)
            counts_plots(*args, pixelsize=0.25, size=size, **common_kwargs)
    if do_tsmap:
        for size in [5,10]:
            tsmap_plots(*args, tsmap_pixelsize=0.1, size=size, new_sources=new_sources, **common_kwargs)

    roi.toRegion('%s/region_%s_%s.reg'%(datadir,hypothesis, name))

def pointlike_analysis(roi, name, hypothesis, max_free,
                       seddir='seds', datadir='data', 
                       localize=False,
                       fit_extension=False, 
                       cutoff=False,
                       model1=None,
                      ):
    """ emin + emax used for computing upper limits. """
    print 'Performing Pointlike analysis for %s' % hypothesis

    for dir in [seddir, datadir]: 
        if not os.path.exists(dir): os.makedirs(dir)

    print_summary = lambda: roi.print_summary(galactic=True, maxdist=10)
    print_summary()

    emin, emax = get_full_energy_range(roi)

    print roi

    def fit(just_prefactor=False, just_source=False, fit_bg_first=False):
        """ Convenience function incase fit fails. """
        try:
            if just_prefactor:
                fit_prefactor(roi, name) 
            elif just_source:
                fit_only_source(roi, name)
            elif fit_bg_first:
                roi.fit(fit_bg_first=True)
            else:
                roi.fit()
                # For some reason, one final fit seems to help with convergence and not getting negative TS values *shurgs*
                roi.fit() 
        except Exception, ex:
            print 'ERROR spectral fitting pointlike: ', ex
            traceback.print_exc(file=sys.stdout)
        print_summary()

    fit(just_prefactor=True)
    fit(fit_bg_first=True)
    fit() 

    frozen  = freeze_far_away(roi, roi.get_source(name).skydir, max_free)

    if localize:
        try:

            print 'About to Grid localize'
            grid=GridLocalize(roi,which=name,
                              update=True,
                              size=0.5, pixelsize=0.1)
            print_summary()

            paranoid_localize(roi, name, update=True)
        except Exception, ex:
            print 'ERROR localizing pointlike: ', ex
            traceback.print_exc(file=sys.stdout)

    if fit_extension:
        init_flux = roi.get_model(which=name).i_flux(emin,emax)
        try:
            before_state = PointlikeState(roi)

            roi.fit_extension(which=name)
            paranoid_localize(roi, name, update=True)

        except Exception, ex:
            print 'ERROR extension fitting pointlike: ', ex
            traceback.print_exc(file=sys.stdout)


    unfreeze_far_away(roi, frozen)

    fit()

    print 'Making pointlike SED'
    sed = PointlikeSED(roi, name, verbosity=True)
    sed.save('%s/sed_pointlike_%s_%s.yaml' % (seddir,hypothesis,name))
    sed.plot('%s/sed_pointlike_%s_%s.png' % (seddir,hypothesis,name)) 

    print_summary()

    p = source_dict(roi, name)

    pul = PointlikePowerLawUpperLimit(roi, name, emin=emin, emax=emax, cl=.95)
    p['powerlaw_upper_limit']=pul.todict()
    cul = PointlikeCutoffUpperLimit(roi, name, Index=1.7, Cutoff=3e3, b=1, cl=.95)
    p['cutoff_upper_limit']=cul.todict()

    if cutoff:
        try:
            tc = PointlikeCutoffTester(roi,name, model1=model1, verbosity=True)
            p['test_cutoff']=tc.todict()
            tc.plot(sed_results='%s/sed_pointlike_%s_%s.yaml' % (seddir,hypothesis,name),
                    filename='%s/test_cutoff_pointlike_%s_%s.png' % (plotdir,hypothesis,name))
        except Exception, ex:
            print 'ERROR plotting cutoff test:', ex
            traceback.print_exc(file=sys.stdout)


    roi.toXML(filename="%s/srcmodel_pointlike_%s_%s.xml"%(datadir, hypothesis, name))
 
    roi.save('roi_%s_%s.dat' % (hypothesis,name))

    return p


def gtlike_analysis(roi, name, hypothesis, max_free,
                    seddir='seds', datadir='data', plotdir='plots',
                    upper_limit=False, cutoff=False, 
                    model1=None):
    print 'Performing Gtlike crosscheck for %s' % hypothesis

    for dir in [seddir, datadir, plotdir]: 
        if not os.path.exists(dir): os.makedirs(dir)

    frozen  = freeze_far_away(roi, roi.get_source(name).skydir, max_free)
    gtlike=Gtlike(roi)
    unfreeze_far_away(roi, frozen)

    global like
    like=gtlike.like

    like.tol = 1e-1 # I found that the default tol '1e-3' would get the fitter stuck in infinite loops

    import pyLikelihood as pyLike
    like.setFitTolType(pyLike.ABSOLUTE)

    emin, emax = get_full_energy_range(like)

    print 'About to fit gtlike ROI'

    print summary(like, maxdist=10)

    paranoid_gtlike_fit(like, niter=3)

    print 'Done fiting gtlike ROI'
    print summary(like, maxdist=10)

    like.writeXml("%s/srcmodel_gtlike_%s_%s.xml"%(datadir, hypothesis, name))

    r=source_dict(like, name)

    if upper_limit:
        pul = GtlikePowerlawUpperLimit(like, name, emin=emin, emax=emax, cl=.95, delta_log_like_limits=10)
        r['powerlaw_upper_limit'] = pul.todict()
        cul = GtlikeCutoffUpperLimit(like, name, Index=1.7, Cutoff=3e3, b=1, cl=.95)
        r['cutoff_upper_limit'] = cul.todict()

    if all_energy(emin,emax):
        bf = GtlikeBandFitter(like, name, bin_edges=one_bin_per_dec(emin,emax), verbosity=True)
        bf.plot('%s/bandfit_gtlike_%s_%s.png' % (plotdir,hypothesis,name))
        bf.save('%s/bandfit_gtlike_%s_%s.yaml' % (datadir,hypothesis,name))

    def sed(kind,**kwargs):
        try:
            print 'Making %s SED' % kind
            sed = GtlikeSED(like, name, always_upper_limit=True, verbosity=True, **kwargs)
            sed.plot('%s/sed_gtlike_%s_%s.png' % (seddir,kind,name)) 
            sed.save('%s/sed_gtlike_%s_%s.yaml' % (seddir,kind,name))
        except Exception, ex:
            print 'ERROR computing SED:', ex
            traceback.print_exc(file=sys.stdout)

    if all_energy(emin,emax):
        sed('1bpd_%s' % hypothesis,bin_edges=one_bin_per_dec(emin,emax))
        sed('2bpd_%s' % hypothesis,bin_edges=two_bin_per_dec(emin,emax))
        sed('4bpd_%s' % hypothesis,bin_edges=four_bin_per_dec(emin,emax))
    elif high_energy(emin,emax):
        sed('4bpd_%s' % hypothesis,bin_edges=np.logspace(4,5.5,7))
        sed('2bpd_%s' % hypothesis,bin_edges=np.logspace(4,5.5,4))
    elif higher_energy(emin,emax):
        sed('4bpd_%s' % hypothesis,bin_edges=np.logspace(4.5,5.5,5))
        sed('2bpd_%s' % hypothesis,bin_edges=np.logspace(4.5,5.5,3))
    else:
        sed(hypothesis)

    if cutoff:
        try:
            tc = GtlikeCutoffTester(like,name, model1=model1, verbosity=True)
            r['test_cutoff']=todict()
            tc.plot(sed_results='%s/sed_gtlike_2bpd_%s_%s.yaml' % (seddir,hypothesis,name),
                    filename='%s/test_cutoff_gtlike_%s_%s.png' % (plotdir,hypothesis,name))
        except Exception, ex:
            print 'ERROR plotting cutoff test:', ex
            traceback.print_exc(file=sys.stdout)

    return r
    
def save_results(results, filename): 
    open(filename,'w').write(yaml.dump(tolist(results)))

def import_module(filename):
    """ import a python module from a pathname. """
    import imp
    return imp.load_source('module',filename)

