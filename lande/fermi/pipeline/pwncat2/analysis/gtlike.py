# Not entirely sure why, but pyLikelihood
# will get a swig::stop_iteraiton error
# unless it is imported first.
from lande.fermi.likelihood.roi_gtlike import Gtlike

import traceback
import sys
import os
import numpy as np

import yaml

from lande.fermi.likelihood.fit import paranoid_gtlike_fit
from lande.fermi.likelihood.save import source_dict, get_full_energy_range
from lande.fermi.likelihood.limits import GtlikePowerLawUpperLimit, GtlikeCutoffUpperLimit


from lande.fermi.likelihood.cutoff import GtlikeCutoffTester
from lande.fermi.likelihood.bandfitter import GtlikeBandFitter
from lande.fermi.likelihood.printing import summary
from lande.fermi.spectra.gtlike import GtlikeSED
from lande.fermi.likelihood.free import freeze_far_away, unfreeze_far_away

from . binning import all_energy, high_energy, higher_energy, one_bin_per_dec, two_bin_per_dec, four_bin_per_dec

    


def gtlike_analysis(roi, name, hypothesis, max_free,
                    seddir, datadir, plotdir,
                    upper_limit=False, cutoff=False, 
                    model1=None):
    print 'Performing Gtlike crosscheck for %s' % hypothesis

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

    paranoid_gtlike_fit(like, niter=3, verbosity=4)

    print 'Done fiting gtlike ROI'
    print summary(like, maxdist=10)

    like.writeXml("%s/srcmodel_gtlike_%s_%s.xml"%(datadir, hypothesis, name))

    r=source_dict(like, name)

    #upper_limit_kwargs=dict(delta_log_like_limits=10)
    upper_limit_kwargs=dict()

    if upper_limit:
        pul = GtlikePowerLawUpperLimit(like, name, emin=emin, emax=emax, cl=.95,
                                       upper_limit_kwargs=upper_limit_kwargs,
                                       verbosity=4)
        r['powerlaw_upper_limit'] = pul.todict()
        cul = GtlikeCutoffUpperLimit(like, name, Index=1.7, Cutoff=3e3, b=1, cl=.95,
                                     upper_limit_kwargs=upper_limit_kwargs,
                                     override_model=model1,
                                     verbosity=4)
        r['cutoff_upper_limit'] = cul.todict()

    if all_energy(emin,emax):
        try:
            bf = GtlikeBandFitter(like, name, bin_edges=one_bin_per_dec(emin,emax), 
                                  upper_limit_kwargs=upper_limit_kwargs,
                                  verbosity=4)
            bf.plot('%s/bandfit_gtlike_%s_%s.png' % (plotdir,hypothesis,name))
            bf.save('%s/bandfit_gtlike_%s_%s.yaml' % (datadir,hypothesis,name))
        except Exception, ex:
            print 'ERROR computing bandfit:', ex
            traceback.print_exc(file=sys.stdout)

    def sed(kind,**kwargs):
        try:
            print 'Making %s SED' % kind
            sed = GtlikeSED(like, name, always_upper_limit=True, 
                            verbosity=4, 
                            upper_limit_kwargs=upper_limit_kwargs,
                            **kwargs)
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
            tc = GtlikeCutoffTester(like,name, model1=model1, verbosity=4)
            r['test_cutoff']=tc.todict()
            tc.plot(sed_results='%s/sed_gtlike_2bpd_%s_%s.yaml' % (seddir,hypothesis,name),
                    filename='%s/test_cutoff_gtlike_%s_%s.png' % (plotdir,hypothesis,name))
        except Exception, ex:
            print 'ERROR plotting cutoff test:', ex
            traceback.print_exc(file=sys.stdout)

    return r
