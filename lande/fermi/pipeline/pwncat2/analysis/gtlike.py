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
