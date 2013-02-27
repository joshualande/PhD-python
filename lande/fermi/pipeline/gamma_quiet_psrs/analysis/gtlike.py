# Not entirely sure why, but pyLikelihood
# will get a swig::stop_iteraiton error
# unless it is imported first.
from lande.fermi.likelihood.roi_gtlike import Gtlike

import numpy as np

from lande.fermi.likelihood.fit import paranoid_gtlike_fit
from lande.fermi.likelihood.save import source_dict
from lande.fermi.likelihood.limits import GtlikePowerLawUpperLimit

from lande.fermi.likelihood.printing import summary
from lande.fermi.spectra.gtlike import GtlikeSED

def gtlike_analysis(pipeline, roi, name, hypothesis, upper_limit):
    print 'Performing Gtlike crosscheck for %s' % hypothesis

    gtlike=Gtlike(roi, savedir='savedir' if pipeline.cachedata else None)
    like=gtlike.like

    print 'About to fit gtlike ROI'

    print summary(like, maxdist=10)

    paranoid_gtlike_fit(like, verbosity=4)

    print 'Done fiting gtlike ROI'
    print summary(like, maxdist=10)

    like.writeXml("%s/srcmodel_gtlike_%s_%s.xml"%(pipeline.dirdict['data'], hypothesis, name))

    r=source_dict(like, name)

    upper_limit_kwargs=dict()

    if upper_limit:
        pul = GtlikePowerLawUpperLimit(like, name, cl=.95, verbosity=4)
        r['powerlaw_upper_limit'] = pul.todict()

    def sed(kind,**kwargs):
        print 'Making %s SED' % kind
        s = GtlikeSED(like, name, always_upper_limit=True, 
                        verbosity=4, 
                        upper_limit_kwargs=upper_limit_kwargs,
                        **kwargs)
        s.plot('%s/sed_gtlike_%s_%s.png' % (pipeline.dirdict['seds'],kind,name)) 
        s.save('%s/sed_gtlike_%s_%s.yaml' % (pipeline.dirdict['seds'],kind,name))

    sed('1bpd_%s' % hypothesis,bin_edges=[10**2,10**3,10**4,10**5.5])
    sed('2bpd_%s' % hypothesis,bin_edges=np.logspace(2,5.5,8))
    if not pipeline.fast:
        sed('4bpd_%s' % hypothesis,bin_edges=np.logspace(2,5.5,15))

    return r
