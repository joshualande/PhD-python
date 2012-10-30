import traceback
import sys
import os

from uw.like.roi_state import PointlikeState

from lande.fermi.likelihood.fit import fit_prefactor, fit_only_source
from lande.fermi.likelihood.save import source_dict, get_full_energy_range
from lande.fermi.likelihood.limits import PointlikePowerLawUpperLimit

from lande.fermi.likelihood.localize import GridLocalize, paranoid_localize,MinuitLocalizer
from lande.fermi.spectra.pointlike import PointlikeSED
from lande.fermi.likelihood.free import freeze_far_away, unfreeze_far_away

def pointlike_analysis(roi, name, hypothesis, dirdict,
                       localize=False,
                       fit_extension=False, 
                      ):
    print 'Performing Pointlike analysis for %s' % hypothesis

    print_summary = lambda: roi.print_summary(galactic=True, maxdist=10)
    print_summary()

    def fit(just_prefactor=False, fit_bg_first=False):
        """ Convenience function incase fit fails. """
        try:
            if just_prefactor:
                fit_prefactor(roi, name) 
            else:
                roi.fit(fit_bg_first=fit_bg_first)
        except Exception, ex:
            print 'ERROR spectral fitting pointlike for hypothesis %s:' % hypothesis, ex
            traceback.print_exc(file=sys.stdout)
        print_summary()

    fit(just_prefactor=True)
    fit(fit_bg_first=True)
    fit() 

    if localize:
        paranoid_localize(roi, name, verbosity=4)

    if fit_extension:
        roi.fit_extension(which=name)
        paranoid_localize(roi, name)

    fit()

    print 'Making pointlike SED for hypothesis %s' % hypothesis
    sed = PointlikeSED(roi, name, verbosity=4)
    sed.save('%s/sed_pointlike_4bpd_%s_%s.yaml' % (dirdict['seds'],hypothesis,name))
    sed.plot('%s/sed_pointlike_4bpd_%s_%s.png' % (dirdict['seds'],hypothesis,name)) 

    print_summary()

    p = source_dict(roi, name)

    pul = PointlikePowerLawUpperLimit(roi, name, cl=.95, verbosity=4)
    p['powerlaw_upper_limit']=pul.todict()

    roi.toXML(filename="%s/srcmodel_pointlike_%s_%s.xml"%(dirdict['data'], hypothesis, name))
 
    roi.save('roi_%s_%s.dat' % (hypothesis,name))

    return p
