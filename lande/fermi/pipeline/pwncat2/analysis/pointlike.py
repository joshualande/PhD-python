import traceback
import sys
import os

from uw.like.roi_state import PointlikeState

from lande.fermi.likelihood.fit import fit_prefactor, fit_only_source
from lande.fermi.likelihood.save import source_dict, get_full_energy_range
from lande.fermi.likelihood.limits import PointlikePowerLawUpperLimit, PointlikeCutoffUpperLimit


from lande.fermi.likelihood.localize import GridLocalize, paranoid_localize
from lande.fermi.likelihood.cutoff import PointlikeCutoffTester
from lande.fermi.spectra.pointlike import PointlikeSED
from lande.fermi.likelihood.free import freeze_far_away, unfreeze_far_away

def pointlike_analysis(roi, name, hypothesis, max_free,
                       seddir, datadir, plotdir,
                       localize=False,
                       fit_extension=False, 
                       cutoff=False,
                       model1=None,
                      ):
    """ emin + emax used for computing upper limits. """
    print 'Performing Pointlike analysis for %s' % hypothesis

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
    sed = PointlikeSED(roi, name, verbosity=4)
    sed.save('%s/sed_pointlike_4bpd_%s_%s.yaml' % (seddir,hypothesis,name))
    sed.plot('%s/sed_pointlike_4bpd_%s_%s.png' % (seddir,hypothesis,name)) 

    print_summary()

    p = source_dict(roi, name)

    pul = PointlikePowerLawUpperLimit(roi, name, emin=emin, emax=emax, cl=.95, verbosity=4)
    p['powerlaw_upper_limit']=pul.todict()
    cul = PointlikeCutoffUpperLimit(roi, name, Index=1.7, Cutoff=3e3, b=1, cl=.95, verbosity=4)
    p['cutoff_upper_limit']=cul.todict()

    if cutoff:
        try:
            tc = PointlikeCutoffTester(roi,name, model1=model1, verbosity=4)
            p['test_cutoff']=tc.todict()
            tc.plot(sed_results='%s/sed_pointlike_%s_%s.yaml' % (seddir,hypothesis,name),
                    filename='%s/test_cutoff_pointlike_%s_%s.png' % (plotdir,hypothesis,name))
        except Exception, ex:
            print 'ERROR plotting cutoff test:', ex
            traceback.print_exc(file=sys.stdout)


    roi.toXML(filename="%s/srcmodel_pointlike_%s_%s.xml"%(datadir, hypothesis, name))
 
    roi.save('roi_%s_%s.dat' % (hypothesis,name))

    return p
