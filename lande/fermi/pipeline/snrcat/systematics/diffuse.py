import numpy as np

from skymaps import DiffuseFunction

from uw.utilities import keyword_options

from lande.fermi.data.gullidiffuse import get_gulli_diffuse
from lande.fermi.likelihood.diffuse import is_significant, merge_diffuse
from lande.fermi.likelihood.save import get_background



class GulliDiffuseSetup(object):

    defaults = (
        ('verbosity', False, 'Make lots of noise'),
        ('fraction', 0.03, 'Fraction to compare model to observed counts'),
        ('mergefile', None, 'name of merged file'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, roi, diffuse_kwargs, **kwargs):
        """ diffuse_kwargs passed into get_gulli_diffuse
        """
        keyword_options.process(self, kwargs)
        self.roi = roi
        self.diffuse_kwargs = diffuse_kwargs

    def update_roi(self):
        roi = self.roi

        if self.verbosity:
            print "Replacing current diffuse file with Gulli's diffuse files"

        if self.verbosity:
            print 'Deleting previous background sources:'
        for source in get_background(roi):
            if self.verbosity:
                print ' .. Deleting source %s' % source
            roi.del_source(source)

        if self.verbosity:
            print 'Getting Gulli Diffuse:'
        diffuse_sources = get_gulli_diffuse(verbosity=self.verbosity, **self.diffuse_kwargs)

        if self.verbosity:
            print 'Adding Gulli Diffuse files to ROI:'
        for source in diffuse_sources:
            if self.verbosity:
                print ' .. adding file %s to the ROI' % source.name
            roi.add_source(source)

        galatic_sources = [i for i in diffuse_sources if isinstance(i.dmodel[0],DiffuseFunction)]
        assert len(galatic_sources) == len(diffuse_sources) - 1, "Exactly one isotropic diffue source"

        free_galatic_sources = [i for i in diffuse_sources if np.any(i.smodel.free)]


        if self.verbosity:
            print 'Testing significance of free Galactic sources (to see if they need to be merged):'
        insignificant = [i for i in free_galatic_sources if not is_significant(roi, i.name, self.fraction, self.verbosity)]

        if len(insignificant) == 0:
            print 'No insignificant sources, returning!'
            return
        elif len(insignificant) == 0:
            print 'Only one insignificant source (%s). Freezing it and returning' % insignificant[0].name
            roi.modify(which=insignificant[0].name, free=False)
            return
        else:
            if self.verbosity:
                print 'Multiple insignificant sources in ROI. Merging them!'
                print 'Insignificant free Galactic diffuse sources are:'
                for i in insignificant:
                    print ' .. %s' % i.name

            if self.verbosity:
                print 'Deleting (and the merging) the insignificant free Galactic sources'
            for source in insignificant:
                if self.verbosity:
                    print ' .. Deleting (and then merging) %s' % source.name
                    roi.del_source(source.name)

            merged = merge_diffuse(insignificant, mergefile=self.mergefile, verbosity=self.verbosity)
            if self.verbosity:
                print 'Adding merged Galactic diffuse source %s to ROI:' % (merged.name)
            roi.add_source(merged)

            if self.verbosity:
                print 'Testing if the merged Galactic diffuse source is significant'
            if is_significant(roi, merged.name, self.fraction, self.verbosity):
                if self.verbosity:
                    print 'Keeping free the merged Galactic diffuse sources!'
            else:
                if self.verbosity:
                    print 'Freezing the insignificant merged Galactic diffuse sources!'
                roi.modify(which=merged.name, free=False)

    def __del__(self):
        pass
