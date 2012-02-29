from os.path import join, exists
import pyfits
import numpy as np
import warnings
import collections
from tempfile import NamedTemporaryFile
import os

from uw.like.pointspec import SpectralAnalysis

class SpectralAnalysisCache(SpectralAnalysis):
    """ Acts just like SpectralAnalysis,
        but creates a cached ft1 file which is smaller. 
        
        This is helpful if you start with an all-sky ft1file
        but want to quickly make counts maps in pointlike. """

    def __init__(self, data_specification, cachedir, **kwargs):

        # n.b. bin the ft1 file before changing it
        super(SpectralAnalysisCache,self).__init__(data_specification, **kwargs)

        if not exists(cachedir):
            os.makedirs(cachedir)

        r = self.maxROI + 5
        l,b=self.roi_dir.l(),self.roi_dir.b()

        name = os.path.basename(data_specification.ft1files[0])
        name = name.replace('.fits','_cached_l_%g_b_%g_r_%g.fits' % (l,b,r))

        cachefile = join(cachedir,name)

        cache_ft1(ft1files = data_specification.ft1files, 
                  skydir = self.roi_dir,
                  radius = r,
                  cachefile=cachefile)

        self.ft1files = cachefile
        self.dataspec.ft1files = cachefile
        self.pixeldata.ft1files = cachefile


def cache_ft1(ft1files, cachefile, skydir, radius, emin=10, emax=1e6):
    """ This function runs gtselect to cache an ft1 file for easy use"""

    if exists(cachefile):
        try:
            # Sometimes the cachefile has problems with it
            # (the previous cache failed but still created
            # a file). So do some tests to make sure the
            # file looks good. If not delete and recreate it.

            # make warnigns raise exceptions
            warnings.filterwarnings('error')
            x=pyfits.open(cachefile)

            # test to see if GTI & Photons are in file
            x['GTI']
            x['EVENTS']

            if np.all(x['EVENTS'].data.field('TIME')==0):
                raise Exception("All TIMES in file are 0")
            if np.all(x['GTI'].data.field('START')==0):
                raise Exception("All GTIs in START are 0")
            if np.all(x['GTI'].data.field('STOP')==0):
                raise Exception("All GTIs in STOP are 0")

            if not x['EVENTS'].header.has_key('NDSKEYS'):
                raise Exception("EVENTS header does not contain NDSKEYS entry")

        except Exception, ex:
            print 'Error reading file %s.' % cachefile
            print 'Error: ',ex
            print 'Deleting and recaching file.'
            os.remove(cachefile)

        warnings.resetwarnings()

    if not exists(cachefile): 

        print 'Caching ft1 file. Saving to %s' % cachefile

        if len(ft1files) == 1:
            infile=ft1files[0]
        else:
            temp=NamedTemporaryFile(delete=True)
            temp.write('\n'.join(ft1files))
            temp.seek(0)
            infile='@%s' % temp.name


        import GtApp
        GtApp.GtApp("gtselect",'dataSubselector').run(
                infile=infile,
                outfile=cachefile,
                ra=skydir.ra(),
                dec=skydir.dec(),
                rad=radius,
                tmin=0, tmax=0,
                emin=emin,
                emax=emax,
                zmax=180)

