from os.path import expandvars

import numpy as np
import pyfits

from skymaps import SkyDir

from lande.utilities.save import loaddict

class RadioPSRLoader(object):

    def __init__(self, radiopsr_data, bigfile):

        self.radiopsr_data = loaddict(radiopsr_data)
        self.bigfile = pyfits.open(expandvars(bigfile))


    def _get_bigfile(self,psr):
        bf=self.bigfile['PULSARS_BIGFILE'].data
        index = np.where(bf['PSRJ'] == psr)[0]
        assert len(index)==1
        index=index[0]
        return bf[index]

    def get_skydir(self,psr):
        bf = self._get_bigfile(psr)
        ra = bf['RAJD']
        dec = bf['DECJD']
        return SkyDir(ra, dec)


    def get_ft1(self,psr):
        from glob import glob
        return glob("/afs/slac/g/glast/groups/catalog/P7_V4_SOURCE/pass7.3_pre_source_merit_*_pass7.4_source_z100_t90_cl0.fits")

    def get_ft2(self,psr):
        return "/afs/slac/g/glast/groups/catalog/P7_V4_SOURCE/ft2_2years.fits"

    def get_ltcube(self,psr):
        return "/afs/slac/g/glast/groups/catalog/P7_V4_SOURCE/ltcube_24m_pass7.4_source_z100_t90_cl0.fits"

    def get_binfile(self,psr):
        return '/nfs/slac/g/ki/ki03/lande/fermi/radiopsrs/lat_data/temp/binned.fits'


