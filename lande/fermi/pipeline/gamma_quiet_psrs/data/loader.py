from os.path import expandvars, join

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
        return self.radiopsr_data[psr]['ft1']

    def get_ft2(self,psr):
        return self.radiopsr_data[psr]['ft2']

    def get_ltcube(self,psr):
        return self.radiopsr_data[psr]['ltcube']

    def get_binfile(self,psr, binsperdec):
        return join(self.radiopsr_data[psr]['binfiledir'], 'binned_%sbpd.fits' % binsperdec)

