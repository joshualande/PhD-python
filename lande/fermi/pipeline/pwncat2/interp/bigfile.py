from os.path import expandvars

import numpy as np
import pyfits

from lande.pysed import units

class PulsarCatalogLoader(object):
    def __init__(self,bigfile_filename,off_peak_auxiliary_filename=None):
        self.bigfile_fits = pyfits.open(expandvars(bigfile_filename))['PULSARS_BIGFILE']

        if off_peak_auxiliary_filename is not None:
            self.off_peak_fits = pyfits.open(expandvars(off_peak_auxiliary_filename))['OFF_PEAK']

    def get_pulsar_classification(self, psr):
        bigfile = self._get_bigfile(psr)
        code=bigfile['PSR_Code']
        return code

    def get_bigfile_psrlist(self):
        return np.char.strip(self.bigfile_fits.data['PSRJ'])

    def get_off_peak_classification(self, psr):

        off_peak = self._get_off_peak(psr)
        classification = off_peak['Classification_OP']
        return classification

    def get_off_peak_psrlist(self):
        assert self.off_peak_fits is not None
        return np.char.strip(self.off_peak_fits.data['PSR'])

    def _get_off_peak(self,psr):
        assert self.off_peak_fits is not None
        return self.off_peak_fits.data[self.off_peak_fits.data['PSR'] == psr][0]

    def _get_bigfile(self,psr):
        return self.bigfile_fits.data[self.bigfile_fits.data['PSRJ'] == psr][0]

    def get_edot(self, psr):
        """ Return edot in units of erg/s """
        bigfile = self._get_bigfile(psr)
        Edot = bigfile['Edot']
        return Edot

    def get_off_peak_eflux(self, psr):

        off_peak = self._get_off_peak(psr)

        eflux = off_peak['EFlux_OP']
        eflux_error = off_peak['Unc_EFlux_OP']
        eflux_ul = off_peak['PowerLaw_EFlux_UL_OP']

        return eflux, eflux_error, eflux_ul

    def get_age(self, psr):
        bigfile = self._get_bigfile(psr)
        return bigfile['Age']

    def get_distance(self, psr):
        bigfile = self._get_bigfile(psr)

        d1 = bigfile['DPSR_1']
        d2 = bigfile['DPSR_2']
        assert np.isnan(d2) # no 2nd estimate when upper limit
        return d1

    def get_luminosity(self, psr):
        """ Returns luminsoity in units of erg/s """

        eflux, eflux_error, eflux_ul = self.get_off_peak_eflux(psr)
        
        bigfile = self._get_bigfile(psr)

        d1 = bigfile['DPSR_1']
        d2 = bigfile['DPSR_2']

        luminosity, luminosity_error_statistical, luminosity_lower_error_systematic, luminosity_upper_error_systematic, luminosity_ul, luminosity_significant = None, None, None, None, None, None

        classification = self.get_off_peak_classification(psr)

        convert = float(units.kpc**2/units.cm**2)

        if d1[0] == '<':
            # luminosity is upper limit
            assert np.isnan(d2) # no 2nd estimate when upper limit

            dist_ul = float(d1[1:])

            if classification == 'L':
                print 'distance UL + flux UL'
                luminosity_ul = eflux_ul*4*np.pi*dist_ul**2 * convert
                luminosity_significant = False
            elif classification in ['W','U','M','M*']:
                print 'distance UL + flux detection'
                luminosity_ul = (eflux + eflux_error)*4*np.pi*dist_ul**2 * convert
                luminosity_significant = False
            else:
                raise Exception("...")

        else:
            if np.isnan(d2):
                # no 2nd distance

                print 'distance detection + flux detection'
                dist = float(d1)
                d1_lower_error,d1_upper_error=eval(bigfile['e_DPSR_1_stat'])

                if classification == 'L':
                    luminosity_ul = eflux_ul*4*np.pi*(dist + d1_upper_error)**2 * convert
                    luminosity_significant = False
                elif classification in ['W','U','M','M*']:
                    luminosity = eflux*4*np.pi*dist**2 * convert
                    luminosity_error_statistical = eflux_error*4*np.pi*dist**2 * convert

                    #josh_method=False;xian_method=True
                    josh_method=True;xian_method=False

                    if josh_method:
                        print 'USING josh_method'
                    if xian_method:
                        print 'USING xian_method'

                    assert (josh_method and not xian_method) or (xian_method and not josh_method)

                    if josh_method:
                        luminosity_lower_error_systematic = luminosity - eflux*4*np.pi*(dist-d1_lower_error)**2 * convert
                        luminosity_upper_error_systematic = eflux*4*np.pi*(dist+d1_upper_error)**2 * convert - luminosity

                    if xian_method:
                        lum = lambda _f, _d: _f*4*np.pi*_d**2*convert

                        def xians_method(dist_err, npts=1000):
                            d_samples = np.random.normal(dist,dist_err,npts)
                            all_luminosities=lum(eflux,d_samples)
                            luminosity_error = np.std(all_luminosities)
                            return luminosity_error

                        luminosity_lower_error_systematic = xians_method(d1_lower_error)
                        luminosity_upper_error_systematic = xians_method(d1_upper_error)
                        if luminosity_lower_error_systematic > luminosity:
                            print 'Hitting weird edge case for psr=%s' % psr
                            luminosity_lower_error_systematic = 0.9*luminosity

                    luminosity_significant = True

                    if psr in [ 'J0205+6449','J0534+2200']:
                        print psr
                        print ' * -> eflux',eflux
                        print ' * -> eflux_error',eflux_error,eflux_error/eflux
                        print ' * -> dist',dist
                        print ' * -> d1_lower_error',d1_lower_error,d1_lower_error/dist
                        print ' * -> d1_upper_error',d1_upper_error,d1_upper_error/dist
                        print ' * -> luminosity',luminosity
                        print ' * -> luminosity_error_statistical',luminosity_error_statistical,luminosity_error_statistical/luminosity
                        print ' * -> luminosity_lower_error_systematic',luminosity_lower_error_systematic,luminosity_lower_error_systematic/luminosity
                        print ' * -> luminosity_upper_error_systematic',luminosity_upper_error_systematic,luminosity_upper_error_systematic/luminosity
                else:
                    raise Exception("...")
            else:
                raise Exception('two distances. Not sure what to do')


        return luminosity, luminosity_error_statistical, luminosity_lower_error_systematic, luminosity_upper_error_systematic, luminosity_ul, luminosity_significant
