from os.path import expandvars

class PulsarCatalogLoader(object):
    def __init__(self,bigfile_filename,off_peak_auxiliary_filename):
        self.bigfile_fits = pyfits.open(expandvars(bigfile_filename))['PULSARS_BIGFILE']
        self.off_peak_fits = pyfits.open(expandvars(off_peak_auxiliary_filename))['OFF_PEAK']

    def get_off_peak_classification(self, psr):

        off_peak = self._get_off_peak(psr)
        classification = off_peak['Classification']
        return classification

    def get_off_peak_psrlist(self):
        return np.char.strip(self.off_peak_fits.data['PSR']).tolist()

    def _get_off_peak(self,psr):
        return self.off_peak_fits.data[self.off_peak_fits.data['PSR'] == psr][0]

    def _get_bigfile(self,psr):
        return self.bigfile_fits.data[self.bigfile_fits.data['PSRJ'] == self.name_mapper(psr)][0]

    def name_mapper(self,psr):
        """ Convert from off-peak pulsar naming convention to BigFile naming convention. """
        return psr.replace('PSRJ','J')

    def get_edot(self, psr):
        """ Return edot in units of erg/s """
        bigfile = self._get_bigfile(psr)
        Edot = bigfile['Edot']
        return Edot

    def get_off_peak_eflux(self, psr):

        off_peak = self._get_off_peak(psr)

        eflux = off_peak['EFlux']
        eflux_error = off_peak['EFlux_error']
        eflux_ul = off_peak['PowerLaw_EFlux_UL']

        return eflux, eflux_error, eflux_ul

    def get_luminosity(self, psr):
        """ Returns luminsoity in units of erg/s """

        eflux, eflux_error, eflux_ul = self.get_off_peak_eflux(psr)
        
        bigfile = self._get_bigfile(psr)

        d1 = bigfile['DPSR_1']
        d2 = bigfile['DPSR_2']

        luminosity, luminosity_lower_error, luminosity_upper_error, luminosity_ul, luminosity_significant = None, None, None, None, None

        classification = self.get_off_peak_classification(psr)

        convert = float(units.kpc**2/units.cm**2)

        if d1[0] == '<':
            # luminosity is upper limit
            assert np.isnan(d2) # no 2nd estimate when upper limit

            dist_ul = float(d1[1:])

            if classification == 'Upper_Limit':
                print 'distance UL + flux UL'
                luminosity_ul = eflux_ul*4*np.pi*dist_ul**2 * convert
                luminosity_significant = False
            else:
                print 'distance UL + flux detection'
                luminosity_ul = (eflux + eflux_error)*4*np.pi*dist_ul**2 * convert
                luminosity_significant = False

        else:
            if np.isnan(d2):
                # no 2nd distance

                print 'distance detection + flux detection'
                dist = float(d1)
                d1_lower_error,d1_upper_error=eval(bigfile['e_DPSR_1_stat'])

                if classification == 'Upper_Limit':
                    luminosity_ul = eflux_ul*4*np.pi*(dist + d1_upper_error) * convert
                    luminosity_significant = False
                else:
                    
                    luminosity = eflux*4*np.pi*dist**2 * convert
                    luminosity_lower_error = luminosity_upper_error = np.sqrt(
                        ((eflux_error)*4*np.pi*dist**2)**2 + 
                        (eflux*4*np.pi*2*dist*d1_lower_error)**2
                    ) * convert
                    luminosity_upper_error = luminosity_upper_error = np.sqrt(
                        ((eflux_error)*4*np.pi*dist**2)**2 + 
                        (eflux*4*np.pi*2*dist*d1_upper_error)**2
                    ) * convert

                    luminosity_significant = True
            else:
                raise Exception('two distances. Not sure what to do')


        return luminosity, luminosity_lower_error, luminosity_upper_error, luminosity_ul, luminosity_significant
