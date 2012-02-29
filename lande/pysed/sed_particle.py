""" This code defines various particle spectra.

    Author: Joshua Lande <joshualande@gmail.com>
"""
from . sed_spectrum import Spectrum
import numpy as np

from . sed_integrate import logsimps
from . import sed_config
from . import sed_units as u

class ParticleSpectrum(Spectrum):
    """ This class to represents a spectrum of particles with total energy total_energy. 
    
        The __call__ function returns dn/de, the number of particles
        per unit energy.  This defaults to be in units of 1/erg, but
        can be set by the parameter units_string.

        Note that the total energy must have units of energy^2*units_string.
        By default, units_string has units of 1/erg so total_energy must have
        energy units. If units_string is instead 1/erg/second, then total_energy
        must be the total energy per unit time.
    """

    vectorized = True

    def __init__(self,total_energy, emin, emax, units_string='1/erg', *args, **kwargs):
        """ Normalize total energy output. """
        self.emin = float(emin/u.erg)
        self.emax = float(emax/u.erg)
        self._units_string = units_string
        self.norm = 1

        self.init(*args,**kwargs)

        self.norm=float(total_energy/self.integrate(units=True, e_weight=1))

    def integrate(self, units=True, e_weight=0):
        integral=logsimps(lambda e: e**(e_weight)*self(e, units=False),
                          self.emin,self.emax,per_decade=sed_config.PER_DECADE)
        return integral*(u.erg**(e_weight+1)*self.units() if units else 1)

    def units_string(self): return self._units_string

class PowerLaw(ParticleSpectrum):
    """ A simple powerlaw spectral model defined
        as dN/dE = N0*(E/E_0)**-gamma

        >>> p = PowerLaw(total_energy = 2e48*u.erg, index=2.6,
        ...              emin=1e-6*u.eV,emax=1e14*u.eV)
        >>> print u.repr(p.integrate(e_weight=1,units=True),'erg')
        2e+48 erg
    """

    def init(self, index, e_scale=u.GeV):
        self.index = index
        self.e_scale = float(e_scale/u.erg)

    def _spectrum(self, energy):
        """ Returns number of particles per unit energy [1/erg]. """
        return self.norm*(energy/self.e_scale)**(-self.index)


class PowerLawCutoff(ParticleSpectrum):

    def init(self, index, e_cutoff, e_scale=u.GeV):

        self.index = index
        self.e_cutoff = float(e_cutoff/u.erg)
        self.e_break = float(e_break/u.erg)
        self.e_scale = float(e_scale/u.erg)

    def _spectrum(self, energy):
        """ Returns number of particles per unit energy [1/erg]. """
        return self.norm*(energy/self.e_scale)**(-self.index)*np.exp(-energy/self.e_cutoff)


class SmoothBrokenPowerLaw(ParticleSpectrum):
    """ A smoothed broken power-law particle distribution.

        This formula is taken from the fermi-LAT publication
        on W51C: http://arxiv.org/abs/0910.0908

        but the hardcoded value 2 from the paper is settable
        as the parameter beta.
    """
    def init(self, index1, index2, e_break, e_scale, beta):
        self.index1 = index1
        self.index2 = index2
        self.e_break = float(e_break/u.erg)                                                                                                                               
        self.e_scale = float(e_scale/u.erg)                                                                                                                                 
        self.beta = beta

    def _spectrum(self, energy):
        """ Returns number of particles per unit energy [1/erg]. """
        return self.norm*(energy/self.e_scale)**(-self.index1)*(1 + (energy/self.e_break)**self.beta)**(-(self.index2-self.index1)/self.beta)

class BrokenPowerLawCutoff(ParticleSpectrum):

    def init(self, index1, index2, e_cutoff, e_break):
        self.index1 = index1
        self.index2 = index2
        self.e_cutoff = float(e_cutoff/u.erg)
        self.e_break = float(e_break/u.erg)
        self.e_scale = float(u.eV/u.erg)

    def _spectrum(self, energy):
        """ Return number of particles per unit energy [1/erg]. """
        return self.norm*((energy/self.e_scale)**-self.index1)/(1.+(energy/self.e_break)**(-self.index1+self.index2))*np.exp(-energy/self.e_cutoff)

if __name__ == "__main__":
    import doctest
    doctest.testmod()

