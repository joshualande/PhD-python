""" sed_thermal.py

    Author: Joshua Lande <joshualande@gmail.com>
"""
import numpy as np
from scipy import integrate

from . sed_integrate import logsimps
from . sed_spectrum import Spectrum
from . import sed_config
from . import sed_units as u

class ThermalSpectrum(Spectrum):

    vectorized = True

    def __init__(self, energy_density, kT=None, T=None):
        """ A thermal spectrum has the sameself):
            spectral shape as the blackbody
            spectrum but has an arbitrarily 
            normalizable energy density.

            The thermal spectrum is

            n(E) = 15*U/(pi*kT)^4*E^2/(exp(E/kT)-1)

            where 
              * n(E) is the number of photons per unit energy per unit volume,
              * U is the total energy per unit volume.
              * kT is the temperature of the photons

            This formula is equation 33 from Sturner et al 1997
            http://iopscience.iop.org/0004-637X/490/2/619/pdf/35841.pdf
        
            Input can be either 'kT' in energy units or
            'T' in temperature units.

            For example, in XXX et al, the infrared photon
            field has temperature kT=3e-3 eV and energy
            density U=0.9 eV/cm^3

                >>> infrared=ThermalSpectrum(kT=3e-3*u.eV, energy_density=0.9*u.eV/u.cm**3)

            To convince yourself that this code correctly normalized
            the spectrum, you can explicity integrate E*dN/dE = total energy per unit volume:

                >>> print u.repr(infrared.integrate(units=True,e_weight=1),'eV/cm^3','%.2f')
                0.90 eV/cm^3
             """
        if kT is not None: kT = kT
        elif T is not None: kT = u.boltzmann*kwargs.pop('T')
        else: raise Exception("kT or k must be passed to ThermalSpectrum")

        self.kT = float(kT/u.erg)

        # function is essentially 0 outside of this energy range.
        self.emin=1e-4*self.kT
        self.emax=1e2*self.kT

        # equation 33 in Sturner et al 1997
        # Note, prefactor*E^2/(exp(E/kT)-1) has units
        # of photons/energy/volume, so prefactor has units
        # of photons/energy^3/volume.
        self.pref = 15*energy_density/(np.pi*kT)**4
        self.pref = float(self.pref/(u.erg**-3*u.cm**-3))

    @staticmethod
    def occupation_number(x):
        """ This is equation 1.49 in R&L. """
        return 1/(np.exp(x)-1)

    def _spectrum(self, energy):
        """ Return the energy density in units of [1/erg/cm^-3]."""
        return self.pref*energy**2*self.occupation_number(energy/self.kT)

    @staticmethod                                                                                                                                                           
    def units_string(): return '1/erg/cm^3'

    def integrate(self, units=True, e_weight=0):
        """ Integrate the thermal spectrum from emin to emax.
            
            Returns the integral in units of [erg^e_weight/cm^-3] """
        int = logsimps(lambda e: e**e_weight*self(e, units=False), self.emin, self.emax, sed_config.PER_DECADE)
        return int*(u.erg**(e_weight+1)*self.units() if units else 1)

class BlackBody(ThermalSpectrum):

    @staticmethod
    def compute_energy_density(kT):
        """ Comparing the formula for a blackbody spectrum
            with prefactor 

                pref = 8pi/(hc)^3

            to the fomrula for a general thermal spectrum:

                pref = 15*U/(pi*kT)^4,

            we find that for a blackbody spectrum,
            we have a thermal spectrum with

                U = (8*pi/(hc)^3)*(pi*kT)^4/15. """
        h=u.planck
        c=u.speed_of_light
        pi=np.pi
        return (8*pi/(h*c)**3)*((pi*kT)**4/15)


    def __init__(self,kT=None,T=None):
        """ Implement a blackbody spectrum.

            The formula for the blackbody spectrum is 
            
            n(E)=((8pi)/(hc)^3)*E^2/(exp(E/kT)-1)

            where 
              * n(E) is the number of photons per unit energy per unit volume,
              * kT is the temperature of the photons

            This formula is on the top of page 208 in R&L
        """
        if kT is not None: kT = kT
        elif T is not None: kT = u.boltzmann*T
        else: raise Exception("kT or k must be passed to BlackBody")

        energy_density=BlackBody.compute_energy_density(kT)
        super(BlackBody,self).__init__(energy_density=energy_density, kT=kT)


class CMB(BlackBody):
    """ The CMB is a blackbody spectrum with temperature 2.725K.

        Note, the energy density for a CMB spectrum is 0.26 eV/cm^3:
        
        >>> cmb = CMB()
        >>> print u.repr(cmb.integrate(units=True,e_weight=1),'eV/cm^3','%.2f')
        0.26 eV/cm^3
    """
    def __init__(self): super(CMB,self).__init__(T=2.725*u.kelvin)



if __name__ == "__main__":
    import doctest
    doctest.testmod()

