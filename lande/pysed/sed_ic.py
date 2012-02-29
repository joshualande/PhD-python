""" Module to compute inverse compton radiation
    for a given electron and photon spectrum.
    
    Author: Joshua Lande <joshualande@gmail.com>
"""
import numpy as np

from . sed_integrate import dbllogsimps
from . sed_spectrum import Spectrum
from . import sed_config
from . import units as u

class InverseCompton(Spectrum):
    """ The inverse compton radiation an electron spectrum
        and photon spectrum. """

    # default energy range = all energies
    vectorized = False

    def __init__(self, electron_spectrum, photon_spectrum):
        print 'The IC code needs to be validated and the formulas inspected + documented'

        self.electron_spectrum = electron_spectrum
        self.photon_spectrum = photon_spectrum

        self.mc2 = float(u.electron_mass*u.speed_of_light**2/u.erg)

        # this formula is basically 7.28a in R&L with the difference that

        # prefactor has units cm^3 s^-1
        self.pref = 2*np.pi*u.r0**2*u.speed_of_light
        self.pref = float(self.pref/(u.cm**3*u.second**-1))

    @staticmethod
    def F(q,gamma_e):
        """ This is equation 2.48 in Blumenthal & Gould. """
        return 2*q*np.log(q)+(1+2*q)*(1-q) + 0.5*(gamma_e*q)**2*(1-q)/(1+gamma_e*q)

    def _spectrum(self, scattered_photon_energy):
        """ Calculates the inverse compton spectrum expected
            from a sinle electron and an arbitrary photon spectrum. 
            
            Returns [ph/s/scattered photon energy]. """


        def integrand(electron_energy, target_photon_energy):
            """ electron_energy in erg
                target_photon_energy in erg
            
                return [ph/s/incident photon energy/scattered photon energy] 
                for a single electron
                
                in units [s^-1 erg^-3]
            """

            electron_gamma = electron_energy/self.mc2

            gamma_e = 4*target_photon_energy*electron_energy/(self.mc2)**2
            q=scattered_photon_energy/(electron_energy*gamma_e*(1-scattered_photon_energy/electron_energy))

            # Note about units:
            #  photon_spectrum has units 'photons/erg/cm^3'
            #  electron spectrum has units 'electrons/erg'
            #  F is unitless
            # so the integrand (pref*photon_spectrum*electron_spectrum*F)
            # has units (cm^3 s^-1 erg^-1) * (ph erg^-1 cm^-3) * (el erg^-1) * (1) = ph s^-1 erg^-3

            kinematically_allowed=(q<=1)&(q>=1./(4*electron_gamma**2))

            return np.where(kinematically_allowed,
                            self.pref*
                            electron_gamma**-2*target_photon_energy**-1*
                            self.photon_spectrum(target_photon_energy, units=False)*
                            self.electron_spectrum(electron_energy, units=False)*
                            self.F(q,gamma_e),
                            0)

        # Return the integrand integrated over photon and electron energy.
        # Note, integrand is in units of s^-1 erg^-3 so the twice
        # integration over energy gets the total number of emitted photons
        # per unit time per unit energy [s^-1 erg-^1]

        return dbllogsimps(
            f=integrand,
            xmin = self.electron_spectrum.emin,
            xmax = self.electron_spectrum.emax,
            ymin = self.photon_spectrum.emin,
            ymax = self.photon_spectrum.emax,
            per_decade = sed_config.PER_DECADE)

    @staticmethod
    def units_string(): return '1/s/erg'

