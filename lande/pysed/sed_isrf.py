""" This file provides a class for parsing the GALPROP 
    Interstellar Radiation Field (ISRF) for computing
    the photon fields important for Inverse Compton
    of photons in astrophysical objects.

    Author: Joshua Lande <joshualande@gmail.com>
"""
from math import pi, exp
import pyfits
import numpy as np
from scipy.special import lambertw

from . import units as u
from . sed_thermal import ThermalSpectrum
from sed_helper import argmax_unit

class ISRF(object):

    def __init__(self, isrf):
        """ isrf is the GALPROP ISRF fits mapcube.

            The file I am using is 
                
                MilkyWay_DR0.5_DZ0.1_DPHI10_RMAX20_ZMAX5_galprop_format.fits.gz

            and can be found at 

                http://galprop.stanford.edu/resources.php?option=data
        """

        self.isrf = pyfits.open(isrf)[0]


        for number, line in [
            [18,'Units micron eV cm^-3 micron^-1'],
            [19,"R:z:log10(wavelength):component"],
            [20,"kpc:kpc:log10(micron):integer"],
            [21,"Component 1 = optical"],
            [22,"Component 2 = infrared"],
            [23,"Component 3 = CMB"]]:

            if self.isrf.header[number] != line: raise Exception("Unrecognized header for ISRF.")

    def R_to_index(self, R):
        R_internal = float(R/u.kpc)

        h = self.isrf.header
        start, delt, num = h['CRVAL1'], h['CDELT1'], h['NAXIS1']

        index = (R_internal - start)/delt
        return index
        

    def z_to_index(self, z):
        z_internal = float(z/u.kpc)

        h = self.isrf.header
        start, delt, num = h['CRVAL2'], h['CDELT2'], h['NAXIS2']
        
        index = (z_internal - start)/delt
        return index

    def get_wavelength(self):
        """ Get the wavelengths in the mapcube. """

        h = self.isrf.header
        start, delt, num = h['CRVAL3'], h['CDELT3'], h['NAXIS3']

        wavelength = 10**(start + delt*np.arange(num))
        return u.tosympy(wavelength, u.micron)

    def get_energy(self):
        wavelength = self.get_wavelength()
        
        f=wavelength.applyfunc(lambda w: u.planck*u.speed_of_light/w)

        return wavelength.applyfunc(lambda w: u.planck*u.speed_of_light/w)


    def get(self, component, R, z):
        """ Get the ISRF for a given compoenent and a given galactic position. 
        
            Returns photon density per unit energy (photons/volume/energy). """

        R_index = self.R_to_index(R)
        z_index = self.z_to_index(z)

        print 'for now, clip. This is BAD. Should be an interpolation.'
        R_index = R_index
        z_index = z_index

        radiation = self.isrf.data[component,:,z_index,R_index]

        radiation = u.tosympy(radiation,u.eV*u.cm**-3)

        energy = self.get_energy()
        inverse_energy = energy.applyfunc(lambda x: x**-1)

        # convert from energy output per unit energy (energy/cm^3)
        # to photons per unit energy (ph/cm^3/energy)
        # by dividing by two factors of energy.
        print 'This dividing by 2 factors of energy needs to be validated!'
        radiation = radiation.multiply_elementwise(inverse_energy)
        radiation = radiation.multiply_elementwise(inverse_energy)

        return radiation

        
    def get_optical(self, *args, **kwargs): return self.get(0, *args, **kwargs)
    def get_infrared(self, *args, **kwargs): return self.get(1, *args, **kwargs)
    def get_CMB(self, *args, **kwargs): return self.get(2, *args, **kwargs)

    def estimate(self, *args, **kwargs):
        """ Estimate a thermal sectrum
            based upon the peak energy in
            the GALPROP spectrum. """

        energy = self.get_energy()
        intensity = self.get(*args, **kwargs)

        def max_energy_to_temperature(energy):
            """ Converts the peak energy to the peak. 

                This formula is derived noting that the funcional dependence
                is in the term x^2/(Exp[x]-1) where x=E/kT.

                From mathematatica, we find:

                In[0]:= ArgMax[x^2/(Exp[x] - 1), x]
                Out[13]= 2 + ProductLog[-(2/E^2)]

                Which is approximately 1.59362

                So the peak energy is at x ~ 1.6, or kT~E/1.6

                Implementation Note: ProductLog in Mathematica is 
                the Lambert W function (scipy.special.lambertw) """
            print 'fix documentation'
            return energy/float(2 + lambertw(-2/exp(2)))


        # find peak in spectrum
        argmax = argmax_unit(intensity)

        max_energy = energy[argmax]
        max_intensity = intensity[argmax]

        kT = max_energy_to_temperature(max_energy)

        energy_density = ((pi*kT)**4/(15*max_energy**2))*(exp(float(max_energy/kT))-1)*max_intensity

        return ThermalSpectrum(
            energy_density=energy_density,
            kT=kT)

    def estimate_optical(self, *args, **kwargs): return self.estimate(0, *args, **kwargs)
    def estimate_infrared(self, *args, **kwargs): return self.estimate(1, *args, **kwargs)
    def estimate_CMB(self, *args, **kwargs): return self.estimate(2, *args, **kwargs)
