""" Code that calculates the synchrotron radiation
    from a given electron spectrum and magnetic
    field.

    Author: Joshua Lande <joshualande@gmail.com>
"""
from numpy import pi, sqrt,inf,sin
from scipy import integrate,special

from . sed_spectrum import Spectrum
from . sed_integrate import halfdbllogsimps
from . sed_cache import FunctionCache
from . import sed_config
from . import sed_units as u

class Synchrotron(Spectrum):
    """ Calculates the syncrotron power radiated
        by a spectrum of electrons in a magnetic field. 

        This treatment follows section 6.2 in R&L and especially  
        the formulas 6.17c, 6.18, and 6.31c
    
        The spectrum is in unit of energy per unit time
        per unit frequency emitted by a single electron.

        Something about averaging over pitch angle
        and 

        """


    # default energy range = all energies
    vectorized = False

    def _F(x):
        """ This is F(x) defined in equation 6.31c in R&L.
            
            F(x) = x*int(K_5/3(x)dx) where the integral goes from x to infinity.
        
            for some reason, special.kv(5/3,1e10) is NaN, not 0 ???
            for now, just clip the function above 1e5 to be 0. 
            
            This function can be evaluated in mathematica using the following command

                F[x_] := N[x*Integrate[BesselK[5/3, y], {y, x, Infinity}]]

            From mathematica, we find that

                      x         F(x)
                  ----- ------------
                    0.1     0.818186
                      1     0.651423
                     10  0.000192238
                    100            0

            Comparing our function to the Mathematica integral, we find

                >>> np.allclose(Synchrotron.F([0.1,1,10,100]), [0.818186, 0.651423, 0.000192238,0], rtol=1e-4, atol=1e-4)
                True

            Note, this function is _F so that the docstring will get executed.
        """
        if x>1e5: return 0
        return x*integrate.quad(lambda j: special.kv(5./3,j),x,inf)[0]
    F=FunctionCache(_F, xmin=0, xmax=20, npts=1000, fill_value=0)

    def __init__(self, electron_spectrum, magnetic_field):
        print 'The IC code needs to be validated and the formulas inspected + documented.'

        self.electron_spectrum = electron_spectrum

        e = u.electron_charge
        B = magnetic_field
        m = u.electron_mass
        c = u.speed_of_light

        # prefactor from Sturner et al 1997 formula 22:
        # Returns power/energy in units of erg/s/erg
        self.pref = sqrt(3)*e**3*B/(u.planck*m*c**2)

        self.pref = float(self.pref/(u.erg*u.second**-1*u.erg**-1))

        # This is formula 6.17 in R&L except for the gamma**2
        omega_c = 3*e*B/(2*m*c)
        self.energy_c_pref = float(u.hbar*omega_c/u.erg)

        self.mc2 = m*c**2
        self.mc2_in_erg = float(self.mc2/u.erg)

    def _spectrum(self, photon_energy):
        """ return total power per emitted per unit energy by
            the spectrum of electrons.
            
            Integrate the power by a single electron over the spectrum of electrons. """

        def integrand(electron_energy, theta):
            """ The syncrotron radiation integrand
                as a function of 

                * electron energy, measured in erg
                * the pitch angle theta, measured in radians.
            """
            electron_gamma = electron_energy/self.mc2_in_erg
            sin_theta = sin(theta)

            energy_c = self.energy_c_pref*electron_gamma**2*sin_theta

            # power_per_energy in units of erg/s/erg
            power_per_energy = self.pref*sin_theta**2*Synchrotron.F(photon_energy/energy_c)
            # divide by photon energy to get photons/energy for a single
            # electron (in units of ph/erg/s)
            photons_per_energy = power_per_energy/photon_energy

            # return [ph/s/photon energy/electron energy] = 
            #             [ph/s/electron/photon energy]*[number of electrons/electron energy]
            return photons_per_energy*self.electron_spectrum(electron_energy, units=False)

        # integrate [ph/s/photon energy/electron energy] over electron energies
        # Returns is photon flux per energy per time [1/erg/s]
        emin = self.electron_spectrum.emin
        emax = self.electron_spectrum.emax

        # integrate the function over the electron distribution and over
        # pitch angle.
        # Note: integrate energy in log space, angle in linear space.
        return halfdbllogsimps(integrand, 
                               xmin=emin, xmax=emax, 
                               ymin=0, ymax=pi/2,
                               x_per_decade=sed_config.PER_DECADE,
                               y_npts=10)

    def energy_loss(self, energy):
        """ Returns the energy loss due to synchrotron radiation
            in units of erg s^-1. """
        raise Exception("Not implemented yet...")

    @staticmethod
    def units_string(): return '1/erg/s'

if __name__ == "__main__":
    import doctest
    doctest.testmod()

