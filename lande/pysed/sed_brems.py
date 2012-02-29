""" Calculates the gamma ray emisison predicted
    by Bremsstrahlung radiation using the prescription from
    Barington et al 1999 (http://arxiv.org/abs/astro-ph/9810158)

    Author: Joshua Lande <joshualande@gmail.com>
"""
import numpy as np
from numpy import sqrt, log

from . sed_cross_section import CrossSection
from . sed_spectrum import Spectrum
from . sed_integrate import logsimps
from . sed_relativity import gamma_to_beta
from . import sed_config
from . import units as u

# Energy of an electron at rest in erg
electron_rest_energy_erg = float(u.electron_mass*u.speed_of_light**2/u.erg)

class BremsEPCrossSection(CrossSection):

    r""" Computes the photon spectrum due to electron ion Bremsstrahlung. 

        The electron-ion cross ection is caused by
        an electron scattering off of a roughly stationary atom.

        The cross section for an electron of energy $\gamma m_e c^2$
        (and velocity \beta c) to scatter off an ion of nuclear charge
        $Ze$ and produce a photon of energy E=\hbar\omega
        is given in the non-relativistic Born approximation
        by the Bethe-Heitler cross section

            \frac{d\sigma}{d\omega} = 
                \tfrac{16}{3} Z^2 
                \frac{\alpha r_0^2}{\omega}
                \frac{1}{\beta^2}
                \log \frac{\beta + \beta'}{\beta - \beta'}

        This is equation 15-99 in Jauch and Rohrlich (1980)
        but was apparently originally derived in in (Heitler 1954)

        In this formula, 
        \alpha is the fine strucutre constant (\alpha=e^2/(\hbar c)),
        r0 is the classical electron radius (r_0=e^2/(m_e c^2), and 
        \beta'c is the final electron velocity. Since this formula
        is only true in the non-relativistic limit:

            0.5 m_e v'^2 = 0.5 m_e v^2 - E_\gamma

        and from this we find that \beta' is
        
            \beta' = \sqrt{\beta^2 - 2 E_\gamma/m_e c^2}

        So the total cross section, differential in photon energy, is 

            \frac{d\sigma}{dE_\gamma} = 
                \tfrac{16}{3} Z^2 \alpha r_0^2
                \frac{1}{E_\gamma} \frac{1}{\beta^2}
                \log \frac{\beta + \beta'}{\beta - \beta'}

        Although this formula is only true in the non-relativistic limit,
        Baring et al 1998 section 3.2 justifies using this formula for
        electron-ion scattering for all electron energies.

        Note, from kinematic constrains, we must have 

            1/2 m_e beta^2 c^2 > E_\gamma

        which corresponds to the constraing:

                beta^2 > 2*E_\gamma/m_e c^2

        Implementation note, this object calculates exclusively the 
        electron-proton cross section (Z=1)
    """

    prefactor = (16./3)*u.alpha*u.r0**2
    prefactor = float(prefactor/u.cm**2)

    def __call__(self, electron_energy, photon_energy):
        """ Assumes electron_energy and photon_energy both in erg.
        
            Return the electron proton Bremsstrahlung cross section 
            in units of cm^-2 erg^-1. """

        gamma = electron_energy / electron_rest_energy_erg
        beta = gamma_to_beta(gamma)
        beta_prime = sqrt(beta**2 - 2*photon_energy/electron_rest_energy_erg)

        sigma = self.prefactor*\
                (1/photon_energy)*(1/beta**2)*\
                log((beta+beta_prime)/(beta-beta_prime))

        allowed = (beta**2 >=2*photon_energy/electron_rest_energy_erg)

        return np.where(allowed, sigma, 0)


class BremsEECrossSectionNR(CrossSection):
    """ Non-relativisitc electron-electron Bremsstrahlung
        radiation cross section. 
        
        The non-relativistic electron-electron Bremsstrahlung
        cross section is an approximation defined in A5 of the appendix
        of Baring et al 1998.
        """

    prefactor = (4./15)*u.r0**2*u.alpha
    prefactor = float(prefactor/u.cm**2)

    @staticmethod
    def B(gamma):
        return 1 + 0.5*(gamma**2-1)

    @staticmethod
    def C(gamma,x):
        beta = gamma_to_beta(gamma)
        g,b=gamma,beta
        return 10*x*g*b*(2+g*b)/(1+x**2*(g**2-1))

    @staticmethod
    def F(gamma,x):
        B = BremsEECrossSectionNR.B(gamma)
        C = BremsEECrossSectionNR.C(gamma,x)
        return B*(17 - 3.*x**2/(2-x)**2 - C)*sqrt(1-x) + \
                (12*(2-x) - 7*x**2/(2-x) - 3*x**4/(2-x)**3)*log((1+sqrt(1-x))/sqrt(x))

    def __call__(self, electron_energy,photon_energy):
        """ Returns a non-relativistic approximation to the
            electron-electron Bremsstrahlung (formula A5 of Baring et al 1998). 

            both electron_energy and photon_energy assuemd to be in units of erg.
            Return is in units of cm**2 erg**-1 """

        epsilon_gamma = photon_energy/electron_rest_energy_erg
        electron_gamma = electron_energy/electron_rest_energy_erg

        x = 4*epsilon_gamma/(electron_gamma**2-1)
        F = BremsEECrossSectionNR.F(electron_gamma, x)

        # Note, the formula in Barington et al 1998 has a factor 
        # of (1/epsilon_gamma) which makes the units of the formula incorrect
        # cm**2 instead of cm**2 erg**-1. I think instead, it is better to 
        # have a factor of (1/electron_energy) in the formula

        sigma = self.prefactor*(1/electron_energy)*F

        # Apply kinematic constraint from equation A5 of Barington et al 1998
        allowed = (0 < epsilon_gamma) & (epsilon_gamma < 0.25*(electron_gamma**2-1))
        
        return np.where(allowed, sigma, 0)

class BremsEECrossSectionRel(CrossSection):
    """ Relativisitc electron-electron Bremsstrahlung
        radiation cross section. """

    prefactor = u.r0**2*u.alpha
    prefactor = float(prefactor/u.cm**2)

    @staticmethod
    def sigma_1(epsilon, gamma):
        r""" Formula A2 in Barington et al 1999
        
            With the exception that I am pretty sure the
            first \frac{1}{epsilon_\gamma} should instead be
            \frac{1}{E_\gamma}. """
        e,g = epsilon, gamma
        photon_energy = E = epsilon * electron_rest_energy_erg
        r2alpha = BremsEECrossSectionRel.prefactor
        return (4*r2alpha/E)*(1 + (1./3 - e/g)*(1-e/g))*(log(2*g*(g-e)/e)-0.5)

    @staticmethod
    def sigma_2(epsilon, gamma):
        r""" Formula A3 in Barington et al 1999.
            
            With the exception that I am pretty sure the
            first \frac{1}{epsilon_\gamma} should instead be
            \frac{1}{E_\gamma}. """
        e,g = epsilon, gamma
        photon_energy = E = epsilon * electron_rest_energy_erg
        r2alpha = BremsEECrossSectionRel.prefactor

        low = 16*(1 - e + e**2)*log(g/e) - 1/e**2 + 3/e - 4 + 4*e - 8*e**2 \
                -2*(1-2*e)*log(1-2*e)*(1/(4*e**3) - 1/(2*e**2) + 3/e - 2 + 4*e)

        high = (2/e)*((4-1/e+1/(4*e**2))*log(2*g) - 2 + 2/e - 5/(8*e**2))

        return (r2alpha/(3*E))*np.where(e<=0.5, low, high)

    @staticmethod
    def A(epsilon, gamma):
        """ Formula A4 in Barington et al 1999. """
        e,g = epsilon, gamma
        return 1 - (8./3)*((g-1)**0.2/(g+1))*(e/g)**(1./3)

    def __call__(self, electron_energy,photon_energy):
        """ Formula A1 in Barington et al 1999. """

        epsilon_gamma = photon_energy/electron_rest_energy_erg
        electron_gamma = electron_energy/electron_rest_energy_erg

        sigma_1 = BremsEECrossSectionRel.sigma_1(epsilon_gamma, electron_gamma)
        sigma_2 = BremsEECrossSectionRel.sigma_2(epsilon_gamma, electron_gamma)

        A = BremsEECrossSectionRel.A(epsilon_gamma, electron_gamma)

        sigma = (sigma_1+sigma_2)*A

        # apply kinematic constraint, electron kinetic
        # energy must be greater than photon energy.
        allowed = electron_gamma -1 > epsilon_gamma

        return np.where(allowed, sigma, 0)

class BremsEECrossSection(CrossSection):
    """ Computes the electron electron Bremsstrahlung cross
        section for all electron eneriges by interplolating
        between the relativistic and non-relativistic 
        regimes. This is described in the Appendix of Barington et al 1999.
    """
    def __init__(self, *args, **kwargs):
        self.cross_section_nr  = BremsEECrossSectionNR(*args, **kwargs)
        self.cross_section_rel = BremsEECrossSectionRel(*args, **kwargs)

        self.two_mev_erg = float(2*u.MeV/u.erg)

    def __call__(self, electron_energy,photon_energy):
        """ Interpolate between the non-relativistic and
            relativistic regimes. """
        sigma_rel = self.cross_section_rel(electron_energy,photon_energy)
        sigma_nr = self.cross_section_nr(electron_energy,photon_energy)

        return np.where(electron_energy > self.two_mev_erg, sigma_rel, sigma_nr)


class Bremsstrahlung(Spectrum):
    r""" Computes the  Bremsstrahlung radiation using the prescription from
        Barington et al 1999 (http://arxiv.org/abs/astro-ph/9810158)

        This computes the e-e and the e-p Bremsstrahlung radiation for
        a spectrum of electrons hitting an (assumed stationary) density
        of hydrogen and helium particles.

        For an input electron spectrum (dN_e/dE_e) and a given
        number density $n_t$ of stationary target particles, the distribution
        of scattered photons is computed by integrating the yield over
        the electron distribution:

        \frac{dN}{dE_\gamma t} = \int_0^\infty dE_e \beta_e c \times
                                 ((n_P + 4n_{He}) \frac{d\sigma_{e-p}}{d E_e} +
                                  (n_P + 2n_{He}) \frac{d\sigma_{e-e}}{d E_e}) \times
                                 \frac{dN_e}{dE_e}
    """

    vectorized = False

    def __init__(self, electron_spectrum, hydrogen_density, helium_density):
        """ electron_spectrum:  a Spectrum object
                which returns the number of electrons per unit energy 
            hydrogen_density: density of the target hydrogen atoms [number/Volume]
            helium_density: density of the target helium atoms [number/Volume]
        """
        self.electron_spectrum = electron_spectrum

        self.hydrogen_density = float(hydrogen_density/u.cm**-3)
        self.helium_density = float(helium_density/u.cm**-3)

        print r'what to do about divergence as \omega -> 0'

        self.e_p_cross_section = BremsEPCrossSection()
        self.e_e_cross_section = BremsEECrossSection()

        self.speed_of_light_cgs = float(u.speed_of_light/(u.cm*u.seconds**-1))
        

    def _spectrum(self, photon_energy):
        """ Returns Bremsstrahlung due to a distribution of electrons
            in units of erg^-1 s^-1 """

        def integrand(electron_energy):
            gamma = electron_energy/electron_rest_energy_erg
            beta = gamma_to_beta(gamma)

            c = self.speed_of_light_cgs 
            nP = self.hydrogen_density 
            nHe = self.helium_density 

            sigmaEE = self.e_e_cross_section(electron_energy, photon_energy)
            sigmaEP = self.e_p_cross_section(electron_energy, photon_energy)
            dnde = self.electron_spectrum(electron_energy,units=False)

            return beta*c*((nP + 4*nHe)*sigmaEP + (nP + 2*nHe)*sigmaEE)*dnde

        emin = self.electron_spectrum.emin
        emax = self.electron_spectrum.emax

        return logsimps(integrand, emin, emax, sed_config.PER_DECADE)

    @staticmethod 
    def units_string(): return '1/s/erg'
