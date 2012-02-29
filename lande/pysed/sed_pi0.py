""" Calculates the gamma ray emisison predicted
    by pi0 decay.
    for a given electron and photon spectrum
    using the method described in 
    Kamae et al 2006 (http://arxiv.org/abs/astro-ph/0605581)
    .
    

    Author: Joshua Lande <joshualande@gmail.com>
"""
import numpy as np

from . sed_spectrum import Spectrum
from . sed_integrate import logsimps
from . sed_cross_section import CrossSection
from . sed_relativity import gamma_to_beta
from . import sed_config
from . import sed_units as u

class PPCrossSection(CrossSection):
    """ Object to calculate the 
        proton-proton cross section
        for decaying into gammas.

        This module using a parameterization
        of numerical pi0 decay codes described
        in Kamae et al 2006:

            http://arxiv.org/abs/astro-ph/0605581

        And performs the calculation by wraping
        the numerical codes they provide in 
        the swig interface to cparamlib

        To install cparamlib so that it can be accessed
        through python, please visit this very nice page:

            http://homepages.spa.umn.edu/~nkarlsson/cparamlib/
    """

    def __init__(self):
        # module to calculate Pi0 cross section
        from cparamlib.cparamlib import ID_GAMMA
        from cparamlib.ParamModel import ParamModel

        self.param = ParamModel(Tp=0,particle=ID_GAMMA)

        self.erg_to_gev = float(u.erg/u.GeV)
        self.millibarn_to_cm2 = float(u.millibarn/u.cm**2)

    def __call__(self, proton_energy,photon_energy):
        """ Computes the proton proton cross section to decay into a gamma.

            proton_energy is the energy of the incident proton, in units of erg
            photon_energy is the energy of the resultant gamma, in units of erg

            The return cross section is d(sigma)/dE where E is the photon energy,
            sigma is in units of cm**2, and E is in units of erg.

            Implementation Note:

                Sigma_incl_tot returns the photon spectrum including all
                processes.

                sigma_incl_tot returns dsigma/dlog(E) in units of mb.
                The two inputs must be in units of GeV and dlog(E) is
                calculated (I assume) in units of GeV
        """
        

        photon_energy_gev = photon_energy*self.erg_to_gev
        proton_energy_gev = proton_energy*self.erg_to_gev

        # Currently, cparamlib is not vecotrized, so vectorize it here :(
        if isinstance(proton_energy_gev,np.ndarray):
            dsigmadloge = np.asarray([self.param.sigma_incl_tot(photon_energy_gev, i) \
                                      for i in proton_energy_gev])
        else:
            dsigmadloge = self.param.sigma_incl_tot(photon_energy_gev, proton_energy_gev)

        # convert from cross section per log(energy) in units of millibarn
        # to cross section per(energy) in units of cm^2

        dsigmade = self.millibarn_to_cm2*dsigmadloge*(1/photon_energy)
        return dsigmade



class Pi0Decay(Spectrum):
    """ Computes the Pi0 decay flux
        for a spectrum of protons
        hitting a density of particles. """

    # default energy range = all energies
    vectorized = False

    def __init__(self,
                 proton_spectrum, 
                 hydrogen_density,
                 scaling_factor):
        """ 
        
            proton_spectrum:  differental number of protons.
            hydrogen_density: density that input spectrum is hitting
        
            scaling_factor: 
                the unitless factor is introduced to account for healium and heavier nulclei. 
            
                In Abdo et al 2009, the SED of W51C was computed setting
                the factor to 1.85 (http://arxiv.org/abs/0910.0908)

                Mori et al 2009 says that the factor should be set between 1.8 and 2:
                (http://arxiv.org/abs/0903.3260)

        """
        print 'The Pi0-decay code needs to be validated and the formulas inspected + documented'


        self.cross_section = PPCrossSection()

        self.proton_spectrum = proton_spectrum

        self.scaling_factor = scaling_factor

        self.hydrogen_density = hydrogen_density

        self.prefactor = self.scaling_factor*self.hydrogen_density*u.speed_of_light
        self.prefactor = float(self.prefactor/(u.cm**-2*u.seconds**-1))

        self.proton_rest_energy_erg = float(u.proton_mass*u.speed_of_light**2/u.erg)

    def _spectrum(self,photon_energy):
        """ Return spectrum in units of s^-1 erg^-1. """

        def integrand(proton_energy):

            proton_gamma = proton_energy/self.proton_rest_energy_erg

            # in case there are unphysically low proton eneriges (gamma<1), set beta=0
            # This is a bit ugly, but these protons have no cross section so there
            # is no harm in including them in the computation.
            proton_beta = np.where(proton_gamma>=1,gamma_to_beta(proton_gamma),0)

            dnde=self.proton_spectrum(proton_energy, units=False) # protons erg^-1
            dsigmade=self.cross_section(proton_energy, photon_energy) # cm^2 erg^-1

            # Units: (s^-1 erg^-2) = (cm^-2 s^-1) x (protons erg^-1) x (cm^2 erg^-1)
            return self.prefactor*proton_beta*dnde*dsigmade

        emin = self.proton_spectrum.emin
        emax = self.proton_spectrum.emax
        return logsimps(integrand, emin, emax, sed_config.PER_DECADE)

    @staticmethod
    def units_string(): return '1/s/erg'
