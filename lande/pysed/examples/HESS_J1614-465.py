""" Computes an SED of HESS J1640-465 according 
    to a time-dependent SED model described
    in Lemiere et al 2009 (http://arxiv.org/pdf/0910.2652v1)

    Author: Joshua Lande <joshualande@gmail.com>
"""
from pysed.sed_ic import InverseCompton
from pysed.sed_synch import Synchrotron
from pysed.helper import linspace_unit
from pysed.sed_spectrum import Constant
from pysed.sed_spectrum import CompositeSpectrum
from pysed.sed_thermal import CMB,ThermalSpectrum
import pysed.units as u


class InjectionRate(object):
    """ Model the injection spectrum of electrons
        following equation (1), (2), and (4) of Lemiere et al 2009
    """
    def __init__(self,
                 Edot, time,
                 slowing_down_time,
                 breaking_index,
                 powerlaw_index,
                 magnetization,
                 emin):
        r""" Parameters:
                Edot: assumed Pulsar energy injection rate at time time
                time: time at which pulsar injeciton rate is Edot
                slowing_down_time: parameter \Tau-0 in formula (1) of Lemiere et al 2009
                breaking_index: parameter \Tau-0 in formula (1) of Lemiere et al 2009
                emin: minimum electron energy
        """
        self.slowing_down_time = slowing_down_time
        self.breaking_index = breaking_index
        self.powerlaw_index = powerlaw_index
        self.magnetization = magnetization
        self.emin = emin

        self.E0dot = 1
        self.E0dot = Edot/self.total_power(time)

    def total_power(self,time):
        """ Computes the total power injected into the pulsar at
            a given time according to equatino (1) of Lemiere et al 2009.
        """
        E0dot = self.E0dot
        tau0 = self.slowing_down_time
        b = self.breaking_index
        t = time
        return E0dot/(1 + float(t/tau0))**((b+1)/(b-1))

    def get_spectrum(self,time):
        s = self.magnetization
        total_power = Edot = self.total_power(time)

        # from section 3.3 of text
        electron_power = 1./(s+1)*total_power

        # equation (4) in text
        e = u.electron_charge
        c = u.speed_of_light
        emax = np.sqrt(s/(1.+s))*e*np.sqrt(Edot/c)

        return PowerLaw(
            power = electron_power,
            units_string = '1/erg/s',
            emin = self.emin,
            emax = self.emax)

class MagneticField(object):
    pass

class PulsarMagneticField(MagneticField):
    def __init__(self,
                 magnetic_field,
                 time,
                 slowing_down_time,
                 alpha):
        """ Model the evolution of the Pulsar
            magnetic field as a function of time according
            to equation (5) of Lemiere et al 2009
        """

        self.slowing_down_time = slowing_down_time
        self.alpha = alpha

        # set the magnetic field normalization 
        self.initial_magnetic_field = 1
        self.initial_magnetic_field = magnetic_field/self(time)

    def __call__(self,time):
        """ Implement equation (5) of Lemiere et al 2009. """
        B0 = self.initial_magnetic_field
        tau0 = self.slowing_down_time
        alpha = self.alpha
        return B0/(1+float(time/tau0)**alpha)


class PWNEvolution(object):
    """ Class to model the evolution of a population of electrons
        injected into a PWN by a Pulsar.
    """

    def __init__(self, 
                 injection_rate,
                 photon_fields,
                 magnetic_field,
                 pulsar_age,
                 ntsteps,
                 emin, emax):
        """ Parameters:
                injection_rate: object of type InjectionRate
                photon_fields: object of type ThermalSpectrum 
                    (or CompositeSpectrum of ThermalSpectrum objects)
                magnetic_field: object of type MagneticField
                pulsar_age: age of the pulsar
                ntsteps: number of time steps for simulation.
                emin, emax: energy range over which to perform simulaiton.
        """
        zero = Constant(value=0,units_string='1/erg')
        self.electrons = NumericalSpectrum(
            spectrum=zero, 
            emin=emin, 
            emax=emax,
            per_decade = sed_config.PER_DECADE
            )

        time_edges = linspace_unit(0, pulsar_age, ntsteps+1)
        self.times = 0.5*(time_edges[1:] + time_edges[:-1])
        self.dtimes = time_edges[1:] - time_edges[:-1]

        self.evolve()

    def evolve(self):

        electrons = self.electrons

        for time, dtime in zip(self.times,self.dtimes):

            # (a) Calculate additional energy input from pulsar

            pulsar_power = injection_rate.get_spectrum(time)

            electrons += pulsar_power*dtime

            # (b) Calculate new magnetic field

            synch = Synchrotron(electron_spectrum=electrons, 
                                magnetic_field=magnetic_field(time))

            ic = InverseCompton(electron_spectrum=electrons,
                                photon_spectrum=photon_fields)

            # (c) Calculate energy loss due to Synch + IC,

            loss_rate = synch.energy_loss_rate() + ic.energy_loss_rate()

            # (d) Update electron distribution due to energy losses

            electrons -= loss_rate*dtime


if __name__ == '__main__':

    from pysed import sed_config
    #sed_config.PER_DECADE = 100

    # (a) Create initial distrubitons for particles

    pulsar_age = 15*u.kiloyear # table 3 of text
    slowing_down_time = 500*u.year # section 3.3 of text

    injection_rate = InjectionRate(
        Edot = 4e36*u.erg/u.second, # table 3  of text
        time = pulsar_age, 
        slowing_down_time = slowing_down_time, 
        breaking_index=3, # section 3.3 of text
        powerlaw_index = 2.4, # section 3.3 of text
        magnetization = 0.45, # table 3 of text
        emin=50*u.GeV, # table 3 of text
        )

    magnetic_field = PulsarMagneticField(
                 magnetic_field = 6*u.microgauss,
                 time = pulsar_age,
                 slowing_down_time = slowing_down_time,
                 alpha = 0.45 # table 3 of text
        )

    cmb = CMB()
    infrared = ThermalSpectrum(kT=XXX, energy_density=0.3*u.eV*u.cm**-3)
    optical = ThermalSpectrum(kT=XXX, energy_density=0.9*u.eV*u.cm**-3)
    photon_fields = CompositeSpectrum(cmb, infrared, optical)

    # (b) evolve in time

    print 'Not sure if there are enough time steps.'
    ntsteps = 15

    print """Not sure correct way to compute fo particles . Todo: ask people. Here is an interesting reference for particle cooling for IC in Klein Nishena limit
    http://iopscience.iop.org/0004-637X/509/1/212/pdf/37990.pdf. """

    evolution = PWNEvolution(
                 injection_rate = injection_rate,
                 photon_fields = photon_fields,
                 magnetic_field = magnetic_field,
                 pulsar_age = pulsar_age,
                 time_step = ntsteps,
                 emin=1e-9*u.eV, # emin + emax are the range of Figure 9 in text
                 emax=1e16*u.eV)

    # (c) plot current time

