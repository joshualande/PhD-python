""" Compute the SED of W51C using parameters from
    http://arxiv.org/abs/0910.0908

    SED for for hypothesis (c) InverseCompton

    Author: Joshua Lande <joshualande@gmail.com>
"""
from matplotlib import rc
rc('text', usetex=True)
rc('font', family='serif', serif="Computer Modern Romain")


import pylab as P
import numpy as np
from mpl_toolkits.axes_grid1 import AxesGrid
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText

from pysed.sed_particle import SmoothBrokenPowerLaw
from pysed.sed_ic import InverseCompton
from pysed.sed_pi0 import Pi0Decay
from pysed.sed_brems import Bremsstrahlung
from pysed.sed_synch import Synchrotron
from pysed.sed_spectrum import CompositeSpectrum
from pysed.sed_plotting import SEDPlotter
from pysed.sed_thermal import CMB,ThermalSpectrum
import pysed.sed_units as u

np.seterr(all='ignore')


def plot_photon_fields(type, **kwargs):
    """ Plot all the photon fields. """

    fig = P.figure(None,figsize=(5.5,4.5))
    axes = fig.add_subplot(111)

    plot_kwargs = dict(x_units_string='eV',y_units_string='eV/cm^3', e_weight=2)

    for k,v in kwargs.items():
        v.loglog(label=k, axes=axes, **plot_kwargs)
    axes.legend(loc=3)
    fig.savefig('w51C_%s_photon_fields.pdf' % type)

def print_photon_field(**kwargs):
    """ Sanity check, print otu temperature and energy density of the photon fields. """
    plot_field = lambda i: 'kT=%s, E=%s' % (u.repr(i.kT*u.erg,'eV'), u.repr(i.integrate(e_weight=1),'eV*cm^-3'))
    print 'Photon Fields:'
    for k,v in kwargs.items():
        print '\t%s: %s' % (k,plot_field(v))

def plot_electrons(type, electrons):
    electrons.loglog(e_weight=2,
                     x_units_string='MeV', y_units_string='MeV',
                     filename='w51C_%s_electrons.pdf' % type)

def sed(axes, label, type, delta_s, e_break, magnetic_field, hydrogen_density, Wp, We):
    """ Compute the SED for W51C in the Pi0-decay dominated hypothesis
        (a) in table 1 Fig.4 of the text. """

    electrons = SmoothBrokenPowerLaw(
        total_energy=We,
        index1 = 1.5, # p10 in text
        index2 = 1.5 + delta_s, # index2 = index1 + delta_s
        e_break = e_break,
        e_scale = 1*u.GeV,
        beta = 2.0,
        emin = 10*u.MeV, # from footnote to table 1
        emax = u.TeV, # Yasunobu told me this was the emax when I asked
        )

    protons = SmoothBrokenPowerLaw(
        total_energy=Wp,
        index1 = 1.5, # p10 in text
        index2 = 1.5 + delta_s, # index2 = index1 + delta_s
        e_break = e_break,
        e_scale = 1*u.GeV,
        beta = 2.0,
        emin = 10*u.MeV, # from footnote to table 1
        emax = u.TeV, # Yasunobu told me this when I asked
        )

    magnetic_field = magnetic_field 

    distance = 6*u.kiloparsec # from page 8 in text

    hydrogen_density = hydrogen_density
    helium_density = 0.1*hydrogen_density # yasunobu told me this when I asked

    # Photon fields take from the footnote of table 1 in the text
    cmb = CMB()
    infrared = ThermalSpectrum(kT=3e-3*u.eV, energy_density=0.9*u.eV*u.cm**-3)
    optical = ThermalSpectrum(kT=0.25*u.eV, energy_density=0.84*u.eV*u.cm**-3)
    photon_fields = CompositeSpectrum(cmb, infrared, optical)

    # Make some nice diagnostic plots
    plot_photon_fields(type, CMB=cmb, infrared=infrared, optical=optical)
    print_photon_field(CMB=cmb, infrared=infrared, optical=optical)
    plot_electrons(type, electrons)

    # Create the radiation processes

    synch = Synchrotron(electron_spectrum=electrons, 
                        magnetic_field=magnetic_field)

    ic = InverseCompton(electron_spectrum=electrons,
                        photon_spectrum=photon_fields)

    pi0 = Pi0Decay(proton_spectrum=protons, 
                   hydrogen_density = hydrogen_density,
                   scaling_factor = 1.85 # from p10 in the text
                  )

    brems = Bremsstrahlung(
        electron_spectrum = electrons,
        hydrogen_density = hydrogen_density,
        helium_density = helium_density)

    sed = SEDPlotter(
        distance=distance,
        x_units_string='eV',
        y_units_string='erg*cm^-2*s^-1',
        axes=axes,
        )

    # Overlay the Synchrotron and Inverse Compton radiation
    sed.plot(synch, color='black')
    sed.plot(pi0, color='black', dashes=[9,2])
    sed.plot(ic, color='black', dashes=[2,2])
    sed.plot(brems, color='black', dashes=[4,2])

    at = AnchoredText(label, frameon=False, loc=2, prop=dict(size=14))
    axes.add_artist(at)


if __name__ == '__main__':

    from pysed import sed_config
    #sed_config.PER_DECADE = 100

    fig = P.figure(None,(7,6))
    grid=AxesGrid(fig, 111,
                  nrows_ncols = (3, 1),
                  axes_pad = 0.1,
                  aspect=False,
                  share_all=True)

    # Setup up the axes to be the same as figure 4 in the publiation
    for axes in grid:

        axes.set_xlabel(r'E [eV]')
        axes.set_ylabel(r'$\nu f_\nu$ [erg cm$^{-2}$ s$^{-1}$]')

        axes.set_xscale('log')
        axes.set_yscale('log')

        axes.set_xlim(xmin=6e-7, xmax=2e12)
        axes.set_ylim(ymin=2e-13, ymax=2e-10)

        axes.xaxis.set_ticks([1e-6,1e-3, 1e-0, 1e3, 1e6, 1e9, 1e12])
        axes.yaxis.set_ticks([2e-12, 1e-11, 1e-10])

    # Plot the three hypothesis from the text into the three axes
    # The parameters below are taken from Table 1 in the text
    
    print 'TODO: ask Yasunobu about how time dependence is implemented (P10 of text)!!!'

    sed(grid[0], type='pi0',
        label='(a) Pion-decay dominated',
        delta_s = 1.4,
        e_break = 15*u.GeV,
        magnetic_field = 40*u.microgauss,
        hydrogen_density = 10*u.cm**-3,
        Wp = 5.2e50*u.erg,
        We = 0.13e50*u.erg)

    sed(grid[1], type='brems',
        label='(b) Brems dominated',
        delta_s = 1.4,
        e_break = 5*u.GeV,
        magnetic_field = 15*u.microgauss,
        hydrogen_density = 10*u.cm**-3,
        Wp = 0.54e50*u.erg,
        We = 0.87e50*u.erg)

    sed(grid[2], type='ic',
        label='(c) IC dominated',
        delta_s = 2.3,
        e_break = 20*u.GeV,
        magnetic_field = 2*u.microgauss,
        hydrogen_density = 0.1*u.cm**-3,
        Wp = 8.4e50*u.erg,
        We = 11e50*u.erg)


    fig.savefig('w51C_sed.pdf')
