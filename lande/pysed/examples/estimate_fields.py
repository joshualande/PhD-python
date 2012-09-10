""" An example to demonstrate the way in
    which the Interstellar Radiation Field (ISRF)
    can be estiamted at a given place in the galaxy
    from the predictions by GALPROP.

    Typical reference for GALPROP: Moskalenko et al 2006 (arXiv:astro-ph/0511149v2)

    Author: Joshua Lande <joshualande@gmail.com>
"""
import pylab as P
from pysed.sed_isrf import ISRF
from pysed.sed_spectrum import CompositeSpectrum
import pysed.units as u

if __name__ == '__main__':

    # load in ISRF
    isrf = ISRF('MilkyWay_DR0.5_DZ0.1_DPHI10_RMAX20_ZMAX5_galprop_format.fits.gz')

    # define position in galaxy
    kwargs=dict(R=5*u.kpc,z=0)

    # get data form mapcube
    energy = isrf.get_energy()
    infrared=isrf.get_infrared(**kwargs)
    cmb=isrf.get_CMB(**kwargs)
    optical=isrf.get_optical(**kwargs)

    # convert to unitless quanitites + plot
    energy = u.tonumpy(energy, u.eV)
    infrared = u.tonumpy(infrared, u.cm**-3*u.eV**-1)
    cmb = u.tonumpy(cmb, u.cm**-3*u.eV**-1)
    optical = u.tonumpy(optical, u.cm**-3*u.eV**-1)

    plot_kwargs = dict(marker='+', linestyle='none')
    P.loglog(energy, infrared, color='red', label='infrared', **plot_kwargs)
    P.loglog(energy, cmb, color='blue', label='CMB', **plot_kwargs)
    P.loglog(energy, optical, color='green', label='optical', **plot_kwargs)

    P.loglog(energy, infrared + cmb + optical, color='black', label='total', **plot_kwargs)

    axes = P.gca()


    # now, overlay estimated quantities
    plot_kwargs = dict(axes=axes, x_units_string = 'eV', y_units_string='ph/cm^3/eV', dashes=[5,2])

    infrared_est = isrf.estimate_infrared(**kwargs)
    infrared_est.loglog(color='red', **plot_kwargs)

    CMB_est = isrf.estimate_CMB(**kwargs)
    CMB_est.loglog(color='blue', **plot_kwargs)
    print 'CMB kT',u.repr(CMB_est.kT*u.erg/u.boltzmann,'kelvin')

    optical_est = isrf.estimate_optical(**kwargs)
    optical_est.loglog(color='green', **plot_kwargs)

    total = CompositeSpectrum(infrared_est, CMB_est, optical_est)
    total.loglog(color='black', **plot_kwargs)

    P.legend(loc=3)

    P.xlabel('energy (eV)')
    P.ylabel('intensity (ph/cm^3/eV)')

    P.savefig('estimate_fields.pdf')
