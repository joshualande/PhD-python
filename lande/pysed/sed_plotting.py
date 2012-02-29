""" Code related to plotting SEDs.

    Author: Joshua Lande <joshualande@gmail.com>
"""
import numpy as np
import pylab as P

import sympy

from . import units as u

class SEDPlotter(object):
    """ Plot an E^2 dN/dE SED with multiple componets.
    
        Spectra is a dictionary:

        E.G:
            spectra = {'Synchrotron': sed_synch.Synctrotron(...),
                       'Inverse Compton': sed_ic.InverseCompton(...)}
    """

    def __init__(self, 
                 distance,
                 x_units_string = 'eV',
                 y_units_string = 'eV/cm^2/s^1',
                 emin=None,
                 emax=None,
                 axes=None, 
                 fignum=None, 
                 figsize=(5.5,4.5)):

        self.distance=distance
        self.x_units_string = x_units_string
        self.x_units = u.fromstring(x_units_string)
        self.y_units_string = y_units_string 
        self.y_units = u.fromstring(y_units_string)
        self.emin, self.emax = emin, emax

        self.scale = 1/(4*np.pi*distance**2)

        # Here, define the axes
        if axes is None:
            fig = P.figure(fignum,figsize)
            P.clf()
            fig.subplots_adjust(left=0.18,bottom=0.13,right=0.95,top=0.95)
            self.axes = fig.add_subplot(111)

            # Format the axes
            self.axes.set_xlabel('Energy (%s)' % x_units_string)
            self.axes.set_ylabel(r'E$^2$ dN/dE (%s)' % y_units_string)

            if self.emin is None or self.emax is None:
                raise Exception("Either an existing axes must be passed into the class or emin and emax must both be set.")

            self.axes.set_xlim(xmin=float(self.emin/self.x_units), xmax=float(self.emax/self.x_units))
        else:
            self.axes = axes
            if self.emin is None and self.emax is None:
                self.emin, self.emax = [i*self.x_units for i in self.axes.get_xlim()]

    def plot(self, spectra, **kwargs):
        spectra.loglog(emin=self.emin, emax=self.emax,
                 x_units_string = self.x_units_string,
                 y_units_string = self.y_units_string,
                 e_weight=2,
                 scale=self.scale,
                 axes=self.axes,
                 **kwargs)

    def save(self, filename):
        self.axes.figure.savefig(filename)
