""" Code relateded to defining a base spectrum
    class that can be subclassed.

    Author: Joshua Lande <joshualande@gmail.com>
"""
from abc import abstractmethod
from operator import add

import numpy as np
import pylab as P
import sympy

from . import sed_config
from . import units as u
from . helper import logrange,logrange_unit

class Spectrum(object):
    """ A base class which represents some
        physical quanity as a function of energy. 

        All spectra implement a __call__ function
        which can be used to obtain the spectrum
        as a function of energy.
        
        Implementation notes:
            * All function must implement a function
              units_string() which is a string specifying
              the units of the returned spectra

            * All spectra must define the function
              '_spectrum' which takes in energy in units of
              erg and returns the output represented
              in the same units as units_string()

            * If the spectrum function is vectorized
              so that a numpy array of energies can 
              be passed into the spectrum function
              and a numpy array of spectra
              are returned, the value self.vectorized
              must be set to True. Otherwise, the
              value self.vectorized must be set to False.

            * All spectra must define an energy
              range over which the spectra is non-zero.
              This energy range is stored internall
              as self.emin and self.emax and must
              be in units of erg. It is not important
              if the spectrum function evalulated
              at energies outside of the range returns
              non-zero values, since the __call__
              function will inforce these limits for
              all spectra.
    """

    def __call__(self, energy, units=True):
        """ Returns number of particles per unit energy [1/energy]. 
        
            This function vectorized the output if a numpy array
            or Sympy Matix of energies is passed """
        def nonzero(energy):
            return (((energy>=self.emin) if hasattr(self,'emin') else True)&\
                    ((energy<=self.emax) if hasattr(self,'emax') else True))

        if isinstance(energy,np.ndarray) and units==False:
            if self.vectorized:
                return np.where(nonzero(energy),self._spectrum(energy),0)
            else:
                return np.asarray([self(i) for i in energy])

        if isinstance(energy,sympy.Matrix) and units==True:
            if self.vectorized:
                energy = u.tonumpy(energy,u.erg)
                spectrum = np.where(nonzero(energy),self._spectrum(energy),0)
                return u.tosympy(spectrum,self.units())
            else:
                return sympy.Matrix([self(i) for i in energy]).transpose()

        if units: energy = float(energy/u.erg)
        spectrum=np.where(nonzero(energy), self._spectrum(energy), 0)
        return spectrum*(self.units() if units else 1)


    @abstractmethod
    def _spectrum(self,energy):
        raise NotImplementedError("Subclasses must implement this.")

    @abstractmethod
    def units_string(self):
        pass

    def units(self):
        """ Returns the units that __call__ is assumed to be in. """                                                                                                        
        return u.fromstring(self.units_string())                                                                                                                             

    def loglog(self, 
               x_units_string,
               y_units_string,
               emin=None, emax=None, 
               e_weight=0, # weight the function by energy raised to this power
               scale=1, # scale
               x_label=None,
               y_label=None,
               filename=None, fignum=None, 
               axes=None, **kwargs):
        """ Plots the energy spectrum. 

            emin and emax be in energy units.
            
            x_units_string and y_units_string must be strings suitable
            for plotting on the matplotlib axes. """

        if axes is None:
            fig = P.figure(fignum,figsize=(5.5,4.5))
            P.clf()
            fig.subplots_adjust(left=0.18,bottom=0.13,right=0.95,top=0.95)
            axes = fig.add_subplot(111)

            if x_label is None: x_label = x_units_string
            if y_label is None: y_label = ('E$^%s$' % e_weight if e_weight>0 else '') + 'dN/dE (%s)' % y_units_string

            axes.set_xlabel(x_label)
            axes.set_ylabel(y_label)

        if emin is None: 
            if not hasattr(self,'emin'): raise Exception("Emin must be set.")
            emin=self.emin*u.erg
        if emax is None: 
            if not hasattr(self,'emax'): raise Exception("Emax must be set.")
            emax=self.emax*u.erg

        x = logrange_unit(emin,emax,sed_config.PER_DECADE)
        y = scale*self(x)
        for i in range(e_weight): y=y.multiply_elementwise(x)

        x=u.tonumpy(x,u.fromstring(x_units_string))
        y=u.tonumpy(y,u.fromstring(y_units_string))

        axes.loglog(x,y, **kwargs)

        if filename is not None: P.savefig(filename)
        return axes


class NumericalSpectrum(Spectrum):
    """ Takes an analytic spectrum but 
        stores it internally as a numerical array of values.

        PDF is either evaluated at a user specified list
        of energies or uniformly
        
        a parameter per_decade determines how many points to
        store the function at per decade in energy.

        Function is evaluated by interpolating in log-space.

        Useful for non-analytic spectra, such as calcluating
        energy losses of a population over time.

        Could also be useful for a very costly spectrum
        to be cached for repeated evaluation.
    """
    vectorized = True

    def __init__(self, spectrum, energies=None, per_decade=None):


        if energies is not None:
            self.energies = energies
        else:
            if per_decade is None: per_decade = sed_config.PER_DECADE

            if not hasattr(spectrum,'emin') or not hasattr(spectrum,'emax'):
                raise Exception("NumericalSpectrum must be passed a spectrum with a min and max energy.")

            self.energies = logrange(spectrum, emin, emax, self.per_decade)


        self.values = spectrum(self.energies, units=False)

        self.log_energies = np.log10(self.energies)
        self.log_values = np.log10(self.log_values)

        self.log_interp = interp1d(self.log_energies,self.log_values,kind=kind, bounds_error=True)


        self._units_string = specturm.units_string()

    def _spectrum(self,energy):
        return 10**self.log_interp(energy)

    def units_string(self): return self._units_string

class Constant(Spectrum):

    vectorized = True

    def __init__(self, value, units_string):
        """ Implements a constant value over a range of energy. 
            
            Value is the constant value to return
            unit_string is the units to store the value internall with

            For example:

                >>> constant = Constant(u.erg, 'erg')
                >>> print u.repr(constant(u.MeV),'erg','%g')
                1 erg
        """
        self._units_string = units_string
        self.internal = float(value/self.units())
        self.value = value

    def _spectrum(self,energy):
        return self.internal*np.ones_like(energy)

    def units_string(self): return self._units_string


class CompositeSpectrum(Spectrum):
    """ This class represents a linear combination
        of Spectrum objects.

        Example:

            >>> c1 = Constant(u.erg,     'erg')
            >>> c2 = Constant(0.5*u.erg, 'erg')
            >>> c  = CompositeSpectrum(c1,c2)

        The composite spectrum is the sum of the two spectra (1.5 erg)

            >>> print u.repr(c(u.MeV),'erg','%g')
            1.5 erg

    """
    @staticmethod
    def all_same(items): return len(set(items))==1

    def __init__(self, *spectra):

        if not self.all_same([s.units_string() for s in spectra]):
            raise Exception('Error in CompositeSpectrum: all spectra must have the same units.')
        self._units_string = spectra[0].units_string()
        self._units = spectra[0].units()

        self.spectra = spectra
        self.emin = min([i.emin for i in spectra])
        self.emax = max([i.emax for i in spectra])
        self.vectorized = np.all([i.vectorized for i in spectra])
    
    def __call__(self, *args, **kwargs):
        """ Nb, override the __call__ function instead of the spectrum
            object in case the emin-emax energy ranges are inconsistent
            for the different spectra. """
        return reduce(add,[s(*args,**kwargs) for s in self.spectra])

    def units_string(self): return self._units_string

if __name__ == "__main__":
    import doctest
    doctest.testmod()
