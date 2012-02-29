""" Various helper functions I can't think
    of a better place to put.

    Author: Joshua Lande <joshualande@gmail.com>
"""
import numpy as np
from . import sed_units as u

def logrange(min,max,per_decade):
    """ Creates a range of values from min to max
        with per_decade points per logarithmic
        decade of energy. """
    npts = int(np.ceil(per_decade*(np.log10(max)-np.log10(min))))
    x = np.logspace(np.log10(min),np.log10(max), npts+1)
    return x

def linspace_unit(min, max, npts):
    """ Convenience function to compute
        an array of points between
        min and max when min and max
        are united and have the same units.
        
        Example:

        >>> linspace_unit(1*u.cm,4*u.cm,4)
        [0.01*m, 0.02*m, 0.03*m, 0.04*m]
    """
    # make sure numbers have same units
    val = lambda x: x.as_two_terms()[0]
    unit = lambda x: x.as_two_terms()[1]

    assert(unit(min)==unit(max))

    units = unit(min)
    min,max = val(min), val(max)

    return u.tosympy(linspace(float(min),float(max), npts), units)

def logrange_unit(min,max, per_decade):
    """ Convenience function to compute
        an array of points uniformly
        between
        min and max when min and max
        are united and have the same units.
        
        Example:

        >>> logrange_unit(1*u.cm,1e3*u.cm,1)
        [0.01*m, 0.1*m, 1.0*m, 10.0*m]
    """
    # make sure numbers have same units
    val = lambda x: x.as_two_terms()[0]
    unit = lambda x: x.as_two_terms()[1]

    assert(unit(min)==unit(max))

    units = unit(min)
    min,max = val(min), val(max)

    return u.tosympy(logrange(float(min),float(max), per_decade), units)

def argmax_unit(array):
    """ Computes the argmax of a sympy list of numbers with units. 
    
        Example:

        >>> import sympy
        >>> x = sympy.Matrix([1*u.erg, 10*u.erg, 2*u.erg])
        >> argmax_unit(x)
        1
    """
    return np.argmax(array.tolist())

if __name__ == "__main__":
    import doctest
    doctest.testmod()
