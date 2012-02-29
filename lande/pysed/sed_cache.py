""" Code to cache mathematical functions for faster evaluation.

    Author: Joshua Lande <joshualande@gmail.com>
"""
import numpy as np
from scipy.interpolate import interp1d

class FunctionCache(object):
    """ Simple object to Cache a function f between xmin and xmax.

        f: function to cache
        xmin, xmax: range over which to cache function
        npts: number of points in which to cache the funciton
        kind: type of function interpolation. See scipy.interpolate.interpt1d
        bound_error: raise an exception if function evaluated outside of xmin-xmax
        fill_value: default value outside range. 
        
        Example:
        
            >>> f = lambda x: x**2
            >>> F = FunctionCache(f, 0,10) 
            >>> x = np.asarray([0,2,4,8])
            >>> print np.allclose(F(x), x**2)
            True
    """

    def __init__(self, f, xmin, xmax, npts=1000, kind='linear', bounds_error=False, fill_value=0):
        self.F = f
        self.x = np.linspace(xmin,xmax,npts)
        self.y = np.asarray([f(i) for i in self.x])

        self.interp=interp1d(self.x,self.y,kind=kind,bounds_error=bounds_error,fill_value=fill_value)

    def __call__(self,x):
        return self.interp(x)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
