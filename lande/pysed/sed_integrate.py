""" These functions define integration routines 
    necessary through the code.

    Author: Joshua Lande <joshualande@gmail.com>
"""
import numpy as np
from scipy import integrate

from . import units as u
from . helper import logrange

def dbltrapz(f, x, y):
    """ Perform a double trapezon integration of a vectorized function

        f(x,y) where x is an array of x values and y is an array of y 
        values to perform the integral over.
    
        Implementation taken from: 
            http://mail.scipy.org/pipermail/scipy-user/2011-February/028592.html

        I don't really understand it, but it seems to work...
        
        For example, if we wanted to integrate some really random function 

            f(x,y) = e^(-x^2*y)

        from 1<x<5 and 1<y<10.

        Using matheamtica, we find that the integral is ~0.09:
        
            N[Integrate[Integrate[Exp[-(x^2*y)], {x, 1, 5}], {y, 1, 10}], 20]

        Using our new function.

        >>> import numpy as np
        >>> f=lambda x,y: np.exp(-(x**2*y))
        >>> int=dbltrapz(f, np.linspace(1,5,1000), np.linspace(1,10,1000))
        >>> print np.allclose(int,0.089071862226039234609, rtol=1e-5, atol=1e-5)
        True
    """
    yy = y[:,np.newaxis]
    xx = x[np.newaxis,:]
    integrand=f(xx,yy)

    return integrate.trapz(integrate.trapz(integrand, yy, axis=0), x, axis=0)

def dblsimps(f, x, y):
    """ Same as dbltrapz, but performs a 2D simpson integral. 
        
        It should work the same as the other function:

        >>> import numpy as np
        >>> f=lambda x,y: np.exp(-(x**2*y))
        >>> int=dbltrapz(f, np.linspace(1,5,1000), np.linspace(1,10,1000))
        >>> print np.allclose(int,0.089071862226039234609, rtol=1e-5, atol=1e-5)
        True
    """
    yy = y[:,np.newaxis]
    xx = x[np.newaxis,:]
    integrand=f(xx,yy)
    return integrate.simps(integrate.simps(integrand, yy, axis=0), x, axis=0)


def logsimps(f,xmin,xmax, per_decade):
    """ Perform the simpson integral of a function f(x)
        from xmin to xmax evaluationg the function
        uniformly in log space.

        This integration method is useful for typical astrophysical spectral which
        roughly fairly constant over many decades in energy

        per_decade is the number of points per decade in energy
        to evaluate the function over.

        Implementation Note: int f(x) dx = int f(x) x dlog(x). 

        For example, a typical function to integrate is x^-2. To integrate it from 1e-2 to 1e2

        Using mathematica, we find that the integral is ~100:

            N[Integrate[x^-2, {x, .01, 100}]]

        Using our nice function:

            >>> print np.allclose(logsimps(lambda x: x**-2, .01, 100, 1000), 99.989, rtol=1e-5, atol=1e-5)
            True
    """
    x = logrange(xmin, xmax, per_decade)
    y = f(x) * x
    log_x = np.log(x)
    return integrate.simps(y=y, x=log_x)


def dbllogsimps(f,xmin,xmax, ymin, ymax, per_decade):
    """ Perform a simpson integral of f(x,y) where
        both x and y are sampled uniformly in log space.

        Implementation Note: int f(x,y) dx dy = int f(x,y) x*y*dlog(x)*dlog(y)
    """
    x = logrange(xmin, xmax, per_decade)
    y = logrange(ymin, ymax, per_decade)
    integrand = lambda x,y: f(np.exp(x),np.exp(y)) * np.exp(x) * np.exp(y)
    log_x, log_y = np.log(x), np.log(y)
    return dblsimps(f=integrand, x=log_x, y=log_y)

def halfdbllogsimps(f, xmin, xmax, ymin, ymax, x_per_decade, y_npts):
    """ Perform a simpson integral of f(x,y) where
        x is sampled uniformly in log space while y is
        sampled uniformly.

        Implementation Note: int f(x,y) dx dy = int f(x,y) x*dlog(x)*dy
    """
    x = logrange(xmin, xmax, x_per_decade)
    y = np.linspace(ymin, ymax, y_npts)
    integrand = lambda x,y: f(np.exp(x),y)*np.exp(x)
    log_x = np.log(x)
    return dblsimps(f=integrand, x=log_x, y=y)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
