""" Function for generating random quantities.

"""
import math
import numpy as np
from skymaps import SkyDir


def random_on_sphere():
    """ Pick a random SkyDir on the sphere. """
    l=np.random.uniform(0,360)
    b=180*math.acos(np.random.uniform(-1,1))/math.pi-90
    skydir=SkyDir(l,b,SkyDir.GALACTIC)
    return skydir

def random_point_near(skydir,distance):
    """ Generate a random point within a given distance of another point """
    while True:
        x = np.random.uniform(-distance,distance)
        y = np.random.uniform(-distance,distance)
        if np.sqrt(x**2+y**2)<distance:
            break

    if skydir.ra()<0.1 and skydir.dec()<0.1:
        # Rotation fails in this case.
        return SkyDir(x+skydir.ra(),y+skydir.dec())
    else:
        rotated_dir=SkyDir(x,y)
        return DualLocalizer.anti_rotate_equator(rotated_dir,skydir)

def mixed_linear(min,max,num):
    """ Just like np.linspace but the numbers are mixed up
        using the Van der Corput sequence. Handy for getting
        a reasonable sample quickly.
        
        Use the http://en.wikipedia.org/wiki/Van_der_Corput_sequence
        to mix up the numbers. 

        To get the csc package, go to http://pypi.python.org/pypi/csc-utils
        
        """
    from csc import util
    x=util.sampling_sequence(min,max)
    return np.asarray([x.next() for i in range(num)])

def mixed_log(min,max,num):
    return 10**mixed_linear(np.log10(min),np.log10(max),num)

