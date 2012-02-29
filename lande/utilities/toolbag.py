""" This file contains various function which I have found useful, but
    can't think of a better place to put. """
from collections import OrderedDict
import numpy as np
import pyfits as pf
from uw.like.roi_analysis import ROIAnalysis
from uw.like.roi_extended import ExtendedSource
from uw.pulsar.phase_range import PhaseRange


def tolist(x):
    """ convenience function that takes in a 
        nested structure of lists and dictionaries
        and converst all

        (a) numpy arrays into python lists
        (b) numpy strings into python scrings.
        (c) an ordered dict to a dict

        which is conveneint for duming a file to yaml.
    """
    if isinstance(x,list):
        return map(tolist,x)
    elif isinstance(x,dict):
        return dict((tolist(k),tolist(v)) for k,v in x.items())
    elif isinstance(x,np.ndarray) or \
            isinstance(x,np.number):
        return x.tolist()
    elif isinstance(x,np.str):
        return str(x)
    elif isinstance(x,PhaseRange):
        return x.tolist(dense=True)
    elif isinstance(x,OrderedDict):
        return dict(x)
    else:
        return x


def expand_fits(pf,factor,hdu=0):
    """ Create a new fits file where
        each pixel is divided into factor**2
        more pixels all of equal value. 
        
        I am not an expert on Fits files and
        I only think this code works. Also,
        I have no idea if it would work with 3D data (?).
        As a rule of thumb, if you use this function,
        you should probably open the before/after
        in ds9 and blink them for a while to convince
        yourself the function worked.
    """

    h=pf[hdu].header
    d=pf[hdu].data

    h['CDELT1']/=factor
    h['CDELT2']/=factor

    h['NAXIS1']*=factor
    h['NAXIS2']*=factor

    h['CRPIX1']=h['CRPIX1']*factor - factor/2.0 + 0.5
    h['CRPIX2']=h['CRPIX2']*factor - factor/2.0 + 0.5

    larger=list(d.shape)
    larger[0]*=factor
    larger[1]*=factor
    #larger_array = np.empty(larger,dtype=d.dtype)
    larger_array = np.zeros(larger,dtype=d.dtype)
    for i in range(factor):
        for j in range(factor):
            larger_array[i::factor,j::factor] = d

    pf[hdu].data=larger_array



# Taken form http://stackoverflow.com/questions/4126348/how-do-i-rewrite-this-function-to-implement-ordereddict
class OrderedDefaultdict(OrderedDict):
    def __init__(self, *args, **kwargs):
        newdefault = None
        newargs = ()
        if len(args):
            newdefault = args[0]
            if not callable(newdefault) and newdefault != None:
                raise TypeError('first argument must be callable or None')
            newargs = args[1:]
        self.default_factory = newdefault
        super(OrderedDefaultdict, self).__init__(*newargs, **kwargs)

    def __missing__ (self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = value = self.default_factory()
        return value




from datetime import datetime
from dateutil import relativedelta
class StopWatch(object):
    """ Simple object for timing code"""
    def __init__(self):
        self.start_time = self.__lap = datetime.now()
  
    @property
    def lap(self):
        current_time = datetime.now()
        lap = relativedelta.relativedelta(current_time,self.__lap)
        self.__lap = current_time
        return lap

    @property
    def total(self):
        return datetime.now() - self.start_time
    
    def __str__(self):
        ret=[]
        lap=self.lap
        if lap.years > 0: ret.append('%g years' % lap.years)
        if lap.months > 0: ret.append('%g months' % lap.months)
        if lap.days > 0: ret.append('%g days' % lap.days)
        if lap.hours > 0: ret.append('%g hours' % lap.hours)
        if lap.minutes > 0: ret.append('%g minutes' % lap.minutes)
        if lap.seconds > 0: ret.append('%g seconds' % lap.seconds)

        # only add microseconds if necessary
        if len(ret) < 1: ret.append('%g microseconds' % lap.microseconds)

        return 'Lap time is '+', '.join(ret[0:2]) # only print out leading two 



