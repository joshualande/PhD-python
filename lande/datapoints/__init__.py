"""

    TODO: implement lower limits
"""
import numpy as np
import math
import numbers
from uncertainties import ufloat

class DataPoint(object):
    
    @staticmethod
    def from_string(string):
        string = string.strip()
        if string[0] == '$':
            string = string[1:]
        if string[-1] == '$':
            string = string[:-1]

        if '+/-' in string:
            value,error = string.split('+/-')
            return Detection(float(value),float(error))
        elif '\pm' in string:
            value,error = string.split('\pm')
            return Detection(float(value),float(error))

        elif string[0] == '<':
            return UpperLimit(float(string[1:]))
        else:
            raise Exception("Unrecognized string: %s" % string)
            
class Detection(DataPoint):
    def __init__(self,value, error):
        self.value,self.error = value,error
    def __repr__(self):
        return '%s +/- %s' % (self.value,self.error)
    def _repr_html_(self):
        return '$%s \pm %s $' % (self.value, self.error)

    def to_uncertainty(self):
        return ufloat((self.value,self.error))

    @staticmethod
    def from_uncertainty(ufloat):
        return Detection(ufloat.nominal_value, ufloat.std_dev())

    def __pow__(self, other):
        if isinstance(other, numbers.Number):
            if np.isnan(other):
                return float('nan')
            return Detection.from_uncertainty(self.to_uncertainty().__pow__(other))
        elif isinstance(other, Detection):
            return Detection.from_uncertainty(self.to_uncertainty().__pow__(other.to_uncertainty()))
        else:
            raise NotImplemented()

    def __div__(self, other):
        if isinstance(other, numbers.Number):
            if np.isnan(other):
                return float('nan')
            return Detection.from_uncertainty(self.to_uncertainty().__div__(other))
        elif isinstance(other, Detection):
            return Detection.from_uncertainty(self.to_uncertainty().__div__(other.to_uncertainty()))
        else:
            raise NotImplemented()

    def __mult__(self, other):
        if isinstance(other, numbers.Number):
            if np.isnan(other):
                return float('nan')
            return Detection.from_uncertainty(self.to_uncertainty().__mult__(other))
        elif isinstance(other, Detection):
            return Detection.from_uncertainty(self.to_uncertainty().__mult__(other.to_uncertainty()))
        elif isinstance(other, UpperLimit):
            return UpperLimit((self.value + self.error)*other.upper_limit)
        else:
            raise NotImplemented()


    def __add__(self, other):
        if isinstance(other, numbers.Number):
            if np.isnan(other):
                return float('nan')
            return Detection.from_uncertainty(self.to_uncertainty().__add__(other))
        elif isinstance(other, Detection):
            return Detection.from_uncertainty(self.to_uncertainty().__add__(other.to_uncertainty()))
        elif isinstance(other, UpperLimit):
            return UpperLimit(self.value + self.error +other.upper_limit)
        else:
            raise NotImplemented()
    
    def __radd__(self, other):
        print self, other
        self.__add__(other)

    def __eq__(self, other):
         return self.value == other.value and self.error == other.error
        

class UpperLimit(DataPoint):
    def __init__(self,upper_limit):
        self.upper_limit = upper_limit
    def __add__(self, other):
        if isinstance(other, UpperLimit):
            return UpperLimit(self.upper_limit + other.upper_limit)
        elif isinstance(other, Detection):
            return other.__add__(self) 
        else:
            raise NotImplemented()

    def __div__(self, other):
        if isinstance(other, numbers.Number):
            if np.isnan(other):
                return float('nan')
            return UpperLimit(self.upper_limit/other)
        elif isinstance(other, Detection):
            return UpperLimit(self.upper_limit/(other.value - other.error))
        elif isinstance(other, Detection):
            raise NotImplemented()
        
    def __repr__(self):
        return '<%s' % self.upper_limit

    def _repr_html_(self):
        return '$<%s$' % self.upper_limit
        
    def __eq__(self, other):
        return self.upper_limit == other.upper_limit
        

def plot(x, y, xlo=None, xhi=None, **kwargs):
    from lande.utilities.plotting import plot_points

    _y = np.asarray([i.value if isinstance(i,Detection) else np.nan for i in y])
    _y_lower_err = _y_upper_err = np.asarray([i.error if isinstance(i,Detection) else np.nan for i in y])
    _y_ul = np.asarray([i.upper_limit if isinstance(i,UpperLimit) else np.nan for i in y])
    _significant = np.asarray([True  if isinstance(i,Detection) else False for i in y])

    print _y,_y_lower_err,_y_ul,_significant
    
    plot_points(x, _y, xlo=xlo, xhi=xhi, 
                y_lower_err=_y_lower_err, y_upper_err=_y_upper_err, 
                y_ul=_y_ul, significant=_significant,
                **kwargs)




