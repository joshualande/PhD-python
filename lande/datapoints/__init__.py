"""

    TODO: implement lower limits
"""
import re
import numpy as np
import math
import numbers
from uncertainties import ufloat

class DataPoint(object):
    
    @staticmethod
    def from_string(string):
        """

            >>> f = DataPoint.from_string
            >>> f('$0.1 \pm 0.2$')
            0.1 +/- 0.2
            >>> f('$<0.3$')
            <0.3
            >>> f('0.1_{0.1}^{0.2}')
            1.0 +/- 0.15
            >>> f('0.1^{0.1}_{0.2}')
            1.0 +/- 0.15
            >>> f('0.2')
            0.2

        """
        if '$' in string:
            string = string.replace('$','')
        if r'\,' in string:
            string = string.replace(r'\,','')

        string = string.strip()


        m = re.match( r'(.+)\+\/\-(.+)', string)
        if m:
            value,error = m.group(1), m.group(2)
            return Detection(float(value),float(error))
        else:
            m = re.match( r'(.+)\\pm(.+)', string)
            if m:
                value,error = m.group(1), m.group(2)
                return Detection(float(value),float(error))
            else:
                m = re.match( r'.*<(.+)', string)
                if m:
                    limit = m.group(1)
                    return UpperLimit(float(limit))
                else:
                    m = re.match( r'(.+)_\{(.+)\}\^\{(.+)\}.*', string)
                    if m:
                        value, lower, upper = map(float,[m.group(1), m.group(2), m.group(3)])
                        return AssymetricError(value, lower=abs(lower), upper=abs(upper))
                    else:
                        m = re.match( r'(.+)\^{(.+)\}\_\{(.+)\}.*', string)
                        if m:
                            value, lower, upper = map(float,[m.group(1), m.group(2), m.group(3)])
                            return AssymetricError(value, lower=abs(lower), upper=abs(upper))
                        else:
                            try:
                                return float(string)
                            except:
                                raise Exception("Unrecognized string: %s" % string)

class AssymetricError(DataPoint):
    def __init__(self,value, lower, upper):
        self.value,self.lower,self.upper = value,lower,upper
    def __repr__(self):
        return '%s +%s -%s' % (self.value,self.lower,self.upper)
    def _repr_html_(self):
        return '$%g_{-%g}^{+%g}$' % (self.value, self.lower,self.upper)
    @staticmethod
    def significant(): return True
            
class Detection(DataPoint):
    def __init__(self,value, error):
        self.value,self.error = value,error
    def __repr__(self):
        return '%s +/- %s' % (self.value,self.error)
    def _repr_html_(self):
        return '$%g \pm %g $' % (self.value, self.error)
    @staticmethod
    def significant(): return True

    def to_uncertainty(self):
        return ufloat((self.value,self.error))

    @staticmethod
    def from_uncertainty(ufloat):
        return Detection(ufloat.nominal_value, ufloat.std_dev())

    def __pow__(self, other):
        if isinstance(other, numbers.Number):
            if np.isnan(other):
                return float('nan')
            else:
                return Detection.from_uncertainty(self.to_uncertainty().__pow__(other))
        elif isinstance(other, Detection):
            return Detection.from_uncertainty(self.to_uncertainty().__pow__(other.to_uncertainty()))
        else:
            raise NotImplemented()

    def __div__(self, other):
        if isinstance(other, numbers.Number):
            if np.isnan(other):
                return float('nan')
            else:
                return Detection.from_uncertainty(self.to_uncertainty().__div__(other))
        elif isinstance(other, float):
            return Detection(self.value/other,self.error/other)
        elif isinstance(other, Detection):
            return Detection.from_uncertainty(self.to_uncertainty().__div__(other.to_uncertainty()))
        elif isinstance(other, UpperLimit):
            return LowerLimit(self.value/other.upper_limit)
        else:
            raise NotImplemented()

    def __mul__(self, other):
        if isinstance(other, numbers.Number):
            if np.isnan(other):
                return float('nan')
            else:
                return Detection.from_uncertainty(self.to_uncertainty().__mul__(other))
        elif isinstance(other, Detection):
            return Detection.from_uncertainty(self.to_uncertainty().__mul__(other.to_uncertainty()))
        elif isinstance(other, UpperLimit):
            return UpperLimit(self.value*other.upper_limit)
        elif isinstance(other, LowerLimit):
            return LowerLimit(self.value*other.lower_limit)
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
            return UpperLimit(self.value+other.upper_limit)
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
    @staticmethod
    def significant(): return False

    def __div__(self, other):
        if isinstance(other, numbers.Number):
            if np.isnan(other):
                return float('nan')
            return UpperLimit(self.upper_limit/other)
        elif isinstance(other, Detection):
            return UpperLimit(self.upper_limit/(other.value))
        elif isinstance(other, Detection):
            raise NotImplemented()

    def __mul__(self, other):
        if isinstance(other, numbers.Number):
            if np.isnan(other):
                return float('nan')
            else:
                return UpperLimit(self.upper_limit*other)
        else:
            raise NotImplemented()
        
    def __repr__(self):
        return '<%s' % self.upper_limit

    def _repr_html_(self):
        return '$<%g$' % self.upper_limit
        
    def __eq__(self, other):
        return self.upper_limit == other.upper_limit

class LowerLimit(DataPoint):
    def __init__(self,lower_limit):
        self.lower_limit = lower_limit
    @staticmethod
    def significant(): return False

    def __repr__(self):
        return '>%s' % self.lower_limit

    def _repr_html_(self):
        return '$>%g$' % self.lower_limit
        
    def __eq__(self, other):
        return self.lower_limit == other.lower_limit

def plot(x, y, 
         axes=None, 
         ul_fraction=0.4,
         log_clipping=False, **kwargs):

    import pylab as P
    import matplotlib.lines as mlines

    plot_kwargs = dict(linestyle='none', capsize=0)
    plot_kwargs.update(kwargs)

    ul_kwargs = plot_kwargs.copy()
    for k in ['capsize', 'elinewidth', 'marker']:
        if k in ul_kwargs: ul_kwargs.pop(k)

    if axes is None:
        axes = P.gca()

    for _x, _y in zip(x,y):
        if isinstance(_y,Detection):

            axes.errorbar([_x], [_y.value],
                          yerr=[_y.error],
                          **plot_kwargs)

        elif isinstance(_y,AssymetricError):
            axes.errorbar([_x], [_y.value],
                          yerr=[[_y.lower], [_y.upper]],
                          **plot_kwargs)

        elif isinstance(_y,UpperLimit):

            axes.errorbar([_x], [_y.upper_limit],
                          yerr=[ [ul_fraction*_y.upper_limit], [0] ],
                          **plot_kwargs)

            axes.plot([_x], [(1-ul_fraction)*_y.upper_limit],
                      marker=mlines.CARETDOWN,
                      **ul_kwargs)

        elif isinstance(_y,float):

            axes.errorbar([_x], [_y],
                          **plot_kwargs)

        elif isinstance(_y, LowerLimit):

            axes.errorbar([_x], [_y.lower_limit],
                          yerr=[ [0], [ul_fraction*_y.lower_limit] ],
                          **plot_kwargs)

            axes.plot([_x], [(1+ul_fraction)*_y.lower_limit],
                      marker=mlines.CARETUP,
                      **ul_kwargs)

        else:
            raise Exception("...")

if __name__ == "__main__":
    import doctest
    doctest.testmod()
