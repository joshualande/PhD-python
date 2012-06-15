import copy
import sys
from argparse import ArgumentParser
from collections import OrderedDict

from StringIO import StringIO
import asciitable

from . tools import parse_strip_known_args

def get_confluence():
    parser = ArgumentParser()
    parser.add_argument("--confluence", action="store_true", default=False)
    args=parse_strip_known_args(parser)
    return args.confluence


def confluence_table(table_dict, units=None, **kwargs):
    """ Format a dictionary into a confluence table. 
        Use the DefaultDict if you care about the order of the table.
    
        >>> d=OrderedDict()
        >>> d['name']=['bright','dim']
        >>> d['flux']=['large','small']
        >>> print confluence_table(d)
        ||   name ||  flux ||
        || bright  | large  |
        ||    dim  | small  |

        >>> print confluence_table(d, units=dict(flux='[ergs]'))
        ||   name || flux [ergs] ||
        || bright  |       large  |
        ||    dim  |       small  |


    """

    if units is not None:
        for k,v in units.items():
            table_dict[k + ' %s' % v] = table_dict.pop(k)

    outtable=StringIO()

    asciitable.write(table_dict, outtable,
                 Writer=asciitable.FixedWidth,
                 names=table_dict.keys(),
                 **kwargs)
    t=outtable.getvalue()
    t=t.replace(' |','  |')
    t=t.strip().split('\n')
    t[0]=t[0].replace(' |','||')
    t=['|'+i for i in t]
    return '\n'.join(t)


def latex_table(table_dict, units=None, latexdict=None, **kwargs):
    r""" Format a dictionary into a latex table.
         Use the DefaultDict if you care about the order of the table.

            >>> d=OrderedDict()
            >>> d['name']=['bright','dim']
            >>> d['flux']=['large','small']
            >>> print latex_table(d, units=dict(flux='[ergs]'))
            \begin{deluxetable}{cc}
            \tablehead{\colhead{name} & \colhead{flux}\\ \colhead{ } & \colhead{[ergs]}}
            \startdata
            bright & large \\
            dim & small \\
            \enddata
            \end{deluxetable}


    """

    if latexdict is None: 
        latexdict=dict()

    if units is not None:
        latexdict['units']=units

    outtable=StringIO()
        
    asciitable.write(table_dict, outtable, 
                     Writer=asciitable.AASTex,
                     names=table_dict.keys(),
                     latexdict=latexdict,
                     **kwargs
                    )
    t=outtable.getvalue()
    return t.strip()

class TableFormatter(object):
    def __init__(self, confluence=False, precision=2):
        """ Default is to format numbers for latex tables.
            confluence=True will format numbers of confluence tables. 
            
            A few examples:

                >>> conf_format=TableFormatter(confluence=True)
                >>> latex_format=TableFormatter(confluence=False)

                >>> print conf_format.value(-1, precision=2)
                \-1.00
                >>> print latex_format.value(-1, precision=2)
                -1.00

                >>> print conf_format.error(10,2, precision=0)
                10 +/- 2
                >>> print latex_format.error(10,2, precision=0)
                $10 \pm 2$

                >>> print conf_format.ul(-1, precision=1)
                <\-1.0
                >>> print latex_format.ul(-1, precision=1)
                $<-1.0$
        """
        self.confluence = confluence
        self.precision = precision
    @staticmethod
    def fix_negative(a):
        return a.replace('-',r'\-')
    def value(self,a,precision=None):
        if a is None: return "None"
        if precision is None: precision=self.precision
        ret = '%.*f' % (precision,a)
        if self.confluence:
            ret = TableFormatter.fix_negative(ret)
        return ret
    def error(self,a,b,precision=None):
        if precision is None: precision=self.precision
        if self.confluence:
            f1=TableFormatter.fix_negative('%.*f' % (precision,a))
            f2=TableFormatter.fix_negative('%.*f' % (precision,b))
            return '%s +/- %s' % (f1,f2)
        else:
            return '$%.*f \pm %.*f$' % (precision,a,precision,b)
    def ul(self,ul,precision=None):
        if precision is None: precision=self.precision
        ret=r'<%.*f' % (precision,ul)
        if not self.confluence:
            ret = '$%s$' % ret
        if self.confluence:
            ret=TableFormatter.fix_negative(ret)
        return ret

    @property
    def nodata(self):
        return '' if self.confluence else r'\nodata'
    
if __name__ == "__main__":
    import doctest
    doctest.testmod()
