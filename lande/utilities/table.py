import copy
import sys
from argparse import ArgumentParser
from collections import OrderedDict

from StringIO import StringIO
import asciitable

from . tools import parse_strip_known_args

def fixed_width_table(table_dict, table_kwargs=dict(bookend=False, delimiter=None), indent=None):
    outtable=StringIO()
    asciitable.write(table_dict, outtable,
                     Writer=asciitable.FixedWidth,
                     names=table_dict.keys(),
                     **table_kwargs)
    t=outtable.getvalue()

    # remove final newline
    assert t[-1] == '\n'
    t=t[:-1]
    if indent is None:
        return t

    t=t.split('\n')
    return '\n'.join(['%s%s' % (indent,i) for i in t])

def get_table_type():
    parser = ArgumentParser()
    parser.add_argument("--confluence", default='latex', action="store_const", const='confluence')
    args=parse_strip_known_args(parser)
    return args.confluence

def t2t_table(table_dict, units=None, **kwargs):
    """ Create a t2t format table. 

        this is required by t2t for tables
        see for example: http://txt2tags.org/markup.html

        Example:
            >>> d=OrderedDict()
            >>> d['name']=['bright','dim']
            >>> d['flux']=['large','small']
            >>> print t2t_table(d)
            ||  name |  flux |
            | bright | large |
            |    dim | small |

            >>> print t2t_table(d, units=dict(flux='[ergs]'))
            ||  name | flux [ergs] |
            | bright |       large |
            |    dim |       small |
     """
    if units is not None:
        for k,v in units.items():
            table_dict[k + ' %s' % v] = table_dict.pop(k)

    outtable=StringIO()

    asciitable.write(table_dict, outtable,
                     Writer=asciitable.FixedWidth,
                     names=table_dict.keys(),
                     **kwargs
                    )
    t=outtable.getvalue()

    t='||' + t[2:]
    return t.strip()


def confluence_table(table_dict, units=None, **kwargs):
    """ Format a dictionary into a confluence table. 
        Use the DefaultDict if you care about the order of the table.
    
        Example:
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
    def __init__(self, table_type, precision=2):
        """ Default is to format numbers for latex tables.
            confluence=True will format numbers of confluence tables. 
            
            A few examples:

                >>> conf_format=TableFormatter(table_type='confluence')
                >>> latex_format=TableFormatter(table_type='latex')

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
        assert table_type in ['confluence', 'latex']

        self.table_type = table_type
        self.precision = precision
    @staticmethod
    def fix_negative(a):
        return a.replace('-',r'\-')
    def value(self,a,precision=None):
        if a is None: return "None"
        if precision is None: precision=self.precision
        ret = '%.*f' % (precision,a)
        if self.table_type == 'confluence':
            return TableFormatter.fix_negative(ret)
        elif self.table_type == 'latex':
            return ret
    def error(self,a,b,precision=None):
        if precision is None: precision=self.precision
        if self.table_type == 'confluence':
            f1=TableFormatter.fix_negative('%.*f' % (precision,a))
            f2=TableFormatter.fix_negative('%.*f' % (precision,b))
            return '%s +/- %s' % (f1,f2)
        elif self.table_type == 'latex':
            return '$%.*f \pm %.*f$' % (precision,a,precision,b)
    def ul(self,ul,precision=None):
        if precision is None: precision=self.precision
        ret=r'<%.*f' % (precision,ul)
        if self.table_type == 'latex':
            ret = '$%s$' % ret
        elif self.table_type == 'confluence':
            ret=TableFormatter.fix_negative(ret)
        return ret

    @property
    def nodata(self):
        if self.table_type == 'confluence':
            return ''
        elif self.table_type == 'latex':
            return ''
    
if __name__ == "__main__":
    import doctest
    doctest.testmod()
