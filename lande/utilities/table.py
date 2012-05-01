import copy
from collections import OrderedDict

from StringIO import StringIO
import asciitable


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
        ||   name ||   flux ||
        ||         | [ergs]  |
        || bright  |  large  |
        ||    dim  |  small  |

    """

    if units is not None:
        table_dict = copy.deepcopy(table_dict)
        for k in table_dict:
            table_dict[k].insert(0,units.pop(k,''))

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
    
if __name__ == "__main__":
    import doctest
    doctest.testmod()
