import copy
from collections import OrderedDict

from StringIO import StringIO
import asciitable


def confluence_table(dict, units=None, **kwargs):
    """ Format a dictionary (or a DefaultDict) into a confluence table. 
    
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
        dict = copy.deepcopy(dict)
        for k in dict:
            dict[k].insert(0,units.pop(k,''))

    outtable=StringIO()

    asciitable.write(dict, outtable,
                 Writer=asciitable.FixedWidth,
                 names=dict.keys(),
                 **kwargs)
    t=outtable.getvalue()
    t=t.replace(' |','  |')
    t=t.strip().split('\n')
    t[0]=t[0].replace(' |','||')
    t=['|'+i for i in t]
    return '\n'.join(t)


def latex_table(dict, units=None, latexdict=None, **kwargs):

    if latexdict is None: 
        latexdict=dict()

    if units is not None:
        latexdict['units']=units

    outtable=StringIO.StringIO()
        
    asciitable.write(table, outtable, 
                     Writer=asciitable.AASTex,
                     names=table.keys(),
                     **kwargs
                    )
    t=outtable.getvalue()
    return t
    
if __name__ == "__main__":
    import doctest
    doctest.testmod()
