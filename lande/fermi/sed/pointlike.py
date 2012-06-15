from collections import defaultdict

import yaml
import numpy as np

from lande.utilities.tools import tolist


def pointlike_sed_to_dict(bandflux):
    results = defaultdict(lambda:defaultdict(list))

    results['Energy']['Units'] = 'MeV'
    results['dNdE']['Units'] = 'ph/cm^2/s/MeV'

    results['Significant'] = []
    for r in bandflux.rec:
        
        results['Energy']['Lower'].append(r.elow)
        results['Energy']['Upper'].append(r.ehigh)
        results['Energy']['Value'].append(np.sqrt(r.elow*r.ehigh))

        # Undo scaling in the bandflux recarray
        fac = r.elow*r.ehigh*bandflux.scale_factor

        if r.flux > 0:
            results['Significant'].append(True)
            results['dNdE']['Value'].append(r.flux/fac)
            results['dNdE']['Average_Error'].append((r.uflux/fac - r.lflux/fac)/2)
            results['dNdE']['Lower_Error'].append((r.flux-r.lflux)/fac)
            results['dNdE']['Upper_Error'].append((r.uflux-r.flux)/fac)
            results['dNdE']['Upper_Limit'].append(np.nan)
        else:
            results['Significant'].append(False)
            results['dNdE']['Value'].append(np.nan)
            results['dNdE']['Error'].append(np.nan)
            results['dNdE']['Upper_Limit'].append(r.uflux/fac)

    return tolist(results)

def pointlike_sed_to_yaml(bandflux, filename):
    """ Save a pointlike SED to a yaml file.

        Usage: 
            
            bf = roi.plot_sed(which=which, filename='plot.png')
            pointlike_sed_to_yaml(bf, filename='data.yaml')
    """
    d=pointlike_sed_to_dict(bandflux)
    open(filename,'w').write(yaml.dump(d))
