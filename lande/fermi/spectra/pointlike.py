from collections import defaultdict

import yaml
import numpy as np

from uw.utilities import keyword_options
from uw.like.sed_plotter import BandFlux

from lande.pysed import units
from lande.utilities.tools import tolist

from lande.fermi.likelihood.save import spectrum_to_dict

from . sed import SED

class PointlikeSED(SED):

    defaults = SED.defaults + (
        ('merge', True, 'merge edge bins'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, roi, name, **kwargs):
        keyword_options.process(self, kwargs)
        self.roi = roi
        self.name = name
        
        bf = BandFlux(self.roi, which=self.name, merge=self.merge, scale_factor=1)
        results = PointlikeSED.pointlike_sed_to_dict(bf, flux_units=self.flux_units, energy_units=self.energy_units)

        model = roi.get_model(name)
        results['spectrum'] = spectrum_to_dict(model)

        super(PointlikeSED,self).__init__(results, **keyword_options.defaults_to_kwargs(self, SED))
        

    @staticmethod
    def pointlike_sed_to_dict(bandflux, flux_units='erg', energy_units='MeV'):
        results = defaultdict(lambda:defaultdict(list))
        
        ce=lambda e: units.convert(e,'MeV',energy_units)
        cp = lambda e: units.convert(e,'1/MeV','1/%s' % flux_units)

        results['Energy']['Units'] = energy_units
        results['dNdE']['Units'] = 'ph/cm^2/s/%s' % flux_units

        results['Significant'] = []
        for r in bandflux.rec:
            
            results['Energy']['Lower'].append(ce(r.elow))
            results['Energy']['Upper'].append(ce(r.ehigh))
            results['Energy']['Value'].append(ce(np.sqrt(r.elow*r.ehigh)))

            # Undo scaling in the bandflux recarray
            fac = r.elow*r.ehigh*bandflux.scale_factor

            if r.flux > 0:
                results['Significant'].append(True)
                results['dNdE']['Value'].append(cp(r.flux/fac))
                results['dNdE']['Average_Error'].append(cp((r.uflux/fac - r.lflux/fac)/2))
                results['dNdE']['Lower_Error'].append(cp((r.flux-r.lflux)/fac))
                results['dNdE']['Upper_Error'].append(cp((r.uflux-r.flux)/fac))
                results['dNdE']['Upper_Limit'].append(np.nan)
            else:
                results['Significant'].append(False)
                results['dNdE']['Value'].append(np.nan)
                results['dNdE']['Average_Error'].append(np.nan)
                results['dNdE']['Lower_Error'].append(np.nan)
                results['dNdE']['Upper_Error'].append(np.nan)
                results['dNdE']['Upper_Limit'].append(cp(r.uflux/fac))

        return tolist(results)


def pointlike_sed_to_yaml(bandflux, filename):
    """ helper function to save a pointlike SED (generated using the roi.plot_sed function)
        to a yaml file. Generally, it is perferable to use PointlikeSED
        to generate SEDs.

        Usage: 
            
            bf = roi.plot_sed(which=which, filename='plot.png')
            pointlike_sed_to_yaml(bf, filename='data.yaml')
    """
    d=PointlikeSED.pointlike_sed_to_dict(bandflux)
    open(filename,'w').write(yaml.dump(d))
