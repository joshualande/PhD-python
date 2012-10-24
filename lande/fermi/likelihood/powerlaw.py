import numpy as np
import pylab as P
from scipy.optimize import fmin

from uw.utilities import keyword_options
from uw.like.Models import PowerLaw

from lande.pysed import units

from . basefit import BaseFitter
from . save import spectrum_to_dict
from . specplot import SpectralAxes, SpectrumPlotter

class PowerLawApproximator(BaseFitter):

    defaults = BaseFitter.defaults + (
        ('npoints',1000,'number of points in fit'),
        ('e0',None,'scale for power law'),
        ('energy_units', 'MeV', 'default units to plot energy flux (y axis) in.'),
        ('flux_units',  'erg', 'default units to plot energy (x axis) in'),
)

    @keyword_options.decorate(defaults)
    def __init__(self, input_model, emin, emax, **kwargs):
        """ Create an approximate power law spectrum. """

        raise Exception("This code doesn't work yet. I think you need the exposure to do the fit correctly.")
        self.input_model = input_model
        self.emin = emin
        self.emax = emax

        keyword_options.process(self, kwargs)

        self._calculate()

    def _calculate(self):
        self.results = dict()

        energies = np.logspace(np.log10(self.emin),np.log10(self.emax),self.npoints)

        if self.e0 is None: 
            self.e0=np.sqrt(self.emin*self.emax)

        self.results['input_model'] = spectrum_to_dict(self.input_model)
        self.results['dnde'] = dnde = self.input_model(energies)

        self.pl_model = PowerLaw(e0=self.e0)
        self.pl_model.set_flux(self.input_model.i_flux(emin=self.emin,emax=self.emax),
                    emin=self.emin,emax=self.emax)

        def residuals(args):
            norm,index=args
            self.pl_model['norm']=norm
            self.pl_model['index']=index
            dnde_pl = self.pl_model(energies)
            #return np.sum((np.log(dnde) - np.log(dnde_pl))**2)
            print (np.log10(dnde)-np.log10(dnde_pl))**2
            return np.sum((np.log(dnde) - np.log(dnde_pl))**2)
            #return np.sum((dnde - dnde_pl)**2)

        best_norm,best_index=fmin(residuals,[self.pl_model['norm'],self.pl_model['index']])
        self.pl_model['norm']=best_norm
        self.pl_model['index']=best_index

        self.results['pl_model'] = spectrum_to_dict(self.pl_model)


    def plot(self,filename=None,axes=None,fignum=None,figsize=(4,4)):
        if axes is None:
            fig = P.figure(fignum,figsize)
            axes = SpectralAxes(fig=fig, 
                                rect=(0.22,0.15,0.75,0.8),
                                flux_units=self.flux_units,
                                energy_units=self.energy_units)
            fig.add_axes(axes)
            axes.set_xlim_units(self.emin*units.MeV, self.emax*units.MeV)

        sp=SpectrumPlotter(axes=axes)
        sp.plot(self.results['input_model'], label='input')
        sp.plot(self.results['pl_model'], label='powerlaw')
        
        if filename is not None:
            P.savefig(filename)

    if __name__ == "__main__":
        import doctest
        doctest.testmod()
