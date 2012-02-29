""" Implements a subclass of SED which has nicer features. 

    Author: J. Lande
"""
from collections import defaultdict
from os.path import expandvars

import yaml
import pylab as P
import numpy as np

import pyLikelihood
from SED import SED
_funcFactory = pyLikelihood.SourceFactory_funcFactory()

import lande_units as units
from lande_toolbag import tolist



class SuperSED(SED):
    """ object to make SEDs using pyLikelihood. 
    
        Currently, this object only allows the SED
        points to be the same as the binning in
        the FT1 file. """

    def __init__(self, like, *args, **kwargs):
        """ Additional Parameters
            * flux_units - desired units to quote energy flux (y axis) in.
            * nergy_units - desired units to quote energy (x axis) in. """

        self.energy_units_str = kwargs.pop('energy_units','MeV')
        self.flux_units_str   = kwargs.pop('flux_units','erg')
        self.energy_units     = units.fromstring(self.energy_units_str)
        self.flux_units       = units.fromstring(self.flux_units_str)

        if isinstance(like,dict)or isinstance(like,str):
            self.fromdict(like, *args, **kwargs)
        else:
            super(LandeSED,self).__init__(like, *args, **kwargs)

    @staticmethod
    def dict_to_spectrum(d):
        """ Load back as a pyLikelihood spectrum object
            a spectrum that has been saved by the spectrum_to_string
            object. This undoes the conversion of SED.spectrum_to_dict """
        spectrum=_funcFactory.create(d['name'])
        for k,v in d.items(): 
            if k != 'name' and k[-4:] != '_err': spectrum.getParam(k).setTrueValue(v)
        return spectrum


    def _calculate(self,*args,**kwargs):
        """ Convert all units into sympy arrays after the initial calculation. """

        # easier to represnt upper limits as NaN if you use yaml to load/dump
        self.dnde_ul=np.nan*self.dnde_ul
        self.flux_ul=np.nan*self.flux_ul
        self.eflux_ul=np.nan*self.eflux_ul

        try:
            super(LandeSED,self)._calculate(*args,**kwargs)
            self.crashed = False
        except Exception, ex:
            print 'ERROR computing SED:', ex
            for v in ['dnde', 'dnde_err', 'dnde_ul',
                      'flux', 'flux_err', 'flux_ul',
                      'eflux','eflux_err', 'eflux_ul']:
                self.__dict__[v] *= np.nan

            self.crashed = True
            self.significant = np.zeros_like(self.energy).astype(bool)


        for values, u in [
            [['lower_energy', 'upper_energy', 'energy'], units.MeV],
            [['dnde', 'dnde_err', 'dnde_ul'], units.ph/units.cm**2/units.s/units.MeV],
            [['flux', 'flux_err', 'flux_ul'], units.ph/units.cm**2/units.s],
            [['eflux','eflux_err', 'eflux_ul'], units.MeV/units.cm**2/units.s]]:

            for v in values:
                self.__dict__[v] = units.tosympy(self.__dict__[v], u)

    def todict(self):
        """ Return a dictionary of the SED points with the desired units. """

        if self.crashed: return dict()

        c_energy=lambda x: units.tonumpy(x,self.energy_units).tolist()
        c_dnde=lambda x: units.tonumpy(x,units.ph/units.cm**2/units.s/self.flux_units).tolist()
        c_flux=lambda x: units.tonumpy(x,units.ph/units.cm**2/units.s).tolist()
        c_eflux=lambda x: units.tonumpy(x,self.flux_units/units.cm**2/units.s).tolist()

        return dict(
            Name=self.name,
            Energy=dict(
                Value=c_energy(self.energy),
                Lower=c_energy(self.lower_energy),
                Upper=c_energy(self.upper_energy),
                Units='%s' % self.energy_units_str),
            dNdE=dict(
                Value=c_dnde(self.dnde),
                Error=c_dnde(self.dnde_err),
                Upper_Limit=c_dnde(self.dnde_ul),
                Units='ph/cm^2/s/%s' % self.flux_units_str),
            Ph_Flux=dict(
                Value=c_flux(self.flux),
                Error=c_flux(self.flux_err),
                Upper_Limit=c_flux(self.flux_ul),
                Units='ph/cm^2/s'),
            En_Flux=dict(
                Value=c_eflux(self.eflux),
                Error=c_eflux(self.eflux_err),
                Upper_Limit=c_eflux(self.eflux_ul),
                Units='%s/cm^2/s' % self.flux_units_str),
            Test_Statistic=self.ts.tolist(),
            Significant=self.significant.tolist(),
            Spectrum=SED.spectrum_to_dict(self.spectrum))

    def fromdict(self,d):
        """ Update the internal values in this object
            from a dictionary of SED points. """

        if isinstance(d, str): d = yaml.load(open(expandvars(d)))

        e = lambda x: units.tosympy(x, units.fromstring(d['Energy']['Units']))
        dnde = lambda x: units.tosympy(x, units.fromstring(d['dNdE']['Units']))
        flux = lambda x: units.tosympy(x, units.fromstring(d['Ph_Flux']['Units']))
        eflux = lambda x: units.tosympy(x, units.fromstring(d['En_Flux']['Units']))

        if d.has_key('name'): self.name = d['Name']

        self.lower_energy = e(d['Energy']['Lower'])
        self.upper_energy = e(d['Energy']['Upper'])
        self.energy = e(d['Energy']['Value'])

        self.dnde = dnde(d['dNdE']['Value'])
        self.dnde_err = dnde(d['dNdE']['Error'])
        self.dnde_ul = dnde(d['dNdE']['Upper_Limit'])

        if d.has_key('Ph_flux'):
            self.flux = flux(d['Ph_Flux']['Value'])
            self.flux_err = flux(d['Ph_Flux']['Error'])
            self.flux_ul = flux(d['Ph_Flux']['Upper_Limit'])

        if d.has_key('En_flux'):
            self.eflux = eflux(d['En_Flux']['Value'])
            self.eflux_err = eflux(d['En_Flux']['Error'])
            self.eflux_ul = eflux(d['En_Flux']['Upper_Limit'])

        if d.has_key('Test_Statistic'):
            self.ts = np.asarray(d['Test_Statistic'])
        self.significant = np.asarray(d['Significant'])

        if d.has_key('Spectrum'):
            self.spectrum = LandeSED.dict_to_spectrum(d['Spectrum'])

        self.crashed = False

    def __str__(self):
        results = self.todict()
        return yaml.dump(results)

    @staticmethod
    def _plot_spectrum(spectrum, axes, energy_units, flux_units, npts=100, **kwargs):
        """ Plot a pyLikelihood spectrum onto a matplotlib axes
            whose x axis has units energy_units and y_axes has units
            flux_units/cm^2/s. 
            
            N.B. This function corrects for the fact that gtlike always
            internally represnts the spectrum in units of ph/cm^2/s/MeV. """
        low_lim, hi_lim = axes.get_xlim()
        energies = np.logspace(np.log10(low_lim), np.log10(hi_lim), npts)

        # (a) convert to acutal units
        dnde = units.tosympy(SED.get_dnde(spectrum,energies),units.ph/units.cm**2/units.s/units.MeV)
        e_units = units.tosympy(energies,energy_units)
        # (b) create E^2 dN/dE in acutal units
        e2_dnde = dnde.multiply_elementwise(e_units).multiply_elementwise(e_units)
        # (c) convert to desired units
        e2_dnde = units.tonumpy(e2_dnde,flux_units/units.cm**2/units.s)
        axes.plot(energies, e2_dnde, **kwargs)

    def plot_spectrum(self, spectrum, **kwargs):
        LandeSED._plot_spectrum(spectrum, self.axes, self.energy_units, self.flux_units, **kwargs)

    def plot(self, filename=None,
             axes=None, 
             fignum=None, figsize=(4,4),
             plot_spectral_fit=True,
             data_kwargs=dict(),
             spectral_kwargs=dict(color='red')):
        """ Plot the SED using matpotlib. """

        if axes is None:
            fig = P.figure(fignum,figsize)
            axes = fig.add_axes((0.22,0.15,0.75,0.8))
            self.set_xlim(axes,
                          float(self.lower_energy[0]/self.energy_units),
                          float(self.upper_energy[-1]/self.energy_units))
        self.axes = axes

        if not self.crashed:

            ce = lambda x: units.tonumpy(x, self.energy_units)
            cf = lambda y: units.tonumpy(
                y.multiply_elementwise(self.energy).multiply_elementwise(self.energy),
                self.flux_units/units.cm**2/units.s)

            SED._plot_points(
                x=ce(self.energy), 
                xlo=ce(self.lower_energy), 
                xhi=ce(self.upper_energy), 
                y=cf(self.dnde),
                y_err=cf(self.dnde_err),
                y_ul=cf(self.dnde_ul),
                significant=self.significant,
                energy_units=self.energy_units_str,
                flux_units=self.flux_units_str,
                axes=axes, **data_kwargs)

            if plot_spectral_fit and hasattr(self,'spectrum'):
                # kind of ugly, but spectrum is ph/cm^2/s/MeV
                # and it gets mutlplied by energy_units**2,
                # so we need to multiple overall spectrum by
                # flux_units*MeV/energy_units**2
                self.plot_spectrum(self.spectrum, **spectral_kwargs)

        if filename is not None: P.savefig(filename)
        return axes
