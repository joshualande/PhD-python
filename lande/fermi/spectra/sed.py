""" Implements a subclass of SED which has nicer features. 

    Author: J. Lande
"""
from collections import defaultdict
from os.path import expandvars

import yaml
import pylab as P
import numpy as np

import pyLikelihood

# overload SED namespace for our base object
from SED import SED as BaseGtlikeSED
_funcFactory = pyLikelihood.SourceFactory_funcFactory()

from uw.like.Models import Model

from lande.pysed import units
from lande.utilities.tools import tolist
from uw.utilities import keyword_options


class SED(object):
    """ Base object for plotting SEDs.

        The input must be XXX
        Additional input

        Something about optional significant flag, what to do with asymetrical errors.

    """
    defaults= (
        ('energy_units', 'MeV', 'default units to plot energy flux (y axis) in.'),
        ('flux_units',  'erg', 'default units to plot energy (x axis) in'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, results, *args, **kwargs):
        keyword_options.process(self, kwargs)

        self.energy_units_str = self.energy_units
        self.flux_units_str   = self.flux_units

        self.energy_units     = units.fromstring(self.energy_units_str)
        self.flux_units       = units.fromstring(self.flux_units_str)

        self.fromdict(results)


    @staticmethod
    def dict_to_spectrum(d):
        """ Load back as a pyLikelihood spectrum object
            a spectrum that has been saved by the spectrum_to_string
            object. This undoes the conversion of BaseGtlikeSED.spectrum_to_dict """
        spectrum=_funcFactory.create(d['name'])
        for k,v in d.items(): 
            if k != 'name' and k[-4:] != '_err': spectrum.getParam(k).setTrueValue(v)
        return spectrum


    def fromdict(self,results):
        """ Update the internal values in this object
            from a dictionary of SED points. """

        if isinstance(results, str): results = yaml.load(open(expandvars(results)))

        edict = results['Energy']
        fdict = results['dNdE']

        e = lambda x: units.tosympy(x, units.fromstring(edict['Units']))
        dnde = lambda x: units.tosympy(x, units.fromstring(fdict['Units']))

        if results.has_key('Name'): self.name = results['Name']

        self.energy = e(edict['Value'])

        if 'Lower' in edict and 'Upper' in edict:
            self.lower_energy = e(edict['Lower'])
            self.upper_energy = e(edict['Upper'])
            self.has_assymetric_errors = True
        else:
            self.has_assymetric_errors = False

        self.dnde = dnde(fdict['Value'])
        self.dnde_err = dnde(fdict['Average_Error'])

        if 'Lower_Error' in fdict and 'Upper_Error' in fdict:
            # assymetric errors
            self.dnde_lower_err = dnde(fdict['Lower_Error'])
            self.dnde_upper_err = dnde(fdict['Upper_Error'])
            self.has_energy_errors = True
        else:
            self.has_energy_errors = False

        # get limits, otherwise assume all significant
        if 'Upper_Limit' in fdict and 'Significant' in results:
            self.dnde_ul = dnde(fdict['Upper_Limit'])
            self.significant = np.asarray(results['Significant'])

            self.has_upper_limits=True
        else:
            self.has_upper_limits=False

        if 'spectrum' in results:
            self.spectrum = results['spectrum']
            self.has_spectrum=True
        else:
            self.has_spectrum=False


    @staticmethod
    def get_dnde(spectrum,energies):
        if isinstance(spectrum,pyLikelihood.Function):
            return SED.get_dnde(spectrum,energies)
        elif isinstance(spectrum,Model):
            return spectrum(energies)
        else:
            raise Exception("Unrecognized type %s for spectrum." % type(spectrum))

    @staticmethod
    def _plot_spectrum(spectrum, axes, energy_units, flux_units, npts=100, **kwargs):
        """ Plot a pyLikelihood spectrum onto a matplotlib axes
            whose x axis has units energy_units and y_axes has units
            flux_units/cm^2/s. 
            
            N.B. This function corrects for the fact that gtlike always
            internally represnts the spectrum in units of ph/cm^2/s/MeV. """
        low_lim, hi_lim = axes.get_xlim()
        energies = np.logspace(np.log10(low_lim), np.log10(hi_lim), npts)

        e_units = units.tosympy(energies,energy_units)

        energies_mev = units.tonumpy(e_units, units.MeV)

        # (a) convert to acutal units. gtlike spectra take energy in MeV, return flux in ph/cm^2/s/MeV
        dnde = units.tosympy(SED.get_dnde(spectrum,energies_mev),units.ph/units.cm**2/units.s/units.MeV)
        # (b) create E^2 dN/dE in acutal units
        e2_dnde = dnde.multiply_elementwise(e_units).multiply_elementwise(e_units)
        # (c) convert to desired units
        e2_dnde = units.tonumpy(e2_dnde,flux_units/units.cm**2/units.s)
        axes.plot(energies, e2_dnde, **kwargs)

    def plot_spectrum(self, spectrum, axes, **kwargs):
        if axes is None: 
            axes=self.axes
        SED._plot_spectrum(spectrum, axes, self.energy_units, self.flux_units, **kwargs)

    def plot(self, filename=None,
             axes=None, 
             fignum=None, figsize=(4,4),
             plot_spectral_fit=True,
             data_kwargs=dict(),
             spectral_kwargs=dict(color='red',zorder=1.9)):
        """ Plot the SED using matpotlib. """

        if axes is None:
            fig = P.figure(fignum,figsize)
            axes = fig.add_axes((0.22,0.15,0.75,0.8))
            
            if self.has_assymetric_errors:
                lower = float(self.lower_energy[0]/self.energy_units)
                upper = float(self.upper_energy[-1]/self.energy_units)
                axes.set_xlim(lower,upper)

        self.axes = axes

        ce = lambda x: units.tonumpy(x, self.energy_units)
        cf = lambda y: units.tonumpy(
            y.multiply_elementwise(self.energy).multiply_elementwise(self.energy),
            self.flux_units/units.cm**2/units.s)

        if plot_spectral_fit and self.has_spectrum:
            # kind of ugly, but spectrum is ph/cm^2/s/MeV
            # and it gets mutlplied by energy_units**2,
            # so we need to multiple overall spectrum by
            # flux_units*MeV/energy_units**2
            self.plot_spectrum(self.spectrum, axes=axes, **spectral_kwargs)

        BaseGtlikeSED._plot_points(
            x=ce(self.energy),
            xlo=ce(self.lower_energy) if self.has_energy_errors else None,
            xhi=ce(self.upper_energy) if self.has_energy_errors else None,
            y=cf(self.dnde),
            y_lower_err=cf(self.dnde_lower_err) if self.has_assymetric_errors else cf(self.dnde_err),
            y_upper_err=cf(self.dnde_upper_err) if self.has_assymetric_errors else cf(self.dnde_err),
            y_ul=cf(self.dnde_ul) if self.has_upper_limits else None,
            significant=self.significant if self.has_upper_limits else np.ones(len(self.energy),dtype=bool),
            energy_units=self.energy_units_str,
            flux_units=self.flux_units_str,
            axes=axes, **data_kwargs)

        if filename is not None: P.savefig(filename)
        return axes



class GtlikeSED(SED,BaseGtlikeSED):
    """ object to make SEDs using pyLikelihood. 
    
        Currently, this object only allows the SED
        points to be the same as the binning in
        the FT1 file. """

    def __init__(self, like, *args, **kwargs):
        """ Additional Parameters
            * flux_units - 
            * nergy_units - . """

        self.energy_units_str = kwargs.pop('energy_units','MeV')
        self.flux_units_str   = kwargs.pop('flux_units','erg')
        self.energy_units     = units.fromstring(self.energy_units_str)
        self.flux_units       = units.fromstring(self.flux_units_str)

        BaseGtlikeSED.__init__(self, like, *args, **kwargs)

    def _calculate(self,*args,**kwargs):
        """ Convert all units into sympy arrays after the initial calculation. """

        # easier to represnt upper limits as NaN if you use yaml to load/dump
        self.dnde_ul=np.nan*self.dnde_ul
        self.flux_ul=np.nan*self.flux_ul
        self.eflux_ul=np.nan*self.eflux_ul

        try:
            super(GtlikeSED,self)._calculate(*args,**kwargs)
            self.crashed = False
        except Exception, ex:
            print 'ERROR computing SED:', ex
            for v in ['dnde', 'dnde_err', 'dnde_lower_err', 'dnde_upper_err', 'dnde_ul',
                      'flux', 'flux_err', 'flux_ul',
                      'eflux','eflux_err', 'eflux_ul']:
                self.__dict__[v] *= np.nan

            self.crashed = True
            self.significant = np.zeros_like(self.energy).astype(bool)


        for values, u in [
            [['lower_energy', 'upper_energy', 'energy'], units.MeV],
            [['dnde', 'dnde_err', 'dnde_lower_err', 'dnde_upper_err', 'dnde_ul'], units.ph/units.cm**2/units.s/units.MeV],
            [['flux', 'flux_err', 'flux_ul'], units.ph/units.cm**2/units.s],
            [['eflux','eflux_err', 'eflux_ul'], units.MeV/units.cm**2/units.s]]:

            for v in values:
                self.__dict__[v] = units.tosympy(self.__dict__[v], u)

    def __str__(self):
        results = self.todict()
        return yaml.dump(results)

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
                Average_Error=c_dnde(self.dnde_err),
                Lower_Error=c_dnde(self.dnde_lower_err),
                Upper_Error=c_dnde(self.dnde_upper_err),
                Upper_Limit=c_dnde(self.dnde_ul),
                Units='ph/cm^2/s/%s' % self.flux_units_str),
            Ph_Flux=dict(
                Value=c_flux(self.flux),
                Average_Error=c_flux(self.flux_err),
                Upper_Limit=c_flux(self.flux_ul),
                Units='ph/cm^2/s'),
            En_Flux=dict(
                Value=c_eflux(self.eflux),
                Average_Error=c_eflux(self.eflux_err),
                Upper_Limit=c_eflux(self.eflux_ul),
                Units='%s/cm^2/s' % self.flux_units_str),
            Test_Statistic=self.ts.tolist(),
            Significant=self.significant.tolist(),
            Spectrum=BaseGtlikeSED.spectrum_to_dict(self.spectrum))

