from os.path import expandvars

import numpy as np

from SED import SED

from uw.like.Models import PowerLaw

import lande_units as units

class TeVSED(object):
    """ General SED-like object
        for loading in publically
        avaliable HESS-like SED data.

        This data can typically be found at:
            http://www.mpi-hd.mpg.de/hfm/HESS/pages/publications.
    """
    def __init__(self, sed_file, energy_units_str='MeV', flux_units_str='erg'):
        self.sed_file = sed_file

        self.energy_units_str = energy_units_str
        self.flux_units_str   = flux_units_str
        self.energy_units     = units.fromstring(self.energy_units_str)
        self.flux_units       = units.fromstring(self.flux_units_str)

        self.__load__(sed_file)

    def __load__(self, sed_file):
        """ Parse a HESS-style plain text SED file.
            Note, there is no real consistency to
            these files, so just try many
            examples. """

        lines = open(expandvars(sed_file)).readlines()

        # strip whitespace + newline characters
        lines = [line.strip() for line in lines]

        # remove blank lines + comments
        lines = [line for line in lines if line is not '' and line[0] != '#' ]

        # easier to parse '/TeV/cm^2/s' than '/TeV cm^2 s'
        lines = [line.replace('[/TeV cm^2 s]','[/TeV/cm^2/s]') for line in lines]

        # split on empty spaces
        lines = [ line.split() for line in lines ]

        energy_units='[TeV]'
        flux_units='[/TeV/cm^2/s]'

        header = lines[0]
        units_line = lines[1]

        for u in units_line:
            if u not in [energy_units, flux_units]:
                raise Exception("%s is an unrecognized unit" % u)

        # get all the data lines
        data = lines[2:]
        
        # convert to float + transpose rows + convert to numpy
        data = [ map(float,i) for i in data]
        data = zip(*data)
        data = [np.asarray(i,dtype=float) for i in data]

        if header == ['Energy','Flux','Flux','Error_low','Flux','Error_high']:
            # This is like SEDs from the galactic plane survey
            #  * http://www.mpi-hd.mpg.de/hfm/HESS/pages/publications/auxiliary/ApJ_636.html

            self.energy, self.dnde, dnde_lower, dnde_upper = data
            self.dnde_err = (dnde_lower + dnde_upper)/2

            # lines are energy, flux, flux_lower, flux_upper

        elif header == ['Energy','interval','Mean','energy','Flux','Flux','Error']:
            # This is like the SED for HESS J1303-631
            #  * http://www.mpi-hd.mpg.de/hfm/HESS/pages/publications/auxiliary/AA439_1013.html
            print 'data',data

            self.lower_energy, self.upper_energy, self.energy, self.dnde, self.dnde_err = data

        elif header == ['Mean','energy','Flux','Flux','Error']:
            # This is like the SED for Vela X
            #  * http://www.mpi-hd.mpg.de/hfm/HESS/pages/publications/auxiliary/VelaX_auxinfo.html

            self.energy, self.dnde, self.dnde_err = data

        elif header == ['Energy','interval','Mean','energy','Flux','Flux','Error']:
            # This is like the SED for MSH 15-52
            #  * http://www.mpi-hd.mpg.de/hfm/HESS/pages/publications/auxiliary/AA435_L17.html
            self.lower_energy, self.upper_energy, self.dnde, self.dnde_err = data

        elif header == ['Mean','energy','Flux','Flux','Error','(','-','/','+',')']:
            # This is like the SED for the Kookaburra
            #  * http://www.mpi-hd.mpg.de/hfm/HESS/pages/publications/auxiliary/kookaburra_auxinfo.html
            raise Exception("not implemented")

        if self.energy is None:
            self.energy = np.sqrt(self.lower_energy*self.upper_energy)
        if self.lower_energy is None or self.upper_energy is None:
            # if only energy is given, assume no range in energy
            self.lower_energy = self.energy
            self.upper_energy = self.energy

        self.significant = np.ones_like(self.energy).astype(bool)

        def estimate_flux(dnde,energy,emin,emax, e_weight):
            """ estimate the emin to emax flux for a source with prefactor
                dnde at the given energy. assuming the source has a
                spectral index of 2. 
                
                Note, for our situation dnde has units [ph/cm^2/s/TeV]
                and energy has units [TeV], but the ouptut of the
                i_flux function is correct. If e_weight=0,
                the return has units [ph/cm^2/s]. If e_weight=1,
                the return has units [TeV/cm^2/s]. """
            model = PowerLaw(index=2, norm=dnde, e0=energy)
            return model.i_flux(emin=emin, emax=emax, e_weight=e_weight)

        estimate_flux = np.vectorize(estimate_flux)

        self.flux = estimate_flux(self.dnde, self.energy, self.lower_energy, self.upper_energy, 0)
        self.eflux = estimate_flux(self.dnde, self.energy, self.lower_energy, self.upper_energy, 1)

        # Note, we can use this shortcut to comptue the errors,
        # But ONLY because the function flux(dnde) is linear in dnde and dnde is the only
        # parameter with an error: flux=dnde*f(gamma,e0) so flux_err=dnde_err*f(gamma,e0)
        self.flux_err = estimate_flux(self.dnde_err, self.energy, self.lower_energy, self.upper_energy, 0)
        self.eflux_err = estimate_flux(self.dnde_err, self.energy, self.lower_energy, self.upper_energy, 1)


        # no upper limits
        self.dnde_ul=np.nan*self.dnde
        self.flux_ul=np.nan*self.flux
        self.eflux_ul=np.nan*self.eflux

        for values, u in [
            [['lower_energy', 'upper_energy', 'energy'], units.TeV],
            [['dnde', 'dnde_err', 'dnde_ul'], units.ph/units.cm**2/units.s/units.TeV],
            [['flux', 'flux_err', 'flux_ul'], units.ph/units.cm**2/units.s],
            [['eflux','eflux_err', 'eflux_ul'], units.TeV/units.cm**2/units.s]]:

            for v in values:
                self.__dict__[v] = units.tosympy(self.__dict__[v], u)

        # all hess points are significant (as far as I can tell)
        self.significant = np.ones(len(self.energy),dtype=bool)
        print 'en',self.energy
        print 'sig',self.significant 

    def plot(self, filename=None,
             axes=None, 
             fignum=None, figsize=(4,4),
             data_kwargs=dict()):
        """ Plot the SED using matpotlib. """

        if axes is None:
            fig = P.figure(fignum,figsize)
            axes = fig.add_axes((0.22,0.15,0.75,0.8))
            self.set_xlim(axes,
                          float(self.lower_energy[0]/self.energy_units),
                          float(self.upper_energy[-1]/self.energy_units))
        self.axes = axes

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

        if filename is not None: P.savefig(filename)
        return axes
