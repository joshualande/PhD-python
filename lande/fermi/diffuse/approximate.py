import numpy as np

from skymaps import IsotropicConstant

from uw.like.Models import FileFunction
from uw.like.roi_diffuse import DiffuseSource

class ApproximateIsotropic(DiffuseSource):

    def __init__(self,
                 name,
                 diffuse_sources,
                 file,
                 emin, emax,
                 skydir,
                 scaling_factor=1):
        """ Approximates a linear combination of
            diffuse sources as a single isotropic
            spectrum. 
            
            Note, probably does not work for extended
            sources. """

        self.name = name
        self.scaling_factor = scaling_factor

        for ds in diffuse_sources:
            if len(ds.dmodel) != 1: raise Exception("dmodels must have length 1")

        self.file  = open(file,'w')
        self.file.write(
            self._make_file(diffuse_sources, skydir, emin, emax, self.scaling_factor)
            )
        self.file.close()

        self.smodel = FileFunction(file=file)
        self.dmodel = [IsotropicConstant()]

    @staticmethod
    def get_dnde(diffuse_sources, skydir, energy, scaling_factor):
        return scaling_factor*sum(ds.smodel(energy)*ds.dmodel[0](skydir, energy) for ds in diffuse_sources)

    @staticmethod
    def _make_file(diffuse_sources, skydir, emin, emax, scaling_factor, npts=1e3):

        energies = np.logspace(np.log10(emin), np.log10(emax), npts)

        fluxes = [ ApproximateIsotropic.get_dnde(diffuse_sources, skydir, energy, scaling_factor) for energy in energies ]

        return '\n'.join(['%s\t%s' % (e,f) for e,f in zip(energies, fluxes)])
