import numpy as np
import re
from os.path import basename

from skymaps import IsotropicConstant

from uw.like.Models import PowerLaw, FileFunction
from uw.like.roi_diffuse import DiffuseSource
from uw.like.pointspec_helpers import get_diffuse_source

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

        
def get_background(*args):
    bg = []
    for source in args:
        if re.search(r'\.fit(s)?(\.gz)?$',source) is not None:
            bg.append(get_diffuse_source('MapCubeFunction',source,'PowerLaw',None,name=basename(source)))
        elif re.search(r'\.txt$',source) is not None:
            bg.append(get_diffuse_source('ConstantValue',None,'FileFunction',source,name=basename(source)))
        else:
            raise Exception("Diffuse Sources must end in .fit, .fits, .fit.gz, .fits.gz, or .txt (file is %s)" % basename(source))

    return bg[0] if len(args)==1 else bg

def get_sreekumar(diff_factor=1, free=(True, False)):

    # use Sreekumar-like defaults
    if diff_factor == 1:
        name = 'Sreekumar Isotropic'
    else:
        name = 'Sreekumar Isotropic x%s' % diff_factor

    free = np.asarray(free).copy()
    model = PowerLaw(index=2.1, free=free)
    model.set_flux(1.5e-5*diff_factor, emin=100, emax=np.inf)

    return DiffuseSource(
        name=name,
        diffuse_model=IsotropicConstant(),
        scaling_model=model)

