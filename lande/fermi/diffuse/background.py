import numpy as np
import re
from os.path import basename

from skymaps import IsotropicConstant

from uw.like.Models import PowerLaw
from uw.like.roi_diffuse import DiffuseSource
from uw.like.pointspec_helpers import get_diffuse_source

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

