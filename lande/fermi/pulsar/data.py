from os.path import expandvars

import pylab as P
import numpy as np

from uw.utilities.fitstools import rad_extract
from uw.utilities.decorators import memoize
from uw.pulsar.phase_range import PhaseRange

@memoize
def get_all_phases(ft1, skydir):
    """ Cache photons = faster """
    ed = rad_extract(expandvars(ft1),skydir,radius_function=180,return_cols=['PULSE_PHASE', 'TIME'])
    return ed


def get_phases_and_times(ft1, skydir, emin, emax, radius):
    ed = get_all_phases(ft1, skydir)
    all_phases = ed['PULSE_PHASE']
    all_times = ed['TIME']

    cut = (ed['ENERGY'] >= emin) & (ed['ENERGY'] <= emax) & (ed['DIFFERENCES'] <= np.radians(radius))
    return all_phases[cut], all_times[cut]

def get_phases(*args, **kwargs):
    return get_phases_and_times(*args, **kwargs)[0]



