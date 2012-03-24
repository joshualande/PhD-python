from GtApp import GtApp
import pyfits

def fix_pointlike_ltcube(ltcube):

    f = pyfits.open(ltcube, mode='update')

    for exposure in ['EXPOSURE', 'WEIGHTED_EXPOSURE']:
        if exposure in [i.name for i in f]:
            h=f[exposure].header

            if 'NDSKEYS' not in h.keys():
                h.update('DSTYP1', 'TIME')
                h.update('DSUNI1', 's')
                h.update('DSVAL1', 'TABLE')
                h.update('DSREF1', ':GTI')
                h.update('NDSKEYS', 1)
    f.flush()

def gtltcube(**kwargs):
    l=GtApp('gtltcube', 'Likelihood')
    l.run(**kwargs)

