import re
from tempfile import NamedTemporaryFile

import numpy as np
import pyfits

from GtApp import GtApp


def gtbin_from_binfile(evfile, outfile, scfile, binfile):
    """ Run gtbin on ft1 file evfile to create outfile outfile.
        but make the outfile have a binning consistent with the reference
        binfile. """

    global b
    b = pyfits.open(binfile)

    h = b['primary'].header

    naxis = h['naxis']

    # This should be true of all binned maps
    assert naxis in [2,3]
    assert abs(h['cdelt1']) == abs(h['cdelt2'])
    binsz=abs(h['cdelt1'])

    # split on 1 or more dashes
    m=re.compile('-+')
    xname,projx=m.split(h['ctype1'])
    yname,projy=m.split(h['ctype2'])

    assert projx == projy
    proj = projx

    if xname == 'GLON' and yname == 'GLAT':
        coordsys = 'GAL'
    elif xname== 'RA' and yname == 'DEC':
        coordsys = 'CEL'
    else:
        raise Exception("Unrecognized coordinate system (%s,%s)" % (xname,yname))

    kwargs = dict(
        nxpix=h['naxis1'],
        nypix=h['naxis2'],
        binsz=binsz,
        xref=h['crval1'],
        yref=h['crval2'],
        axisrot=h['crota2'],  # I hope this is right - J.L.
        coordsys=coordsys,
        proj=proj,
    )


    if naxis == 2 or (naxis == 3 and h['naxis3'] == 1):
        # 2D image or 3D image with 3rd axis only one deep
        kwargs['algorithm'] = 'CMAP'

    elif naxis == 3:

        ebins = b['ENERGIES'].data.field('ENERGY')

        emin = ebins[0]
        emax = ebins[-1]
        enumbins = len(ebins)

        # For now, just crash if ebinalg != 'LOG'
        assert np.allclose(ebins, np.logspace(np.log10(emin), np.log10(emax), enumbins))

        kwargs.update(
            algorithm = 'CCUBE',
            ebinalg = 'LOG',
            emin = emin,
            emax = emax,
            enumbins = enumbins,
        )

    GtApp('gtbin').run(
        evfile=evfile, 
        outfile=outfile, 
        scfile=scfile, 
        **kwargs)
