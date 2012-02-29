import re
from tempfile import NamedTemporaryFile

import pyfits

from GtApp import GtApp

def make_ebinfile(ebins, filename):
    """ Make an ebinfile suitable for gtbin given an array of energy
        bins. Not sure why, but this is the required format. """

    c1=pyfits.Column(name='E_MIN', format='D', array=ebins[:-1])
    c2=pyfits.Column(name='E_MAX', format='D', array=ebins[1:])

    t=pyfits.new_table([c1, c2])
    t.update_ext_name('ENERGYBINS')
    h=pyfits.HDUList([pyfits.PrimaryHDU(), t])

    h.writeto(filename)


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


    if naxis == 2:
        kwargs['algorithm'] = 'CMAP'

    elif naxis == 3:

        ebins = b['ENERGIES'].data.field('ENERGY')

        temp=NamedTemporaryFile(delete=True)
        ebinfile = temp.name
        make_ebinfile(ebins, ebinfile)

        kwargs.update(
            algorithm = 'CCUBE',
            ebinalg = 'FILE',
            ebinfile = ebinfile,
        )

    GtApp('gtbin').run(
        evfile=evfile, 
        outfile=outfile, 
        scfile=scfile, 
        **kwargs)
