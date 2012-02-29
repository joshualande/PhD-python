import numpy as np
import pyfits as pf

def expand_fits(pf,factor,hdu=0):
    """ Create a new fits file where
        each pixel is divided into factor**2
        more pixels all of equal value. 
        
        I am not an expert on Fits files and
        I only think this code works. Also,
        I have no idea if it would work with 3D data (?).
        As a rule of thumb, if you use this function,
        you should probably open the before/after
        in ds9 and blink them for a while to convince
        yourself the function worked.
    """

    h=pf[hdu].header
    d=pf[hdu].data

    h['CDELT1']/=factor
    h['CDELT2']/=factor

    h['NAXIS1']*=factor
    h['NAXIS2']*=factor

    h['CRPIX1']=h['CRPIX1']*factor - factor/2.0 + 0.5
    h['CRPIX2']=h['CRPIX2']*factor - factor/2.0 + 0.5

    larger=list(d.shape)
    larger[0]*=factor
    larger[1]*=factor
    larger_array = np.zeros(larger,dtype=d.dtype)
    for i in range(factor):
        for j in range(factor):
            larger_array[i::factor,j::factor] = d

    pf[hdu].data=larger_array
