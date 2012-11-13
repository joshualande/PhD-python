import copy

import numpy as np
import pyfits as pf

def expand_fits_header(header, factor):
    header = header.copy()
    header['CDELT1']/=factor
    header['CDELT2']/=factor

    header['NAXIS1']*=factor
    header['NAXIS2']*=factor

    header['CRPIX1']=header['CRPIX1']*factor - factor/2.0 + 0.5
    header['CRPIX2']=header['CRPIX2']*factor - factor/2.0 + 0.5
    return header

def expand_fits_data(data, factor):
    larger=list(data.shape)
    larger[0]*=factor
    larger[1]*=factor
    larger_data = np.zeros(larger,dtype=data.dtype)
    for i in range(factor):
        for j in range(factor):
            larger_data[i::factor,j::factor] = data

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
    pf = copy.deepcopy(pf)

    pf[hdu].header = expand_fits_header(pf[hdu].header)
    pf[hdu].data=expand_fits_data(pf[hdu].data)
