import pyfits
import os
from subprocess import check_call
from skymaps import DiffuseFunction

from uw.like.Models import Constant
from uw.like.roi_diffuse import DiffuseSource
from uw.like.roi_monte_carlo import MCModelBuilder

def merge_diffuse(diffuse_sources, scaling_model=None, mergefile=None, verbosity=False, short_name=False , compress=False):
    """ merge diffuse files into one file. """

    for i in diffuse_sources:
        assert MCModelBuilder.isone(i.smodel)
        assert isinstance(i.dmodel[0], DiffuseFunction)

    if not short_name:
        merged_name = 'merged+%s' % ('+'.join([i.name for i in diffuse_sources]))
    else:
        merged_name = 'merged'

    if scaling_model is None:
        scaling_model = Constant()

    filenames = [ds.dmodel[0].name() for ds in diffuse_sources]
    if verbosity:
        print 'Merging files:'
        for file in filenames:
            print ' .. ',file

    if mergefile is None:
        if not compress:
            mergefile = 'merged.fits'
        else:
            mergefile = 'merged.fits.gz'

    
    if compress:
        if mergefile.find(".gz")>0:
            temp_file=mergefile.remove(".gz")
        else:
            temp_file=mergefile
            mergefile=mergefile+".gz"
    

    pf = [pyfits.open(i) for i in filenames]

    first = pf.pop(0)
    for i in pf:
        assert first['PRIMARY'].data.shape == i['PRIMARY'].data.shape
        if first['PRIMARY'].header != i['PRIMARY'].header:
            print 'WARNING: HEADER FILES DISAGREE. PERFORMING MERGE ANYWAY'
        first['PRIMARY'].data += i['PRIMARY'].data

    if not compress:
        first.writeto(mergefile, clobber=True)
    else:
        first.writeto(temp_file, clobber=True)
        
        if os.path.exists(mergefile):
                os.remove(mergefile)
        check_call(['gzip', temp_file])
        
    merged_dmodel=DiffuseFunction(mergefile)

    new_source = DiffuseSource(
        name=merged_name,
        scaling_model=scaling_model,
        diffuse_model=merged_dmodel)
    return new_source

