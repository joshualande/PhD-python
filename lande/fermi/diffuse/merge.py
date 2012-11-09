import pyfits

from skymaps import DiffuseFunction

from uw.like.Models import Constant
from uw.like.roi_diffuse import DiffuseSource
from uw.like.roi_monte_carlo import MCModelBuilder

def merge_diffuse(diffuse_sources, scaling_model=None, mergefile=None, verbosity=False):
    """ merge diffuse files into one file. """

    for i in diffuse_sources:
        assert MCModelBuilder.isone(i.smodel)
        assert isinstance(i.dmodel[0], DiffuseFunction)


    merged_name = 'merged(%s)' % (','.join([i.name for i in diffuse_sources]))

    if scaling_model is None:
        scaling_model = Constant()

    filenames = [ds.dmodel[0].name() for ds in diffuse_sources]
    if verbosity:
        print 'Merging files:'
        for file in filenames:
            print ' .. ',file

    if mergefile is None:
        mergefile = 'merged.fits'

    pf = [pyfits.open(i) for i in filenames]

    first = pf.pop(0)
    for i in pf:
        assert first['PRIMARY'].header == i['PRIMARY'].header
        first['PRIMARY'].data += i['PRIMARY'].data
    first.writeto(mergefile, clobber=True)

    merged_dmodel=DiffuseFunction(mergefile)

    new_source = DiffuseSource(
        name=merged_name,
        scaling_model=scaling_model,
        diffuse_model=merged_dmodel)
    return new_source

