""" This simple object just creates an ROI very quickly. This is useful for debugging,
    but of no other paritcular use. 
    
    Author: Joshua Lande
"""
from os.path import join

from skymaps import SkyDir
from uw.like.pointspec import DataSpecification
from uw.like.pointspec_helpers import PointSource
from uw.like.Models import PowerLaw
from uw.like.roi_monte_carlo import SpectralAnalysisMC
from uw.utilities import keyword_options

from . roi_gtlike import Gtlike


class FastROI(object):
    """ Usage:

            fast = FastROI()

            # pointlike
            roi = fast.get_roi()

            # gtlike
            like = fast.get_like()

    """
    defaults = (
        ('emin', 1e4),
        ('emax', 1e5),
        ('binsperdec', 1),
        ('roi_dir', SkyDir()),
        ('savedir', 'savedir'),
        ('seed', 0),
        ('maxROI', 5),
        ('simtime', 2629743.83, 'simulation time (in seconds)'),
    )



    @keyword_options.decorate(defaults)
    def __init__(self,**kwargs):
        keyword_options.process(self, kwargs)

        model = PowerLaw()
        model.set_flux(1e-8, self.emin, self.emax)

        ps = PointSource(
            name = 'source',
            model = model,
            skydir = self.roi_dir)

        point_sources = [ ps ] 
        diffuse_sources = None

        ds = DataSpecification(
            ft1files = join(self.savedir,'ft1.fits'),
            ft2files = join(self.savedir,'ft2.fits'),
            ltcube = join(self.savedir,'ltcube.fits'),
            binfile = join(self.savedir,'binfile.fits')
        )

        sa = SpectralAnalysisMC(ds,
                                seed=self.seed,
                                emin=self.emin,
                                emax=self.emax,
                                binsperdec=self.binsperdec,
                                event_class=0,
                                roi_dir=self.roi_dir,
                                minROI=self.maxROI,
                                maxROI=self.maxROI,
                                irf='P7SOURCE_V6',
                                use_weighted_livetime=True,
                                savedir=self.savedir,
                                tstart=0,
                                tstop=self.simtime,
                                ltfrac=0.9,
                               )

        roi = sa.roi(roi_dir=self.roi_dir,
                     point_sources = point_sources,
                     diffuse_sources = diffuse_sources)

        self.roi = roi

    def get_roi(self):
        return self.roi

    def get_like(self, *args, **kwargs):

        return Gtlike(self, savedir=self.savedir, *args, **kwargs)
