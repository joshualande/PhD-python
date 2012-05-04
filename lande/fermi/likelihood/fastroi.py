""" This simple object just creates an ROI very quickly. This is useful for debugging,
    but of no other paritcular use. 
    
    Author: Joshua Lande
"""
# this import has to come first (shrugs)
from . roi_gtlike import Gtlike

from os.path import join

from skymaps import SkyDir
from uw.like.pointspec import DataSpecification
from uw.like.pointspec_helpers import PointSource
from uw.like.Models import PowerLaw
from uw.like.roi_monte_carlo import SpectralAnalysisMC
from uw.utilities import keyword_options

from . diffuse import get_sreekumar

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
        ('irf','P7SOURCE_V6'),
        ('flux',1e-7),
        ('isotropic_bg',False, 'simulate on top of an isotropic background'),
        ('conv_type',0)
    )

    @keyword_options.decorate(defaults)
    def __init__(self,**kwargs):
        keyword_options.process(self, kwargs)

        point_sources, diffuse_sources = self.get_sources()


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
                                conv_type=self.conv_type,
                                roi_dir=self.roi_dir,
                                minROI=self.maxROI,
                                maxROI=self.maxROI,
                                irf=self.irf,
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

    def get_sources(self):

        point_sources, diffuse_sources = [], []

        model = PowerLaw()
        model.set_flux(self.flux, self.emin, self.emax)
        ps = PointSource(
            name = 'source',
            model = model,
            skydir = self.roi_dir)
        point_sources.append(ps)

        if self.isotropic_bg:
            ds = get_sreekumar()
            diffuse_sources.append(ds)

        if diffuse_sources == []: diffuse_sources = None

        return point_sources, diffuse_sources

    def get_roi(self):
        return self.roi

    def get_like(self, *args, **kwargs):

        self.gtlike=Gtlike(self.roi, savedir=self.savedir, 
                      fix_pointlike_ltcube=True,
                      *args, **kwargs)
        return self.gtlike.like
