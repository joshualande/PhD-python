""" This simple object just creates an ROI very quickly. This is useful for debugging,
    but of no other paritcular use. 
    
    Author: Joshua Lande
"""
# this import has to come first (shrugs)
from lande.fermi.likelihood.roi_gtlike import Gtlike

import numpy as np
from os.path import join

from skymaps import SkyDir
from uw.like.pointspec import DataSpecification
from uw.like.pointspec_helpers import PointSource
from uw.like.Models import PowerLaw
from uw.like.roi_monte_carlo import SpectralAnalysisMC
from uw.utilities import keyword_options

from lande.fermi.diffuse.background import get_sreekumar
from lande.fermi.data.livetime import fix_pointlike_ltcube

class FastROI(object):
    """ Usage:

            from lande.fermi.testing.fastroi import FastROI
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
        ('tempdir', 'tempdir'),
        ('seed', 0),
        ('maxROI', 5),
        ('simtime', 2629743.83, 'simulation time (in seconds)'),
        ('irf','P7SOURCE_V6'),
        ('flux',1e-7),
        ('isotropic_bg',False, 'simulate on top of an isotropic background'),
        ('nearby_source',False, ''),
        ('conv_type',0),
        ('point_sources',[]),
        ('diffuse_sources',[]),
        ('powerlaw_index',2),
    )

    @keyword_options.decorate(defaults)
    def __init__(self,**kwargs):
        keyword_options.process(self, kwargs)

        if self.point_sources == [] and self.diffuse_sources == []:
            self.point_sources, self.diffuse_sources = self.get_default_sources()

        ltcube = join(self.tempdir,'ltcube.fits')
        ds = DataSpecification(
            ft1files = join(self.tempdir,'ft1.fits'),
            ft2files = join(self.tempdir,'ft2.fits'),
            ltcube = ltcube, 
            binfile = join(self.tempdir,'binfile.fits')
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
                                savedir=self.tempdir,
                                tstart=0,
                                tstop=self.simtime,
                                ltfrac=0.9,
                               )

        roi = sa.roi(roi_dir=self.roi_dir,
                     point_sources = self.point_sources,
                     diffuse_sources = self.diffuse_sources)

        self.roi = roi

        fix_pointlike_ltcube(ltcube)

    def get_default_sources(self):

        point_sources, diffuse_sources = [], []

        model = PowerLaw(index=self.powerlaw_index, e0=np.sqrt(self.emin*self.emax))
        model.set_flux(self.flux, emin=self.emin, emax=self.emax)
        ps = PointSource(
            name = 'source',
            model = model.copy(),
            skydir = self.roi_dir)
        point_sources.append(ps)

        if self.isotropic_bg:
            ds = get_sreekumar()
            diffuse_sources.append(ds)

        if self.nearby_source:
            ps = PointSource(
                name = 'nearby_source',
                model = model.copy(),
                skydir = SkyDir(self.roi_dir.ra(),self.roi_dir.dec()+3)
            )
            point_sources.append(ps)

        return point_sources, diffuse_sources

    def get_roi(self):
        return self.roi

    def get_like(self, *args, **kwargs):

        self.gtlike=Gtlike(self.roi, savedir=self.tempdir, 
                      fix_pointlike_ltcube=True,
                      *args, **kwargs)
        return self.gtlike.like
