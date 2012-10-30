#!/usr/bin/env python
import numpy as np

from uw.like.pointspec import SpectralAnalysis,DataSpecification
from uw.like.pointspec_helpers import get_default_diffuse, PointSource
from uw.like.roi_catalogs import Catalog2FGL
from uw.like.Models import PowerLaw

class RadioPSRROIBuilder(object):
    def __init__(self, radiopsr_loader):
        self.radiopsr_loader = radiopsr_loader

    def build_roi(self, name):

        catalog = Catalog2FGL('$FERMI/catalogs/gll_psc_v05.fit', 
                              latextdir='$FERMI/extended_archives/gll_psc_v05_templates',
                              prune_radius=0,
                              max_free=5,
                              free_radius=5,
                              limit_parameters=True)

        ft1 = self.radiopsr_loader.get_ft1(name)
        ft2 = self.radiopsr_loader.get_ft2(name)
        ltcube = self.radiopsr_loader.get_ltcube(name)
        binfile = self.radiopsr_loader.get_binfile(name)

        roi_dir = self.radiopsr_loader.get_skydir(name)
        roi_size = 10

        ds = DataSpecification(
            ft1files = ft1,
            ft2files = ft2,
            ltcube   = ltcube,
            binfile  = binfile)

        sa = SpectralAnalysis(ds,
                              binsperdec = 4,
                              emin       = 100,
                              emax       = 1000000,
                              irf        = "P7SOURCE_V6",
                              roi_dir    = roi_dir,
                              maxROI     = roi_size,
                              minROI     = roi_size,
                              event_class= 0)

        fit_emin = 1e3
        fit_emax = 10**5.5

        model=PowerLaw(index=2, e0=np.sqrt(fit_emin*fit_emax))
        model.set_limits('index',-5,5)

        ps = PointSource(name=name, model=model, skydir=roi_dir)
        point_sources = [ps]

        diffuse_sources = get_default_diffuse(diffdir="/afs/slac/g/glast/groups/diffuse/rings/2year",
                                   gfile="ring_2year_P76_v0.fits",
                                   ifile="isotrop_2year_P76_source_v0.txt",
                                   limit_parameters=True)

        roi=sa.roi(point_sources=point_sources,
                   diffuse_sources=diffuse_sources,
                   catalogs=catalog,
                   fit_emin=fit_emin, fit_emax=fit_emax)
        return roi
