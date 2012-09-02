#!/usr/bin/env python
import yaml
import numbers
import os
from glob import glob
from os.path import expandvars, join, exists
from tempfile import mkdtemp
import numbers
import shutil

import numpy as np

from skymaps import SkyDir

from uw.like.pointspec import SpectralAnalysis,DataSpecification
from uw.like.pointspec_helpers import get_default_diffuse, PointSource
from uw.like.SpatialModels import Gaussian
from uw.like.roi_catalogs import Catalog2FGL
from uw.like.roi_extended import ExtendedSource
from uw.like.roi_save import load
from uw.like.Models import PowerLaw
from uw.utilities import phasetools
from uw.pulsar.phase_range import PhaseRange
from uw.utilities.parmap import LogMapper,LimitMapper

isnum = lambda x: isinstance(x, numbers.Real)


import cPickle
def load_pwn(filename, savedir=None, **kwargs):
    """ By default, pwn roi's are not loadable. """
    d=cPickle.load(open(expandvars(filename),'r'))

    extra=d['extra']

    ft1=extra['unphased_ft1']
    ltcube=extra['unphased_ltcube']
    phase=PhaseRange(extra['phase'])

    if exists(ft1) and exists(ltcube) and exists(phase):
        return load(filename, **kwargs)

    if savedir is None:
        save_data = False
        savedir=mkdtemp(prefix='/scratch/')
    else:
        if not os.path.exists(savedir):
            os.makedirs(savedir)
        save_data = True

    print 'new savedir is',savedir

    phased_ft1=PWNRegion.phase_ft1(ft1,phase,savedir)
    phased_ltcube=PWNRegion.phase_ltcube(ltcube,phase,savedir)
    binfile=join(savedir,'binned_phased.fits')

    roi = load(filename, ft1files=phased_ft1, ltcube=phased_ltcube, binfile=binfile, **kwargs)

    if not save_data:
        roi.__del__ = lambda x: shutil.rmtree(savedir)

    return roi



class PWNRegion(object):

    def __init__(self, pwndata, savedir=None):

        self.pwndata = pwndata

        if savedir is None: 
            self.savedata=False
            self.savedir=mkdtemp(prefix='/scratch/')
        else:
            self.savedir=savedir
            self.savedata=True
            if os.path.exists(self.savedir):
                for file in glob(join(self.savedir,'*')):
                    os.remove(file)
            else:
                os.makedirs(self.savedir)

    @staticmethod
    def get_background():
        return get_default_diffuse(diffdir="/afs/slac/g/glast/groups/diffuse/rings/2year",
                                   gfile="ring_2year_P76_v0.fits",
                                   ifile="isotrop_2year_P76_source_v0.txt",
                                   limit_parameters=True)

    @staticmethod
    def get_source(name, position, 
                   fit_emin, fit_emax, 
                   extended=False, sigma=None):
        """ build a souce. """
        model=PowerLaw(index=2, e0=np.sqrt(fit_emin*fit_emax))
        # don't limit prefactor. Use oomp limits in gtlike.
        model.set_limits('index',-5,5)
        flux=PowerLaw(norm=1e-11, index=2, e0=1e3).i_flux(fit_emin,fit_emax)
        model.set_flux(flux,emin=fit_emin,emax=fit_emax)

        if extended and sigma != 0:
            if not isnum(sigma): raise Exception("sigma must be set. """)
            return ExtendedSource(
                name=name,
                model=model,
                spatial_model=Gaussian(sigma=sigma, center=position))
        else:
            return PointSource(
                name=name,
                model=model,
                skydir=position)

    @staticmethod
    def get_catalog(**kwargs):
        return Catalog2FGL('$FERMI/catalogs/gll_psc_v05.fit', 
                           latextdir='$FERMI/extended_archives/gll_psc_v05_templates',
                           prune_radius=0,
                           limit_parameters=True,
                           **kwargs)

    @staticmethod
    def phase_ltcube(ltcube,phase,savedir):

        if np.allclose(phase.phase_fraction,1):
            phased_ltcube = ltcube
        else:
            # create a temporary ltcube scaled by the phase factor
            phased_ltcube=join(savedir,'phased_ltcube.fits')
            if not os.path.exists(phased_ltcube):
                phasetools.phase_ltcube(ltcube,phased_ltcube, phase=phase)
        return phased_ltcube

    @staticmethod
    def phase_ft1(ft1,phase,savedir):

        if np.allclose(phase.phase_fraction,1):
            phased_ft1 = ft1
        else:
            # apply phase cut to ft1 file
            phased_ft1 = join(savedir,'ft1_phased.fits')
            if not os.path.exists(phased_ft1):
                phasetools.phase_cut(ft1,phased_ft1,phaseranges=phase.tolist(dense=False))
        return phased_ft1

    def get_roi(self, name, phase, 
                fit_emin, fit_emax, 
                binsperdec,
                extended=False, 
                roi_size=10, 
                catalog_kwargs=dict(),
                **kwargs):
        """ Sets up the ROI for studying a LAT Pulsar in the off pulse. """

        sourcedict=yaml.load(open(self.pwndata))[name]
        ltcube=sourcedict['ltcube']
        pulsar_position=SkyDir(*sourcedict['cel'])
        ft1=sourcedict['ft1']
        ft2=sourcedict['ft2']

        source = PWNRegion.get_source(name, 
                                      position = pulsar_position,
                                      fit_emin = fit_emin, 
                                      fit_emax = fit_emax, 
                                      sigma = 0.1,
                                      extended=extended)
        sources = [source]
        roi_dir = pulsar_position

        phase = PhaseRange(phase)

        point_sources, diffuse_sources = [],[]
        for source in sources:
            if isinstance(source,PointSource):
                point_sources.append(source)
            else:
                diffuse_sources.append(source)

        diffuse_sources += PWNRegion.get_background()

        catalog=PWNRegion.get_catalog(**catalog_kwargs)

        binfile=join(self.savedir,'binned_phased.fits')

        phased_ltcube=PWNRegion.phase_ltcube(ltcube,phase,self.savedir)
        phased_ft1=PWNRegion.phase_ft1(ft1,phase,self.savedir)

        ds = DataSpecification(
            ft1files = phased_ft1,
            ft2files = ft2,
            ltcube   = phased_ltcube,
            binfile  = binfile)

        print 'For now, 4 bins per decade. Eventually, this will have to be better.'
        sa = SpectralAnalysis(ds,
                              binsperdec = binsperdec,
                              emin       = 100,
                              emax       = 1000000,
                              irf        = "P7SOURCE_V6",
                              roi_dir    = roi_dir,
                              maxROI     = roi_size,
                              minROI     = roi_size,
                              event_class= 0)

        roi=sa.roi(point_sources=point_sources,
                   diffuse_sources=diffuse_sources,
                   catalogs=catalog,
                   phase_factor=1,
                   fit_emin=fit_emin, fit_emax=fit_emax,
                   **kwargs)

        print 'bins ',roi.bin_edges

        roi.extra = dict(
            unphased_ft1 = ft1,
            unphased_ltcube = ltcube,
            phase = phase)

        self.roi = roi

        return roi

    def __del__(self):
        if not self.savedata:
            print 'Removing temporary directory',self.savedir
            shutil.rmtree(self.savedir)
