#!/usr/bin/env python
# This import has to come first
from lande.fermi.likelihood.roi_gtlike import Gtlike

import os
from os.path import join
from collections import defaultdict
from argparse import ArgumentParser

import yaml
import numpy as np

from skymaps import SkyDir

from uw.like.SpatialModels import Gaussian
from uw.like.SpatialModels import Gaussian


from lande.utilities.save import loaddict,savedict
from lande.utilities.load import import_module
from lande.utilities.argumentparsing import argparse_to_kwargs

from lande.fermi.pulsar.plotting import plot_phaseogram, plot_phase_vs_time
from lande.fermi.likelihood.tools import force_gradient
from lande.fermi.likelihood.parlimits import all_params_limited
from lande.fermi.likelihood.variability import GtlikeVariabilityTester
from lande.fermi.likelihood.save import pointlike_dict_to_spectrum
from lande.fermi.likelihood.free import freeze_far_away, unfreeze_far_away


from . setup import PWNRegion, load_pwn
from . pointlike import pointlike_analysis
from . gtlike import gtlike_analysis
from . plots import plots, tsmaps

class Pipeline(object):

    @staticmethod
    def get_kwargs():

        parser = ArgumentParser()
        parser.add_argument("--pwndata", required=True)
        parser.add_argument("--pwnphase")
        parser.add_argument("--name", required=True, help="Name of the pulsar")
        parser.add_argument("--emin", default=1e2, type=float)
        parser.add_argument("--emax", default=10**5.5, type=float)
        parser.add_argument("--binsperdec", default=4, type=int)
        parser.add_argument("--use-gradient", default=False, action="store_true")
        parser.add_argument("--no-at-pulsar", default=False, action="store_true")
        parser.add_argument("--no-point", default=False, action="store_true")
        parser.add_argument("--no-extended", default=False, action="store_true")
        parser.add_argument("--no-cutoff", default=False, action="store_true")
        parser.add_argument("--no-savedir", default=False, action="store_true")
        parser.add_argument("--max-free", default=5, type=float)
        parser.add_argument("--modify", required=True)
        args=parser.parse_args()
        return argparse_to_kwargs(args)


    def __init__(self, **kwargs):
        self.__dict__.update(**kwargs)

        force_gradient(use_gradient=self.use_gradient)
        np.seterr(all='ignore')

        self.datadir='data' 
        self.plotdir='plots'
        self.seddir='seds'

        for dir in [self.seddir, self.datadir, self.plotdir]: 
            if not os.path.exists(dir): os.makedirs(dir)

    def main(self):
        do_at_pulsar = not self.no_at_pulsar
        do_point = not self.no_point
        do_extended = not self.no_extended 
        do_cutoff = not self.no_cutoff


        name=self.name
        emin=self.emin
        emax=self.emax

        pwnphase=yaml.load(open(self.pwnphase))[name]
        phase=pwnphase['phase']

        print 'phase = ',phase

        pwndata=yaml.load(open(self.pwndata))[name]

        savedir=None if self.no_savedir else join(os.getenv('PWD'),'savedir')

        results=r=defaultdict(lambda: defaultdict(dict))
        results['name']=name
        results['phase']=phase

        pointlike_kwargs=dict(name=name, max_free=self.max_free) 

        print 'Building the ROI'
        reg=PWNRegion(pwndata=self.pwndata, savedir=savedir)
        roi=reg.get_roi(name=name, phase=phase, 
                        catalog_kwargs=dict(free_radius=5, max_free=self.max_free),
                        fit_emin=emin, fit_emax=emax, 
                        binsperdec=self.binsperdec,
                        extended=False)

        modify = import_module(self.modify)
        roi.extra['phase'] = phase
        roi.extra['new_sources'] = modify.modify_roi(name,roi)
        roi.extra['pwnphase'] = pwnphase

        assert all_params_limited(roi, except_sources=[name])

        model1=modify.cutoff_model1(name)

        if do_at_pulsar:
            r['at_pulsar']['pointlike']=pointlike_analysis(roi, hypothesis='at_pulsar', 
                                                           seddir=self.seddir, datadir=self.datadir, 
                                                           cutoff=do_cutoff, **pointlike_kwargs)
        if do_point:
            r['point']['pointlike']=pointlike_analysis(roi, hypothesis='point', localize=True, 
                                                       seddir=self.seddir, datadir=self.datadir, 
                                                       cutoff=do_cutoff, 
                                                       model1 = model1,
                                                       **pointlike_kwargs)
        if do_extended:
            roi.modify(which=name, spatial_model=Gaussian(sigma=0.1), keep_old_center=True)

            r['extended']['pointlike']=pointlike_analysis(roi, hypothesis='extended', cutoff=False, 
                                                          seddir=self.seddir, datadir=self.datadir, 
                                                          fit_extension=True, 
                                                          **pointlike_kwargs)

        savedict(results,'results_%s_pointlike.yaml' % name)

    def reload_roi(self,hypothesis):
        name = self.name
        roi = load_pwn('roi_%s_%s.dat' % (hypothesis,name))
        return roi

    def gtlike_followup(self, hypothesis):

        name = self.name
        roi = self.reload_roi(hypothesis)

        cutoff = (not self.no_cutoff) and hypothesis in ['at_pulsar', 'point']
        upper_limit = hypothesis=='at_pulsar'
        if cutoff:
            pointlike_results = loaddict('results_%s_pointlike.yaml' % name)
            model1=pointlike_results[hypothesis]['pointlike']['test_cutoff']['model_1']
            model1=pointlike_dict_to_spectrum(model1)
            model1.set_default_limits(oomp_limits=True)
        else:
            model1=None

        results = {hypothesis:{}}
        results[hypothesis]['gtlike']=gtlike_analysis(roi, name=name,
                                                      max_free = self.max_free,
                                                      seddir=self.seddir, datadir=self.datadir, plotdir=self.plotdir,
                                                      hypothesis=hypothesis, 
                                                      upper_limit=upper_limit,
                                                      cutoff=cutoff,
                                                      model1=model1,
                                                     )

        savedict(results,'results_%s_gtlike_%s.yaml' % (name,hypothesis))

    def tsmaps_followup(self, hypothesis):
        name = self.name
        roi = self.reload_roi(hypothesis)

        if not os.path.exists('plots'): os.makedirs('plots')

        pwndata=yaml.load(open(self.pwndata))[name]
        pulsar_position = SkyDir(*pwndata['cel'])
        
        new_sources = roi.extra['new_sources']
        tsmaps(roi, name, hypothesis, new_sources=new_sources, pulsar_position=pulsar_position,
               datadir=self.datadir, plotdir=self.plotdir)

    def plots_followup(self, hypothesis):
        name = self.name
        roi = self.reload_roi(hypothesis)

        if not os.path.exists('plots'): os.makedirs('plots')

        pwndata=yaml.load(open(self.pwndata))[name]
        ft1 = pwndata['ft1']
        pulsar_position = SkyDir(*pwndata['cel'])
        phase = roi.extra['phase']
        pwnphase=roi.extra['pwnphase']

        print 'Making phaseogram'

        plot_kwargs = dict(ft1=ft1, skydir=pulsar_position, phase_range=phase, 
                           emin=pwnphase['optimal_emin'], emax=pwnphase['emax'], radius=pwnphase['optimal_radius'])
        plot_phaseogram(title='Phaseogram for %s' % name, filename='plots/phaseogram_%s.png' % name, **plot_kwargs)
        plot_phase_vs_time(title='Phase vs Time for %s' % name, filename='plots/phase_vs_time_%s.png' % name, **plot_kwargs)

        new_sources = roi.extra['new_sources']
        
        plots(roi, name, hypothesis, new_sources=new_sources, pulsar_position=pulsar_position,
              datadir=self.datadir, plotdir=self.plotdir)

    def variability_followup(self, hypothesis):
        name = self.name
        roi = self.reload_roi(hypothesis)

        roi.print_summary()
        roi.fit(use_gradient=False)
        roi.print_summary()

        frozen  = freeze_far_away(roi, roi.get_source(name).skydir, self.max_free)
        v = GtlikeVariabilityTester(roi,name, nbins=36, 
                              use_pointlike_ltcube=True, refit_background=True, refit_other_sources=True)
        v.plot(filename='plots/variability_%s_hypothesis_%s.pdf' % (name,hypothesis))
        unfreeze_far_away(roi, frozen)

        results = {hypothesis:{'variability':v.todict()}}
        savedict(results,'results_%s_variability_%s.yaml' % (name,hypothesis))

    def extul_followup(self, hypothesis):
        print 'Calculating extension upper limit'

        name = self.name
        roi = self.reload_roi(hypothesis)

        r=roi.extension_upper_limit(which=name, confidence=0.95, spatial_model=Gaussian)
        results = {hypothesis:{'pointlike':{'extension_upper_limit':r}}}
        savedict(results,'results_%s_extul_%s.yaml' % (name,hypothesis))
