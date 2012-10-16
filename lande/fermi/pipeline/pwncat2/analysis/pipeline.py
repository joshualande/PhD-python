#!/usr/bin/env python
# This import has to come first
from lande.fermi.likelihood.roi_gtlike import Gtlike

import sys
import os
from os.path import join
from collections import defaultdict
from argparse import ArgumentParser
import traceback

import yaml
import numpy as np

from skymaps import SkyDir

from uw.like.SpatialModels import Gaussian
from uw.like.SpatialModels import Gaussian
from uw.utilities.parmap import LogMapper


from lande.utilities.save import loaddict,savedict
from lande.utilities.load import import_module
from lande.utilities.argumentparsing import argparse_to_kwargs

from lande.fermi.pulsar.plotting import plot_phaseogram, plot_phase_vs_time
from lande.fermi.likelihood.tools import force_gradient
from lande.fermi.likelihood.parlimits import all_params_limited
from lande.fermi.likelihood.variability import CombinedVariabilityTester
from lande.fermi.likelihood.load import pointlike_dict_to_spectrum
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
        parser.add_argument("--fast", default=False, action="store_true")
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

        if not self.fast:
            emin=self.emin
            emax=self.emax
            binsperdec=self.binsperdec
            free_radius=5
            roi_size=10
            max_free=self.max_free
        else:
            emin=1e4
            emax=1e5
            binsperdec=2
            free_radius=2
            roi_size=5
            max_free=2

        pwnphase=yaml.load(open(self.pwnphase))[name]
        phase=pwnphase['phase']

        print 'phase = ',phase

        pwndata=yaml.load(open(self.pwndata))[name]

        savedir=None if self.no_savedir else join(os.getenv('PWD'),'savedir')

        results=r=defaultdict(lambda: defaultdict(dict))
        results['name']=name
        results['phase']=phase

        pointlike_kwargs=dict(name=name, max_free=max_free) 

        print 'Building the ROI'
        reg=PWNRegion(pwndata=self.pwndata, savedir=savedir)
        roi=reg.get_roi(name=name, phase=phase, 
                        catalog_kwargs=dict(free_radius=free_radius, max_free=max_free),
                        fit_emin=emin, fit_emax=emax, 
                        binsperdec=binsperdec,
                        roi_size=roi_size, 
                        extended=False)

        modify = import_module(self.modify)
        roi.extra['phase'] = phase
        roi.extra['new_sources'] = modify.modify_roi(name,roi)
        roi.extra['pwnphase'] = pwnphase

        savedict(results,'results_%s_general.yaml' % name)

        assert all_params_limited(roi, except_sources=[name])
        mapper0 = roi.get_model(which=name).get_mapper(0)
        assert isinstance(mapper0,LogMapper) or (type(mapper0)==type and issubclass(mapper0,LogMapper))

        cutoff_model=modify.get_cutoff_model(name)
        override_localization=modify.get_override_localization(name)

        if do_at_pulsar:
            hypothesis = 'at_pulsar'
            results=pointlike_analysis(roi, hypothesis=hypothesis, 
                                       seddir=self.seddir, datadir=self.datadir, plotdir=self.plotdir,
                                       cutoff=do_cutoff, **pointlike_kwargs)
            savedict(results,'results_%s_pointlike_%s.yaml' % (name,hypothesis))

        if do_point:
            hypothesis = 'point'
            results=pointlike_analysis(roi, hypothesis=hypothesis, localize=True, 
                                       seddir=self.seddir, datadir=self.datadir, plotdir=self.plotdir,
                                       cutoff=do_cutoff, 
                                       cutoff_model = cutoff_model,
                                       override_localization=override_localization,
                                       **pointlike_kwargs)
            savedict(results,'results_%s_pointlike_%s.yaml' % (name,hypothesis))
        if do_extended:
            hypothesis = 'extended'
            roi.modify(which=name, spatial_model=Gaussian(sigma=0.1), keep_old_center=True)

            results=pointlike_analysis(roi, hypothesis=hypothesis, cutoff=False, 
                                       seddir=self.seddir, datadir=self.datadir, plotdir=self.plotdir,
                                       fit_extension=True, 
                                       **pointlike_kwargs)
            savedict(results,'results_%s_pointlike_%s.yaml' % (name,hypothesis))

    def reload_roi(self,hypothesis, *args, **kwargs):
        name = self.name
        roi = load_pwn('roi_%s_%s.dat' % (hypothesis,name), *args, **kwargs)

        roi.print_summary(galactic=True, maxdist=10)

        return roi

    def gtlike_followup(self, hypothesis):

        name = self.name
        roi = self.reload_roi(hypothesis)

        cutoff = (not self.no_cutoff) and hypothesis in ['at_pulsar', 'point']
        upper_limit = hypothesis=='at_pulsar'
        if cutoff:
            pointlike_results = loaddict('results_%s_pointlike_%s.yaml' % (name,hypothesis))
            cutoff_model=pointlike_results['test_cutoff']['hypothesis_1']['spectrum']
            cutoff_model=pointlike_dict_to_spectrum(cutoff_model)
            cutoff_model.set_default_limits(oomp_limits=True)
        else:
            cutoff_model=None

        results=gtlike_analysis(roi, name=name,
                                max_free = self.max_free,
                                seddir=self.seddir, datadir=self.datadir, plotdir=self.plotdir,
                                hypothesis=hypothesis, 
                                upper_limit=upper_limit,
                                cutoff=cutoff,
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

        modify = import_module(self.modify)
        good_interval = modify.get_variability_time_cuts(name)

        nbins=36

        if good_interval is not None:
            ft1files=roi.sa.pixeldata.ft1files
            earliest_time, latest_time = CombinedVariabilityTester.get_time_range(ft1files)
            bins = b = np.round(np.linspace(earliest_time, latest_time, nbins+1)).astype(int)
            starts  = b[:-1]
            stops = b[1:]

            print 'Initial binning:'
            print ' * starts=',starts
            print ' * stops=',stops

            starts,stops = zip(*[(start,stop) \
                                for (start,stop) in zip(starts,stops) \
                                if good_interval(start,stop)])

            print 'Initial binning:'
            print ' * starts=',starts
            print ' * stops=',stops

            kwargs=dict(tstarts=starts, tstops=stops)
        else:
            kwargs=dict(nbins=nbins)


        frozen  = freeze_far_away(roi, roi.get_source(name).skydir, self.max_free)
        v = CombinedVariabilityTester(roi,name, 
                                      use_pointlike_ltcube=True, refit_background=True, 
                                      refit_other_sources=True,
                                      verbosity=4, **kwargs)
        unfreeze_far_away(roi, frozen)

        results = v.todict()
        savedict(results,'results_%s_variability_%s.yaml' % (name,hypothesis))

        try:
            v.plot(filename='plots/variability_%s_hypothesis_%s.pdf' % (name,hypothesis))
        except Exception, ex:
            print 'ERROR plotting variability tester:', ex
            traceback.print_exc(file=sys.stdout)

    def extul_followup(self, hypothesis):
        print 'Calculating extension upper limit'

        name = self.name
        roi = self.reload_roi(hypothesis)

        r=roi.extension_upper_limit(which=name, confidence=0.95, spatial_model=Gaussian)
        results = {hypothesis:{'pointlike':{'extension_upper_limit':r}}}
        savedict(results,'results_%s_extul_%s.yaml' % (name,hypothesis))
