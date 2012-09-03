#!/usr/bin/env python
# This import has to come first
from . helper import pointlike_analysis,save_results,import_module, gtlike_analysis,plots,save_results,plot_phaseogram,plot_phase_vs_time
from lande.fermi.likelihood.tools import force_gradient

from uw.like.SpatialModels import Gaussian

import os
from os.path import join

from lande.fermi.likelihood.parlimits import all_params_limited

from setup_pwn import PWNRegion
import yaml

from collections import defaultdict


from uw.like.SpatialModels import Gaussian

from lande.fermi.likelihood.variability import VariabilityTester
from lande.fermi.likelihood.save import pointlike_dict_to_spectrum
from lande.utilities.save import loaddict
from lande.fermi.likelihood.free import freeze_far_away, unfreeze_far_away

from setup_pwn import load_pwn




class Pipeline(args):
    def __init__(self, **kwargs):
        self.__dict__.update(**kwargs)

    def main():
        do_at_pulsar = not self.no_at_pulsar
        do_point = not self.no_point
        do_extended = not self.no_extended 

        force_gradient(use_gradient=self.use_gradient)

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
                                                           cutoff=do_cutoff, **pointlike_kwargs)
        if do_point:
            r['point']['pointlike']=pointlike_analysis(roi, hypothesis='point', localize=True, 
                                                       cutoff=do_cutoff, 
                                                       model1 = model1,
                                                       **pointlike_kwargs)
        if do_extended:
            roi.modify(which=name, spatial_model=Gaussian(sigma=0.1), keep_old_center=True)

            r['extended']['pointlike']=pointlike_analysis(roi, hypothesis='extended', cutoff=False, 
                                                          fit_extension=True, 
                                                          **pointlike_kwargs)

        save_results(results,'results_%s_pointlike.yaml' % name)

    def followup(self, hypothesis, followup):

        do_cutoff = not self.no_cutoff

        force_gradient(use_gradient=self.use_gradient)

        name = self.name

        roi = load_pwn('roi_%s_%s.dat' % (hypothesis,name))


        if followup == 'gtlike':

            cutoff = do_cutoff and hypothesis in ['at_pulsar', 'point']
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
                                                          hypothesis=hypothesis, 
                                                          upper_limit=upper_limit,
                                                          cutoff=cutoff,
                                                          model1=model1,
                                                         )

            save_results(results,'results_%s_%s_%s.yaml' % (name,followup,hypothesis))

        if followup in ['plots','tsmaps']:

            if not os.path.exists('plots'): 
                os.makedirs('plots')

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
            overlay_kwargs = dict(pulsar_position=pulsar_position, new_sources=new_sources)
            
            if followup == 'plots':
                plots(roi, name, hypothesis, do_plots=True, do_tsmap=False, **overlay_kwargs)
            elif followup == 'tsmaps':
                plots(roi, name, hypothesis, do_plots=False, do_tsmap=True, **overlay_kwargs)


        elif followup == 'variability':

            roi.print_summary()
            roi.fit(use_gradient=False)
            roi.print_summary()

            frozen  = freeze_far_away(roi, roi.get_source(name).skydir, self.max_free)
            v = VariabilityTester(roi,name, nbins=36, 
                                  use_pointlike_ltcube=True, refit_background=True, refit_other_sources=True)
            v.plot(filename='plots/variability_%s_hypothesis_%s.pdf' % (name,hypothesis))
            unfreeze_far_away(roi, frozen)

            results = {hypothesis:{'variability':v.todict()}}
            save_results(results,'results_%s_%s_%s.yaml' % (name,followup,hypothesis))

        elif followup == 'extul':
            print 'Calculating extension upper limit'

            r=roi.extension_upper_limit(which=name, confidence=0.95, spatial_model=Gaussian)
            results = {hypothesis:{'pointlike':{'extension_upper_limit':r}}}
            save_results(results,'results_%s_%s_%s.yaml' % (name,followup,hypothesis))
