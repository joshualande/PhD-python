#!/usr/bin/env python
# This import has to come first
from . helper import pointlike_analysis,save_results,import_module
from lande.fermi.likelihood.tools import force_gradient

from uw.like.SpatialModels import Gaussian

import os
from os.path import join

from lande.fermi.likelihood.parlimits import all_params_limited

from setup_pwn import PWNRegion
import yaml

from collections import defaultdict


class Pipeline(args):
    def __init__(self,
        name,
        pwndata, pwnphase, modify,
        emin, emax, binsperdec,
        use_gradient,
        no_at_pulsar, no_point, no_extended, no_cutoff, no_savedir,
        max_free):

        self.name=name
        self.pwndata=pwndata
        self.pwnphase=pwnphase
        self.emin=emin
        self.emax=emax
        self.binsperdec=binsperdec
        self.use_gradient=use_gradient
        self.no_at_pulsar=no_at_pulsar
        self.no_point=no_point
        self.no_extended=no_extended
        self.no_cutoff=no_cutoff
        self.no_savedir=no_savedir
        self.max_free=max_free
        self.modify=modify

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

        pointlike_kwself=dict(name=name, max_free=self.max_free) 

        print 'Building the ROI'
        reg=PWNRegion(pwndata=self.pwndata, savedir=savedir)
        roi=reg.get_roi(name=name, phase=phase, 
                        catalog_kwself=dict(free_radius=5, max_free=self.max_free),
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
                                                           cutoff=do_cutoff, **pointlike_kwself)
        if do_point:
            r['point']['pointlike']=pointlike_analysis(roi, hypothesis='point', localize=True, 
                                                       cutoff=do_cutoff, 
                                                       model1 = model1,
                                                       **pointlike_kwself)
        if do_extended:
            roi.modify(which=name, spatial_model=Gaussian(sigma=0.1), keep_old_center=True)

            r['extended']['pointlike']=pointlike_analysis(roi, hypothesis='extended', cutoff=False, 
                                                          fit_extension=True, 
                                                          **pointlike_kwself)

        save_results(results,'results_%s_pointlike.yaml' % name)
