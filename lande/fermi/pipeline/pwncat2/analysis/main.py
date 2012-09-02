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

def main_pipeline(args):
    do_at_pulsar = not args.no_at_pulsar
    do_point = not args.no_point
    do_extended = not args.no_extended 

    force_gradient(use_gradient=args.use_gradient)

    name=args.name
    emin=args.emin
    emax=args.emax

    pwnphase=yaml.load(open(args.pwnphase))[name]
    phase=pwnphase['phase']

    print 'phase = ',phase

    pwndata=yaml.load(open(args.pwndata))[name]

    savedir=None if args.no_savedir else join(os.getenv('PWD'),'savedir')

    results=r=defaultdict(lambda: defaultdict(dict))
    results['name']=name
    results['phase']=phase

    pointlike_kwargs=dict(name=name, max_free=args.max_free) 

    print 'Building the ROI'
    reg=PWNRegion(pwndata=args.pwndata, savedir=savedir)
    roi=reg.get_roi(name=name, phase=phase, 
                    catalog_kwargs=dict(free_radius=5, max_free=args.max_free),
                    fit_emin=emin, fit_emax=emax, 
                    binsperdec=args.binsperdec,
                    extended=False)

    modify = import_module(args.modify)
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
