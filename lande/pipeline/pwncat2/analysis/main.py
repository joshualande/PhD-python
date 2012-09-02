#!/usr/bin/env python

# This import has to come first
from analyze_helper import pointlike_analysis,gtlike_analysis,save_results,\
        plot_phaseogram,plot_phase_vs_time,all_energy,import_module
from uw.pulsar.phase_range import PhaseRange
from lande.fermi.likelihood.tools import force_gradient
from skymaps import SkyDir

from uw.like.SpatialModels import Gaussian

import os
from os.path import join
from glob import glob

from lande.fermi.likelihood.parlimits import all_params_limited

from setup_pwn import PWNRegion
from argparse import ArgumentParser
import yaml

from collections import defaultdict
import numpy as np
np.seterr(all='ignore')


def get_args():

    parser = ArgumentParser()
    parser.add_argument("--pwndata", required=True)
    parser.add_argument("-p", "--pwnphase")
    parser.add_argument("--no-phase-cut", default=False, action="store_true")
    parser.add_argument("-n", "--name", required=True, help="Name of the pulsar")
    parser.add_argument("--emin", default=1e2, type=float)
    parser.add_argument("--binsperdec", default=4, type=int)
    parser.add_argument("--emax", default=10**5.5, type=float)
    parser.add_argument("--use-gradient", default=False, action="store_true")
    parser.add_argument("--no-at-pulsar", default=False, action="store_true")
    parser.add_argument("--no-point", default=False, action="store_true")
    parser.add_argument("--no-extended", default=False, action="store_true")
    parser.add_argument("--no-cutoff", default=False, action="store_true")
    parser.add_argument("--no-savedir", default=False, action="store_true")
    parser.add_argument("--max-free", default=5, type=float)
    parser.add_argument("--modify", required=True)
    args=parser.parse_args()
    return args


if __name__ == '__main__':
    args = get_args()

    do_at_pulsar = not args.no_at_pulsar
    do_point = not args.no_point
    do_extended = not args.no_extended
    do_cutoff = not args.no_cutoff

    force_gradient(use_gradient=args.use_gradient)

    name=args.name
    emin=args.emin
    emax=args.emax

    pwnphase=yaml.load(open(args.pwnphase))[name]
    if args.no_phase_cut:
        phase = PhaseRange(0,1)
    else:
        phase=pwnphase['phase']

    print 'phase = ',phase

    pwndata=yaml.load(open(args.pwndata))[name]

    # nb, $PWD gets nicer paths then os.getcwd()
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
