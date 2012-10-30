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

from uw.like.roi_save import load
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

from lande.fermi.pipeline.radiopsrs.data.loader import RadioPSRLoader

from . setup import RadioPSRROIBuilder
from . pointlike import pointlike_analysis

class Pipeline(object):

    @staticmethod
    def get_kwargs():

        parser = ArgumentParser()
        parser.add_argument("--radiopsr-data", required=True)
        parser.add_argument("--bigfile", required=True)
        parser.add_argument("--name", required=True, help="Name of the pulsar")
        args=parser.parse_args()
        return argparse_to_kwargs(args)


    def __init__(self, **kwargs):
        self.__dict__.update(**kwargs)

        self.input_kwargs = kwargs

        force_gradient(use_gradient=False)
        np.seterr(all='ignore')

        self.dirdict = dict(data='data',plots='plots',seds='seds')

        for dir in self.dirdict:
            if not os.path.exists(dir): os.makedirs(dir)

        self.radiopsr_loader = RadioPSRLoader(self.radiopsr_data, self.bigfile)

    def main(self):

        name=self.name

        rb = RadioPSRROIBuilder(radiopsr_loader=self.radiopsr_loader)
        roi = rb.build_roi(name=name)

        results = self.input_kwargs

        savedict(results,'results_%s_general.yaml' % name)

        hypothesis = 'at_pulsar'
        results=pointlike_analysis(roi, name, hypothesis=hypothesis, dirdict=self.dirdict)
        savedict(results,'results_%s_pointlike_%s.yaml' % (name,hypothesis))

        hypothesis = 'point'
        results=pointlike_analysis(roi, name, hypothesis=hypothesis, localize=True, dirdict=self.dirdict)
        savedict(results,'results_%s_pointlike_%s.yaml' % (name,hypothesis))

        hypothesis = 'extended'
        roi.modify(which=name, spatial_model=Gaussian(sigma=0.1), keep_old_center=True)

        results=pointlike_analysis(roi, name, hypothesis=hypothesis, fit_extension=True, dirdict=self.dirdict)
        savedict(results,'results_%s_pointlike_%s.yaml' % (name,hypothesis))

    def reload_roi(self,hypothesis, *args, **kwargs):
        roi = load('roi_%s_%s.dat' % (hypothesis,name), *args, **kwargs)
        roi.print_summary(galactic=True, maxdist=10)
        return roi