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

from lande.fermi.pipeline.gamma_quiet_psrs.data.loader import RadioPSRLoader

from . setup import RadioPSRROIBuilder
from . pointlike import pointlike_analysis
from . plots import smooth_plots
from . tsmaps import tsmap_plots
from . gtlike import gtlike_analysis

class Pipeline(object):

    @staticmethod
    def get_kwargs():
        parser = ArgumentParser()
        parser.add_argument("--radiopsr-data", required=True)
        parser.add_argument("--bigfile", required=True)
        parser.add_argument("--name", required=True, help="Name of the pulsar")
        parser.add_argument("--fast", default=False, action='store_true')
        parser.add_argument("--cachedata", default=False, action='store_true')
        parser.add_argument("--no-point", default=False, action='store_true')
        parser.add_argument("--no-extended", default=False, action='store_true')
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
        roi = rb.build_roi(name=name, fast=self.fast)

        results = self.input_kwargs

        savedict(results,'results_%s_general.yaml' % name)

        hypothesis = 'at_pulsar'
        results=pointlike_analysis(self, roi, name, hypothesis=hypothesis, upper_limit=True)
        savedict(results,'results_%s_pointlike_%s.yaml' % (name,hypothesis))

        if not self.no_point:
            hypothesis = 'point'
            results=pointlike_analysis(self, roi, name, hypothesis=hypothesis, localize=True)
            savedict(results,'results_%s_pointlike_%s.yaml' % (name,hypothesis))

        if not self.no_extended:
            hypothesis = 'extended'
            roi.modify(which=name, spatial_model=Gaussian(sigma=0.1), keep_old_center=True)
            results=pointlike_analysis(self, roi, name, hypothesis=hypothesis, fit_extension=True)
            savedict(results,'results_%s_pointlike_%s.yaml' % (name,hypothesis))

    def reload_roi(self,hypothesis, *args, **kwargs):
        roi = load('roi_%s_%s.dat' % (hypothesis,self.name), *args, **kwargs)
        roi.print_summary(galactic=True, maxdist=10)
        return roi

    def gtlike_followup(self, hypothesis):
        roi = self.reload_roi(hypothesis)
        results=gtlike_analysis(self, roi, self.name, hypothesis, upper_limit=hypothesis=='at_pulsar')
        savedict(results,'results_%s_gtlike_%s.yaml' % (self.name,hypothesis))

    def plots_followup(self, hypothesis):
        roi = self.reload_roi(hypothesis)
        smooth_plots(self, roi, self.name, hypothesis, size=5)
        smooth_plots(self, roi, self.name, hypothesis, size=10)

    def tsmaps_followup(self, hypothesis):
        roi = self.reload_roi(hypothesis)
        tsmap_plots(self, roi, self.name, hypothesis, size=5)
        tsmap_plots(self, roi, self.name, hypothesis, size=10)

    def irfs_systematics_followup(self, hypothesis):
        # use bracketing irfs
        pass

    def alterante_diffuse_systematics_followup(self, hypothesis):
        # use bracketing irfs
        pass
