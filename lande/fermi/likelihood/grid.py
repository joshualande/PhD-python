import numpy as np
import pylab as P

from uw.utilities import keyword_options
from uw.utilities.parmap import LimitMapper
from uw.like.roi_state import PointlikeState

from . basefit import BaseFitter
from . save import source_dict


class SpectralFitLimited(BaseFitter):

    defaults = BaseFitter.defaults + (
        ('energy_units', 'MeV', 'default units to plot energy flux (y axis) in.'),
        ('flux_units',  'erg', 'default units to plot energy (x axis) in'),
        ('keep_best', True, "keep the best fit"),
        ('fit_kwargs', dict(use_gradient=False), 'kwargs past into roi.fit'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, roi, name, param_name, param_min, param_max, **kwargs):

        self.roi = roi
        self.name = name
        self.param_name = param_name
        self.param_min = param_min
        self.param_max = param_max

        keyword_options.process(self, kwargs)

        self._calculate()

    def _calculate(self):

        self.results=dict()

        self.init_state = PointlikeState(self.roi)

        model = self.roi.get_model(self.name)
        init_mapper = model.get_mapper(self.param_name)
        param_val = model[self.param_name]

        assert param_val >= self.param_min and param_val <= self.param_max

        if isinstance(init_mapper,LimitMapper):
            assert init_mapper.min <= self.param_min and init_mapper.max >= self.param_max


        model.set_mapper(self.param_name,LimitMapper(self.param_min, self.param_max, scale=param_val))

        self.roi.modify(which=self.name, model=model)

        if self.verbosity:
            print 'Before fit'
            self.roi.print_summary()

        self.results['fit_before']=source_dict(self.roi, self.name, energy_units=self.energy_units, flux_units=self.flux_units)

        self.roi.fit(**self.fit_kwargs)

        if self.verbosity:
            print 'After fit'
            self.roi.print_summary()

        self.results['fit_after']=source_dict(self.roi, self.name, energy_units=self.energy_units, flux_units=self.flux_units)

        if self.keep_best:
            model = self.roi.get_model(self.name)
            model.set_mapper(self.param_name,init_mapper)
        else:
            self.init_state.restore(just_spectra=True)


class SpectralGrid(BaseFitter):
    """ Perform a grid over parameters. """

    defaults = BaseFitter.defaults + (
        ('energy_units', 'MeV', 'default units to plot energy flux (y axis) in.'),
        ('flux_units',  'erg', 'default units to plot energy (x axis) in'),
        ('param_vals', None, "List of parameters to grid over. If specified, don't set param_min, param_max, or nparams"),
        ('param_min', None, "min parameter. If set, don't specify param_vals"),
        ('param_max', None, "max parameter. If set, don't specify param_vals"),
        ('nparams', None, "Number of params in grid. If set, don't specify param_vals"),
        ('keep_best', True, "keep the best fit"),
        ('fit_kwargs', dict(use_gradient=False), 'kwargs past into roi.fit'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, roi, name, param_name, **kwargs):

        self.roi = roi
        self.name = name
        self.param_name = param_name

        keyword_options.process(self, kwargs)

        if self.param_vals is not None:
            assert self.param_min is None and self.param_max is None and self.nparams is None
        else:
            assert self.param_min is not None and self.param_max is not None and self.nparams is not None
            self.param_vals = np.linspace(self.param_min, self.param_max, self.nparams)

        self._calculate()

    def _calculate(self):

        roi = self.roi
        name = self.name
        param_name = self.param_name

        self.init_state = PointlikeState(roi)

        self.results = dict(
            name=name,
            param_name=param_name,
            param_vals=self.param_vals,
                            grid=[])

        if self.verbosity:
            print 'Performing grid over parameter %s for source %s' % (name, param_name)

        best_state = None
        self.best_ll = -np.inf

        model = roi.get_model(which=name)
        old_free = model.get_free(param_name)

        for i,p in enumerate(self.param_vals):
            if self.verbosity:
                print 'looping for param %s=%s (%d/%d)' % (param_name, p, i+1,len(self.param_vals))

            self.init_state.restore(just_spectra=True)

            model = roi.get_model(which=name)
            model[param_name]=p
            model.set_free(param_name,False)

            roi.modify(which=name, model=model)

            if self.verbosity:
                roi.print_summary()
            roi.fit(**self.fit_kwargs)
            if self.verbosity:
                roi.print_summary()

            d=source_dict(roi,name, energy_units=self.energy_units, flux_units=self.flux_units)
            self.results['grid'].append(d)

            ll = self.results['grid'][-1]['logLikelihood']

            if ll > self.best_ll:
                self.best_state = PointlikeState(roi)
                self.best_ll = ll
                self.best_d = d

        self.results['best'] = self.best_d

        if self.keep_best:
            self.best_state.restore(just_spectra=True)
            model = roi.get_model(which=name)
            model.set_free(param_name,old_free)
            roi.modify(which=name, model=model)
        else:
            self.init_state.restore(just_spectra=True)

    def plot(self, filename):

        param_vals = self.results['param_vals']
        ll = [i['logLikelihood'] for i in self.results['grid']]

        P.plot(param_vals,ll)
        P.ylabel('logLikelihood')
        P.xlabel(self.results['param_name'])
        P.savefig(filename)


