import traceback
import sys

from uw.utilities import keyword_options
import numpy as np

from skymaps import SkyDir
from skymaps import SkyImage

from uw.like.roi_state import PointlikeState
from uw.like.roi_localize import DualLocalizer
from uw.utilities import rotations

from . fit import fit_prefactor
from . tools import galstr
from . save import logLikelihood,skydirdict,ts_dict
from . base import BaseFitter


def paranoid_localize(roi, name, verbosity=True):

    state = PointlikeState(roi)
    try:
        if verbosity: print 'Trying ROILocalizer'
        roi.localize(which=name, update=True, verbose=verbosity)
        ellipse = roi.get_ellipse()
        ellipse['method'] = 'ROILocalizer'

        roi.get_source(which=name).localization = ellipse
        return ellipse

    except Exception, ex:
        print 'ERROR with ROILocalizer', ex
        traceback.print_exc(file=sys.stdout)
        state.restore()

        try:
            if verbosity: print 'Trying MinuitLocalizer:'
            m=MinuitLocalizer(roi,name, verbosity=True)
            ellipse = m.todict()
        except Exception, ex:
            print 'ERROR with MinuitLocalizer', ex
            traceback.print_exc(file=sys.stdout)
            state.restore()

            ellipse = dict(method='Failed')
            roi.get_source(which=name).localization = ellipse

        return ellipse



class GridLocalize(object):
    """ Simple class evalulates the TS of a source
        along a grid in position (a TS map) and find the 
        maximum TS.
    
        This code can e useful because it might avoid false minima 
        and help converge on the true minimum.
        
        In addition, during the spectral analysis, this function
        will try an alternate default spectral parmaeters for the
        source of interest, so it may be a ltitle more robust
        at really finding the best spectral parmaeters for the source
        """


    defaults = (
        ('size',         2, 'size of image in degrees'),
        ('pixelsize',  0.1, 'size, in degrees, of pixels'),
        ('galactic', False, 'galactic or equatorial coordinates'),
        ('proj',     'ZEA', 'projection name: can change if desired'),
        ('update',    True, 'Update the source of interest with the best fit'),
        ('verbosity',   True, "Print more stuff during fit."),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, roi, name, **kwargs):
        keyword_options.process(self, kwargs)

        self.roi = roi
        self.name = name

        self.__fill__()

    def __call__(self,skydir):
        roi = self.roi
        name = self.name

        roi.modify(which=name,skydir=skydir)

        # always start with same spectral model (better good for convergence)
        self.state.restore(just_spectra=True)

        def fit():
            try:
                # Fit just prefactor to begin with.
                # Helps with convergence
                fit_prefactor(roi,name)

                roi.fit(estimate_errors=False)
            except Exception, ex:
                print 'ERROR spectral fitting pointlike: ', ex
                traceback.print_exc(file=sys.stdout)

        fit()
        ll = logLikelihood(roi)

        best_model = roi.get_model(name).copy()

        self.state.restore(just_spectra=True)

        # Try a new model with a more standard starting value
        # This implementation is better than previous because
        # it works for FileFunction
        new_model = best_model.copy()
        new_model.set_all_parameters(new_model.default_p)

        roi.modify(which=name, model=new_model)
        fit()
        ll_alt = logLikelihood(roi)

        if ll_alt > ll:
            ll = ll_alt
            best_model = roi.get_model(name).copy()


        if self.verbosity: print '-- %s, ll=%.2f, ll-ll_initial=%.2f' % (galstr(skydir), ll, ll-self.ll_initial)


        return ll,best_model

    def __fill__(self):

        roi = self.roi
        name = self.name
        source = roi.get_source(name)

        self.state = PointlikeState(roi)

        self.ll_initial = logLikelihood(roi)

        self.old_quiet = roi.quiet
        roi.quiet = True

        self.init_skydir = source.skydir
        self.init_model = source.model.copy()

        if self.verbosity: 
            print 'Grid localizing around initial position %s' % (galstr(self.init_skydir))

        # here, create skyimage

        self.skyimage = SkyImage(self.init_skydir, '', self.pixelsize, self.size, 1, self.proj, self.galactic, False)
        self.all_dirs = self.skyimage.get_wsdl()

        self.all_ll,self.all_models = zip(*[self(skydir) for skydir in self.all_dirs])

        self.state.restore()

        if self.verbosity: 
            print 'Done Grid localizing, best position=%s, best ll=%.2f' % (galstr(self.best_position), self.best_logLikelihood-self.ll_initial)

        if self.update:
            roi.modify(which = name, 
                       skydir = self.best_position,
                       model = self.best_model)

    @property
    def best_position(self):
        return self.all_dirs[int(np.argmax(self.all_ll))]

    @property
    def best_logLikelihood(self):
        return np.max(self.all_ll)

    @property
    def best_model(self):
        return self.all_models[int(np.argmax(self.all_ll))]
        

class MinuitLocalizer(BaseFitter):
    """ Fit two point sources at the same time. Fit the center of position
        and relative difference since they are more robust parameters.

        This method is only suitable for sources relativly nearby (~<1 degree).

        Note, bandfits is not allowed for fitting because the bandfits
        algorithm does not really work when fitting two really nearby
        sources. """

    defaults = (
        ('fit_kwargs', dict(), 'kwargs into fit function'),
        ('tolerance',    0.01, "Fit tolerance to use when fitting"),
        ('verbosity',   True, "Print more stuff during fit."),
        ('update',      True,"Update the source position after localization"),
    )


    @keyword_options.decorate(defaults)
    def __init__(self, roi, name, **kwargs):
        keyword_options.process(self, kwargs)

        self.roi = roi
        self.name = name

        self.localize()

    def fit(self,p):
        roi=self.roi
        name = self.name
        source = roi.get_source(name)

        x,y = p
        s = SkyDir(x,y)
        rot_back = rotations.anti_rotate_equator(s,self.start)

        self.init_state.restore(just_spectra=True)
        roi.modify(which=self.name,skydir=rot_back)

        roi.fit(estimate_errors=False, **self.fit_kwargs)
        ll=logLikelihood(roi)

        if self.verbosity: print 'd=%s f=%.1e, dist=%.3f logL=%.3f dlogL=%.3f' % \
                (galstr(rot_back), DualLocalizer.print_flux(source,roi),
                 np.degrees(rot_back.difference(self.start)),
                 ll,ll-self.ll_initial)

        return -ll # minimize negative log likelihood

    def localize(self):
        roi=self.roi
        name=self.name

        if self.verbosity: print 'Using MinuitLocalizer to localize source %s' % name

        self.start = roi.get_source(name).skydir

        self.ll_initial=ll_initial=logLikelihood(roi)
        TS_initial=ts_dict(roi,name)
        self.init_state = PointlikeState(roi)

        old_quiet=roi.quiet
        roi.quiet=True

        from uw.utilities.minuit import Minuit
        steps=[0.1,0.1] # expect to fit sources ~ 0.1 degrees away.
        p0 = [0,0]
        m = Minuit(self.fit,
                   p0,
                   tolerance = self.tolerance,
                   maxcalls  = 500,
                   limits    = [[-1,1],[-1,1]],
                   printMode = True, 
                   steps     = steps)

        best_p,fval = m.minimize(method="SIMPLEX")

        if self.verbosity: print 'Computing covariance matrix for source %s' % name 
        self.fit(best_p)
        cov_matrix = m.errors(method="HESSE")

        if self.verbosity: print 'Fitting one last time with best parameters'
        self.fit(best_p)
        best_spatial = roi.get_source(name).skydir
        TS_best = ts_dict(roi,name)
        ll_best = logLikelihood(roi)

        self.results = dict(
            ll_initial = ll_initial,
            ll_best = ll_best,
            delta_TS = 2*(ll_best-ll_initial),
            position = skydirdict(best_spatial),
            TS_best = TS_best,
            TS_initial = TS_initial,
            method = 'MinuitLocalizer'
        )

        sigma_x = np.sqrt(cov_matrix[0,0])
        sigma_y = np.sqrt(cov_matrix[1,1])
        self.results['position_error'] = np.sqrt(sigma_x*sigma_y)

        roi.quiet = old_quiet

        if not self.update: self.init_state.restore()

        # save localization into source.
        roi.get_source(which=name).localization = self.results


if __name__ == "__main__":
    import doctest
    doctest.testmod()
