import traceback
import sys

from uw.utilities import keyword_options
import numpy as np

from uw.like.roi_plotting import ROISlice
from skymaps import SkyImage

from uw.like.roi_state import PointlikeState

from . fit import fit_prefactor
from . tools import galstr


def paranoid_localize(roi, *args, **kwargs):

    state = PointlikeState(roi)
    try:
        roi.localize(*args, **kwargs)
    except Exception, ex:
        print 'ERROR localizing', ex
        traceback.print_exc(file=sys.stdout)
        state.restore()

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
        ('size',     2,     'size of image in degrees'),
        ('pixelsize',0.1,   'size, in degrees, of pixels'),
        ('galactic', False, 'galactic or equatorial coordinates'),
        ('proj',     'ZEA', 'projection name: can change if desired'),
        ('update',   False, 'Update the source of interest with the best fit'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, roi, which, **kwargs):
        keyword_options.process(self, kwargs)

        self.roi = roi
        self.source = roi.get_source(which)
        self.which = self.source.name

        self.__fill__()

    def __call__(self,skydir):
        roi = self.roi
        which = self.which

        roi.modify(which=which,skydir=skydir)

        # always start with same spectral model (better good for convergence)
        self.state.restore(just_spectra=True)

        def fit():
            try:
                # Fit just prefactor to begin with.
                # Helps with convergence
                fit_prefactor(roi,which)

                roi.fit(estimate_errors=False)
            except Exception, ex:
                print 'ERROR spectral fitting pointlike: ', ex
                traceback.print_exc(file=sys.stdout)

        fit()
        ll = -roi.logLikelihood(roi.parameters())

        best_model = roi.get_model(which).copy()

        self.state.restore(just_spectra=True)

        # Try a new model with a more standard starting value
        new_model = best_model.__class__()
        new_model.free = best_model.free.copy()

        roi.modify(which=which, model=new_model)
        fit()
        ll_alt = -roi.logLikelihood(roi.parameters())

        if ll_alt > ll:
            ll = ll_alt
            best_model = roi.get_model(which).copy()


        if not self.old_quiet: print '-- %s, ll=%.2f, ll-ll_0=%.2f' % (galstr(skydir), ll, ll-self.ll_0)


        return ll,best_model

    def __fill__(self):

        roi = self.roi
        self.state = PointlikeState(roi)

        self.ll_0 = -roi.logLikelihood(roi.parameters())

        self.old_quiet = roi.quiet
        roi.quiet = True

        self.init_skydir = self.source.skydir
        self.init_model = self.source.model.copy()

        if not self.old_quiet: print 'Grid localizing around initial position %s' % (galstr(self.init_skydir))

        # here, create skyimage

        self.skyimage = SkyImage(self.init_skydir, '', self.pixelsize, self.size, 1, self.proj, self.galactic, False)
        self.all_dirs = self.skyimage.get_wsdl()

        self.all_ll,self.all_models = zip(*[self(skydir) for skydir in self.all_dirs])

        self.state.restore()

        if not self.old_quiet: print 'Done Grid localizing, best position=%s, best ll=%.2f' % (galstr(self.best_position), self.best_logLikelihood-self.ll_0)

        if self.update:
            roi.modify(which=self.which, 
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
        

class MultiLocalizer():
    """ Object that can simultanously fit the poistion + extension
        of an arbitrary number of point and extended sources. """

    defaults = (
            ('tolerance',        0.01, "Fit tolerance to use when fitting"),
            ('verbose',          True, "Print more stuff during fit.")
    )

    @keyword_options.decorate(defaults)
    def __init__(self, roi, *sources, **kwargs):
        keyword_options.process(self, kwargs)

        print "Warning: this code doesn't really work"

        self.roi=roi

        if len(sources)<2:
            raise Exception("Must be passed in more than one source to localize.")

        for source in sources:
            try:
                roi.get_source(source)
            except:
                raise Exception("Unrecognized source %s" % source)

        self.sources = [roi.get_source(source) for source in sources]

    def update_roi(self,p):
        """ Takes in a list of paramters and updates all
            of the sources in the ROI. """
        for source in self.sources:
            antirotated=DualLocalizer.anti_rotate_equator(SkyDir(p[0],p[1]),self.middle)

            print source.name
            if isinstance(source,PointSource):
                p=p[2:]
                self.roi.modify(which=source,skydir=antirotated)
            else:
                spatial_model=self.roi.get_source(which=source).spatial_model
                nparam=len(spatial_model.p)
                spatial_model.set_parameters(p=p[2:nparam],center=antirotated,absolute=False)
                p=p[nparam:]
                self.roi.modify(which=source,spatial_model=spatial_model,preserve_center=False)
        assert(len(p)==0)

    @staticmethod
    def get_rotated_paramters(sources,middle):
        """ Returns a list of spatial paramters for all the sources (in the rotated
            coordiante system). """
        p=[]
        for source in sources:
            rotated=DualLocalizer.rotate_equator(source.skydir,middle)
            if isinstance(source,ExtendedSource):
                p += [rotated.ra(),rotated.dec()] + source.spatial_model.get_parameters(absolute=False)[2:].tolist()
            else:
                p += [rotated.ra(),rotated.dec()]
        return p

    def fit(self,p):
        print 'update roi'

        self.update_roi(p)

        print 'fit'

        ll=self.roi.fit(estimate_errors=False)

        if ll < self.ll_0:
            prev= [source.model.get_parameters() for source in self.sources]
            for source,p in zip(self.sources,self.init_spectral):
                source.model.set_parameters(p)

            ll_alt=self.roi.fit(estimate_errors=False)

            if ll_alt > ll: 
                ll=ll_alt
            else: 
                for source,p in zip(self.sources,prev):
                    source.model.set_parameters(p)

        if self.verbose: 
            print_str=[]
            for source in self.sources:
                if isinstance(source,PointSource):
                    print_str.append('%s: [%.3f,%.3f] f=%.1e' % \
                            (source.name,source.skydir.l(),source.skydir.b(),
                                DualLocalizer.print_flux(source,self.roi)))
                else:
                    print_str.append('%s: [%s] f=%.1e' % \
                            (source.name,source.spatial_model.full_spatial_string(),
                                DualLocalizer.print_flux(source,self.roi)))
            print_str.append('logL=%.3f' % ll)
            print_str.append('dlogL=%.3f' % (ll-self.ll_0))
            print ', '.join(print_str)

        return -ll # minimize negative log likelihood


    def localize(self):

        # Fit in a rotated coordinate system which is somewhere between all of
        # the sources being fit.
        self.middle = DualLocalizer.approx_mid_point(*[source.skydir for source in self.sources])


        self.init_spectral = [source.model.get_parameters() for source in self.sources]

        p0 = MultiLocalizer.get_rotated_paramters(self.sources,self.middle)

        old_quiet= self.roi.quiet
        self.roi.quiet=True

        steps=reduce(operator.add,
                [[0.1,0.1] if isinstance(source,PointSource) 
                    else source.spatial_model.get_steps().tolist()
                    for source in self.sources])

        self.ll_0=-1*self.roi.logLikelihood(self.roi.parameters())

        m = Minuit(self.fit, p0,
                   tolerance = self.tolerance,
                   maxcalls  = 500,
                   printMode = True, 
                   steps     = steps
                   )

        best_spatial,fval = m.minimize(method="SIMPLEX")

        self.roi.quiet = old_quiet

        return

