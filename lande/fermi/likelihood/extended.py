""" 
    Code to deal with extended sources in pointlike.
"""
from os.path import expandvars

import numpy as np
import yaml
import pywcs
import pyfits

from skymaps import SkyDir

from uw.like.roi_plotting import DegreesFormatter
from uw.like.roi_state import PointlikeState
from uw.like.roi_image import ModelImage
from uw.like.roi_monte_carlo import FitsShrinker
from uw.utilities import keyword_options

from lande.utilities.tools import tolist
from lande.utilities.fits import expand_fits_header

import pylab as P


class ExtensionProfile(object):
    defaults = (
        ("num_points",      15, "Number of poitns to calcualte profile over"),
        ("lower_limit",   None, "smallest sigma for the profile. Default is very smalle xtension."),
        ("upper_limit",   None, "largest sigma for the profile"),
        ("quick",         None, "Default way to calculate TS. Default is taken from roi.TS()"),
        ("use_gradient",  None, "use analytic gradient during spectral fit. Default is taken from roi.fit()."),
        ("fignum",        None, "passed to matplotlib."),
        ("figsize",      (4,4), "size of plot, in inches."),
    )

    @keyword_options.decorate(defaults)
    def __init__(self,roi,which,**kwargs):
        """ Object for calculating TS as a function of sigma. """


        self.roi = roi
        self.which = which

        keyword_options.process(self, kwargs)

        self.source = roi.get_source(which)

        if not hasattr(self.source,'spatial_model'):
            raise Exception("An extension profile can only be calculated for extended sources")

        self.spatial_model=self.source.spatial_model
        if not len(self.spatial_model.p)==3: 
            raise Exception("An extension profile can only be calculated for extended sources with 3 parameters (position + one extension)")

        self.fit_kwargs = dict(estimate_errors=False)
        if self.use_gradient is not None: self.fit_kwargs['use_gradient']=self.use_gradient

        self.ts_kwargs = dict(which=self.which)
        if self.quick is not None: self.ts_kwargs['quick']=self.quick

        self.fill()

    def fill(self):

        roi = self.roi

        state = PointlikeState(roi)

        if not roi.quiet: print 'Calculating extension profile for %s' % self.source.name

        init_p = roi.get_parameters().copy()

        # Keep the TS function quiet
        old_quiet = roi.quiet
        roi.quiet=True

        sigma = self.spatial_model['sigma']
        sigma_err = self.spatial_model.error('sigma')

        upper_limit = min(sigma + max(3*sigma_err,sigma),3) if self.upper_limit is None else self.upper_limit

        # make the bottom point ~ 0.1xfirst point
        lower_limit = float(upper_limit)/self.num_points/10.0 if self.lower_limit is None else self.lower_limit

        self.extension_list=np.linspace(lower_limit,upper_limit,self.num_points)

        self.TS_spectral=np.empty_like(self.extension_list)
        self.TS_bandfits=np.empty_like(self.extension_list)

        roi.setup_energy_bands()

        if not old_quiet: print '%20s %20s %20s' % ('sigma','TS_spectral','TS_bandfits')
        for i,sigma in enumerate(self.extension_list):
            roi.modify(which=self.which, sigma=sigma)

            roi.fit(**self.fit_kwargs)

            params=roi.parameters()
            ll_a=-1*roi.logLikelihood(roi.parameters())

            roi.update_counts(init_p)
            roi.fit(**self.fit_kwargs)
            ll_b=-1*roi.logLikelihood(roi.parameters())
            if ll_a > ll_b: roi.update_counts(params)

            self.TS_spectral[i]=roi.TS(**self.ts_kwargs)
            self.TS_bandfits[i]=roi.TS(bandfits=True,**self.ts_kwargs)

            if not old_quiet: print 'sigma=%.2f ts_spec=%.1f, ts_band=%.1f' % (sigma, self.TS_spectral[i],self.TS_bandfits[i])
        
        state.restore()


    def todict(self):
        d=dict(sigma=self.extension_list,
               TS_spectral=self.TS_spectral,
               TS_bandfits=self.TS_bandfits)

        return tolist(d)

    def save(self, filename):
        file=open(filename,'w')
        file.write(yaml.dump(self.todict()))

    def plot(self,filename):

        fig = P.figure(self.fignum,self.figsize)
        axes = fig.add_subplot(111)
        P.clf()
        P.plot(self.extension_list,self.TS_spectral)
        P.xlabel('Extension')
        P.ylabel('Test Statistic')
        P.gca().xaxis.set_major_formatter(DegreesFormatter)


        P.title('Extension profile %s' % self.source.name)

        P.savefig(expandvars(filename))


class TSExtVsEnergy(object):

    """ Object to make a plot of TSext vs energy. """

    defaults = (
            ('which',           None,  'Source to analyze'),
            ('title',           None, 'Title for the plot'),
            ('fignum',             1,                   ''),
            ('figsize',        (6,4),                   ''),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, roi, **kwargs):
        keyword_options.process(self, kwargs)

        roi.setup_energy_bands()
        self.emin = np.asarray([eb.emin for eb in roi.energy_bands])
        self.emax = np.asarray([eb.emax for eb in roi.energy_bands])

        which = self.which
        old_roi_p   = roi.get_parameters().copy()

        # extended hypothesis

        source=roi.get_source(which='IC443')
        sm = source.spatial_model
        manager,index=roi.mapper(which)
        roi.fit(estimate_errors=False)

        self.ll_ext,self.ll_pt = [],[]


        for eb in roi.energy_bands:
            self.ll_ext.append(
                -sum(band.logLikelihood() for band in eb.bands)
                )


        sm.shrink()
        manager.bgmodels[index].initialize_counts(roi.bands)
        roi.__update_state__()

        roi.fit(estimate_errors=False)

        # point hypothesis

        manager,index=roi.mapper('IC443')
        for eb in roi.energy_bands:
            self.ll_pt.append(
                -sum(band.logLikelihood() for band in eb.bands)
                )

        sm.unshrink()
        manager.bgmodels[index].initialize_counts(roi.bands)

        roi.set_parameters(old_roi_p)
        roi.__update_state__()

        self.ll_ext = np.asarray(self.ll_ext)
        self.ll_pt = np.asarray(self.ll_pt)

        self.ts_ext=2*(self.ll_ext-self.ll_pt)
        self.ts_ext[self.ts_ext<0]=0

    def show(self,filename=None,axes=None):

        if axes is None:
            fig = P.figure(self.fignum,self.figsize)
            axes = fig.get_axes([0.15,0.15,0.7,0.7])
            P.clf()

        ax = self.axes = axes

        ax.set_xlabel(r'Energy (MeV)')
        ax.set_ylabel(r'$\mathrm{TS}_\mathrm{ext}$')
        ax.semilogx(self.emin,self.ts_ext,'k',drawstyle='steps-post')

        if self.title is None: 
            self.title = '$\mathrm{TS}_\mathrm{ext}$ vs Energy'
            self.title += ' for %s' % self.source.name

        ax.set_title(self.title)

        if filename is not None:
            P.savefig(filename)



def save_template_from_header(spatial_model, header, filename, factor=1, **kwargs):
    """ Saves an extended source to a file with the same
        projection and binning as an input pyfits header.

        factor will cause image to be sampled on a grid factor**2 and then downsampled.
        This is useful for function which vary rapidly over the pixels
        of the files file.

        kwargs can be used to pass clobber into writeto
    """
    assert abs(header['CDELT1']) == abs(header['CDELT2'])

    if header['naxis'] == 3:
        header = FitsShrinker.convert_header_2d(header)

    expanded_header = expand_fits_header(header, factor)

    lons,lats = FitsShrinker.pyfits_to_sky(expanded_header)

    galactic = FitsShrinker.is_header_galactic(header)
    coordsys = SkyDir.GALACTIC if galactic else SkyDir.EQUATORIAL

    it = np.nditer([lons, lats, None])
    for lon,lat,d in it:
        skydir=SkyDir(float(lon),float(lat),coordsys)
        d[...] = spatial_model(skydir)
    data = it.operands[2]

    data = ModelImage.downsample(data, factor)

    pixelsize = header['CDELT1']
    normalization = np.sum(data)*np.radians(pixelsize)**2
    normed_data = data/normalization

    hdu = pyfits.PrimaryHDU(header=header, data=data)

    hdulist = pyfits.HDUList([hdu])
    hdulist.writeto(filename, **kwargs)


def fit_extension_fast(roi, which, *args, **kwargs):
    """ Perform an extension fit with all other
        localizable sources in the ROI frozen. 
        
        This implementation will nevertheless
        fit the diffuse background and the source
        of interest during the extension fit 
        so should be a bit faster without too much
        loss of accuracy. """
    name = roi.get_source(which).name

    frozen_sources = dict()
    for other_source in roi.get_sources():
        if np.any(other_source.model.free) and other_source.name != name:
            frozen_sources[other_source.name]=other_source.model.free.copy()
            roi.modify(which=other_source,free=False)

    roi.fit_extension(which=which, *args, **kwargs)
    for other_name,other_free in frozen_sources.items():
        roi.modify(which=other_name,free=other_free)

