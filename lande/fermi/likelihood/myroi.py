""" lande_roi.py contains a ROI subclass which does things the way Joshua Lande likes.

    To create a landeROI object,

    import lande_roi
    sa = SpectralAnalysis(ds, ...)
    roi=lande_roi.LandeROI(sa.roi_from_xml(...))

    Author: Joshua Lande
"""
import numpy as np

from lande.utilities.tools import tolist

from lande_localize import *

# import stuff generally useful
from uw.like.Models import *
from uw.like.roi_analysis import *
from uw.like.roi_image import *
from uw.like.SpatialModels import *
from uw.like.roi_extended import *
from uw.like.roi_monte_carlo import *
from uw.like.roi_localize import *
from uw.like.roi_catalogs import *
from uw.like.roi_plotting import *
from uw.like.pointspec import *
from uw.like.pointspec_helpers import *
from uw.utilities.convolution import *
from uw.utilities.minuit import Minuit
from uw.utilities import keyword_options
from skymaps import *
from os.path import *
from roi_gtlike import *
from lande_plotting import *
from lande_extended import *
from lande_decorators import *

from . tools import galstr

from uw.like import sed_plotter
import pylab as P
import pylab
from pylab import *

from argparse import ArgumentParser
import warnings
import operator
import cPickle
import shutil
import os
import re
import math
import yaml
from tempfile import NamedTemporaryFile
import collections
from cStringIO import StringIO

import numpy
import numpy as np
from numpy import *

from scipy.optimize import fmin
from matplotlib.ticker import FuncFormatter

import pyfits 

from uw.like.pointspec_helpers import get_diffuse_source,FermiCatalog,PointSource,ExtendedSourceCatalog
from uw.like.roi_diffuse import DiffuseSource
from uw.like.roi_analysis import ROIAnalysis,decorate_with
from uw.like.pointspec import DataSpecification,SpectralAnalysis
from uw.like.roi_tsmap import TSCalc,TSCalcPySkyFunction
from uw.like.roi_extended import ExtendedSource,ROIExtendedModelAnalytic
from uw.like.SpatialModels import SpatialModel
from uw.like.Models import Model


from skymaps import SkyDir,SkyImage

import imp
from os.path import expandvars as e
from os.path import join as j

from uw.utilities import colormaps
from uw.utilities.colormaps import *


class Empty: pass


class LandeROI(ROIAnalysis):
    """ Subclass of ROIAnalysis with helper 
        functions that Joshua Lande finds useful. """

    # I dislike the pointlike defaults, so override them in LandeROI.
    localize=modify_defaults(update=True,verbose=True)(ROIAnalysis.localize)
    print_summary=modify_defaults(galactic=True, print_all_ts=True)(ROIAnalysis.print_summary)

    def remove_nan_background_sources(self,which):
        """ kind of a kludge, but nan sources cause problems with gtlike. """
        try:
            source=self.get_source(which)
        except:
            return

        for bg_source in self.get_sources():
            if source != bg_source:
                if np.isnan(source.model.i_flux()):
                    print 'Deleting background source %s because it is NaN' % source.name
                    self.del_source(source)

    @staticmethod 
    def get_point_source(name,catalog):
        return LandeROI.get_point_sources([name],catalog)

    @staticmethod 
    def get_point_sources(names,catalog):
        """ Input is a list of the name of 1FGL sources. What is
            return is a list of point source found from the names. """
        cat=FermiCatalog(catalog,free_radius=180) \
                if not isinstance(catalog,SourceCatalog) else catalog

        point_sources=[]

        for name in names:

            index=np.where(cat.names==name)[0]
            if len(index) < 1:
                raise Exception("Cannot find source %s in the catalog %s" % (name,catalog))
            if len(index) > 1:
                raise Exception("%s has too many counterpats in the catalog %s" % (name,catalog))
            index = int(index)

            point_sources.append(PointSource(cat.dirs[index],
                                 cat.names[index],cat.models[index],
                                 free_parameters=True))

        if len(names)==1:
            return point_sources[0]

        return point_sources

    @staticmethod 
    def get_background(*args):
        bg = []
        for source in args:
            if re.search(r'\.fit(s)?(\.gz)?$',source) is not None:
                bg.append(get_diffuse_source('MapCubeFunction',source,'PowerLaw',None,name=os.path.basename(source)))
            elif re.search(r'\.txt$',source) is not None:
                bg.append(get_diffuse_source('ConstantValue',None,'FileFunction',source,name=os.path.basename(source)))

            else:
                raise Exception("Diffuse Sources must end in .fit, .fits, .fit.gz, .fits.gz, or .txt (file is %s)" % os.path.basename(source))

        return bg[0] if len(args)==1 else bg


    def fit(self,*args,**kwargs):
        try:
            return super(LandeROI,self).fit(*args,**kwargs)
        except Exception, err:
            print '\n\n\n\nERROR FITTING: %s\n\n\n' % (str(err))

    def __init__(self,*args,**kwargs):
        """ Acts just like a ROIAnalysis object, except
            that 
            
            * if only one argument is passed 
                * If it is an ROIAnalysis object, fill in 
                  this object directly from that object. """
        if len(args)==1 and kwargs=={} and \
                isinstance(args[0],ROIAnalysis):
            self.__dict__.update(args[0].__dict__)
        else:
            super(LandeROI,self).__init__(*args,**kwargs)

    def cache_ft1(self,outfile=None):
        """ This function runs gtselect to replace
        """
        emin,emax=min(self.fit_emin),max(self.fit_emax)
        if outfile is None:
            # for some reason, the system pwd command works better then python os.getcwd()
            pwd=os.popen('pwd').readline().strip()
            outfile=os.path.join(pwd,'binned_emin_%g_emax_%g.fits' % (emin,emax))

        if os.path.exists(outfile):
            try:
                # Sometimes the outfile has problems with it
                # (the previous cache failed but still created
                # a file). So do some tests to make sure the
                # file looks good. If not delete and recreate it.

                # make warnigns raise exceptions
                warnings.filterwarnings('error')
                x=pyfits.open(outfile)

                # test to see if GTI & Photons are in file
                x['GTI']
                x['EVENTS']

                if np.all(x['EVENTS'].data.field('TIME')==0):
                    raise Exception("All TIMES in file are 0")
                if np.all(x['GTI'].data.field('START')==0):
                    raise Exception("All GTIs in START are 0")
                if np.all(x['GTI'].data.field('STOP')==0):
                    raise Exception("All GTIs in STOP are 0")

                if not x['EVENTS'].header.has_key('NDSKEYS'):
                    raise Exception("EVENTS header does not contain NDSKEYS entry")

            except Exception, ex:
                print 'Error reading file %s.' % outfile
                print 'Error: ',ex
                print 'Deleting and recaching file.'
                os.remove(outfile)

            warnings.resetwarnings()

        if not os.path.exists(outfile): 

            if not self.quiet: print 'Caching ft1 file. Saving to %s' % outfile

            if isinstance(self.sa.pixeldata.ft1files,collections.Iterable):
                temp=NamedTemporaryFile(delete=True)
                temp.write('\n'.join(self.sa.pixeldata.ft1files))
                temp.seek(0)
                infile='@%s' % temp.name
            else:
                infile=self.sa.pixeldata.ft1files

            import GtApp
            GtApp.GtApp("gtselect",'dataSubselector').run(
                    infile=infile,
                    outfile=outfile,
                    ra=self.roi_dir.ra(),
                    dec=self.roi_dir.dec(),
                    rad=self.sa.maxROI,
                    tmin=0, tmax=0,
                    emin=emin,
                    emax=emax,
                    zmax=180)

        self.sa.ft1files = outfile
        self.sa.ae.ft1files = outfile
        self.sa.pixeldata.ft1files = outfile


    def prune_empty_diffuse_models(self):
        """ Remove diffuse models which predict 0s. """
        del_sources = []
        for i,m in enumerate(self.bgm.diffuse_sources):
            if isinstance(m,ExtendedSource):
                continue
            if np.sum(band.bg_counts[i] for band in self.bands) < 1:
                if not self.quiet: print "Deleting source %s because it predicts < 1 pixel." % (m.name)
                del_sources.append(m)
        if len(del_sources)>0:
            self.del_ps(del_sources)
        else:
            print 'No diffuse sources to prune. All diffuse sources have positive pixels'


    def localize(self,*args,**kwargs):
        try:
            return super(LandeROI,self).localize(*args,**kwargs)
        except Exception, err:
            print '\n\n\n\nERROR LOCALIZING: %s\n\n\n' % (str(err))
            return None

    def fit_extension(self,*args,**kwargs):
        try:
            return super(LandeROI,self).fit_extension(*args,**kwargs)
        except Exception, err:
            print '\n\n\n\nERROR FITTING EXTENSION: %s\n\n\n' % (str(err))
            return None

    def gtlike_followup(self,which=None,output_srcmdl_file=None, **kwargs):

        gtlike=Gtlike(self,**kwargs)
        like=gtlike.like
        like.fit(covar=True)

        if output_srcmdl_file is not None:
            gtlike.like.writeXml(output_srcmdl_file)

        src_info = gtlike.get_gtlike_info_dict(which)
        return src_info
    
    def get_info_dict(self,which):
        """ Return a dictionary of everything Joshua finds
            interesting about the source. """

        src_info={}
        src_info['logLikelihood']={'spectral':float(-self.logLikelihood(self.parameters()))}

        try:
            manager,index=self.mapper(which)
        except:

            try:
                # see if it is 2 point hypothesis
                self.mapper('%s (first)' % which)
                self.mapper('%s (second)' % which)

            except:
                # If source which doesn't exist, return just the log likelihood
                return src_info

            results1=self.get_info_dict('%s (first)' % which)
            results1 = dict((k+'_first' if k!='logLikelihood' else k,v) for k,v in results1.items())

            results2=self.get_info_dict('%s (second)' % which)
            results2 = dict((k+'_second' if k!='logLikelihood' else k,v) for k,v in results2.items())

            results1.update(results2)
            return results1


        if manager == self.dsm:

            ds=self.dsm.diffuse_sources[index]
            if not isinstance(ds,ExtendedSource):
                raise Exception("get_info_dict can only be called on point/extended sources.")

            spatial=ds.spatial_model

            narr=spatial.param_names
            parr,earr=spatial.statistical(absolute=True,two_sided=False)
            parr,harr,larr=spatial.statistical(absolute=True,two_sided=True)
            for n,p,e,l,h in zip(narr,parr,earr,larr,harr):
                src_info['%s' % n]=[float(p),float(e)]
                src_info['%s_err_two_sided' % n]=[float(h),float(l)]

        # Save out TS information
        fit_dir = self.get_source(which).skydir
        src_info['fit_cel'] = [ float(fit_dir.ra()), float(fit_dir.dec()) ]
        src_info['fit_gal'] = [ float(fit_dir.l()),  float(fit_dir.b()) ]

        # Save out likelihood information
        src_info['logLikelihood']['bandfits']=float(self.bandFit(which=which))

        old_quiet=self.quiet; self.quiet=True
        src_info['TS']={'quick':float(self.TS(which=which,quick=True)),
                        'slow':float(self.TS(which=which,quick=False)),
                        'bandfits':float(self.TS(which=which,quick=False,bandfits=True))}
        self.quiet=old_quiet

        if not self.quiet: print 'pointlike TS =',src_info['TS']

        # Save out spectral information
        spectral = manager.models[index]

        for flux_name,emin,emax in [['Flux_100',100,100000],['Flux_1000',1000,100000],
                ['Flux',min(self.bin_edges),max(self.bin_edges)]]:

            flux,flux_err=spectral.i_flux(error=True,emin=emin,emax=emax)
            flux,flux_high,flux_low=spectral.i_flux(error=True,two_sided=True,emin=emin,emax=emax)

            src_info[flux_name]=[float(flux),float(flux_err)]
            src_info['%s_err_two_sided' % flux_name]=[float(flux_high),float(flux_low)]

        narr=spectral.param_names
        parr,earr=spectral.statistical(absolute=True,two_sided=False)
        parr,harr,larr=spectral.statistical(absolute=True,two_sided=True)
        for n,p,e,l,h in zip(narr,parr,earr,larr,harr):
            src_info[str(n)]=[float(p),float(e)]
            src_info['%s_err_two_sided' % n]=[float(h),float(l)]

        src_info['localization_error']=self.get_ellipse()

        return src_info

    def freeze_all_background_sources(self):
        background_sources = [_ for _ in self.dsm.diffuse_sources if not isinstance(_,ExtendedSource)]
        for bgs in background_sources:
            self.modify(which=bgs,free=False)

    def plot_profile(self,which,filename="profile.png", datafile="profile.yaml", **kwargs):
        p=ExtensionProfile(self, which)
        p.plot(filename)
        p.save(datafile)

    def extension_profile(self, which,filename='profile.yaml', **kwargs):
        p=ExtensionProfile(self, which)
        p.save(filename)

    def dual_localize(self,*args,**kwargs):
        return super(LandeROI,self).dual_localize(*args,**kwargs)

    def multi_localize(self,*args,**kwargs):
        dl = MultiLocalizer(self,*args,**kwargs)
        dl.localize()

    @decorate_with(ROIAnalysis.toRegion)
    def toRegion(self,filename='region.reg',**kwargs):
        return super(LandeROI,self).toRegion(filename=filename,**kwargs)

    @decorate_with(ROIAnalysis.toResults)
    def toResults(self,filename='results.reg',**kwargs):
        return super(LandeROI,self).toResults(filename=filename,**kwargs)

    @decorate_with(ROIAnalysis.plot_counts_map)
    def plot_counts_map(self,filename="counts_map.png",**kwargs):
        if not self.quiet: print 'Plotting counts map'
        return super(LandeROI,self).plot_counts_map(filename=filename,**kwargs)

    @decorate_with(ROISmoothedModel)
    def plot_model(self,filename="model_counts.png",**kwargs):
        if not self.quiet: print 'Plotting model counts'
        return super(LandeROI,self).plot_model(filename=filename,**kwargs)

    @decorate_with(ROIAnalysis.plot_counts_spectra)
    def plot_counts_spectra(self,filename="counts_spectra",**kwargs):
        if not self.quiet: print 'Plotting counts spectra'
        return super(LandeROI,self).plot_counts_spectra(filename=filename,**kwargs)

    @decorate_with(ROIAnalysis.plot_slice)
    def plot_slice(self,filename="slice.png",**kwargs):
        if not self.quiet: print 'Plotting slice'
        return super(LandeROI,self).plot_slice(filename=filename,**kwargs)

    @decorate_with(ROIAnalysis.plot_radial_integral)
    def plot_radial_integral(self,filename="radial_integral.png",**kwargs):
        if not self.quiet: print 'Plotting radial'
        return super(LandeROI,self).plot_radial_integral(filename=filename,**kwargs)

    @decorate_with(ROIAnalysis.plot_source)
    def plot_source(self,filename="source.png",**kwargs):
        if not self.quiet: print 'Plotting source'
        return super(LandeROI,self).plot_source(filename=filename,**kwargs)

    @decorate_with(ROIAnalysis.plot_sources)
    def plot_sources(self,filename="sources.png",**kwargs):
        if not self.quiet: print 'Plotting sources'
        return super(LandeROI,self).plot_sources(filename=filename,**kwargs)

    @decorate_with(ROIAnalysis.plot_tsmap)
    def plot_tsmap(self,filename="tsmap.png",**kwargs):
        if not self.quiet: print 'Plotting tsmap'
        return super(LandeROI,self).plot_tsmap(filename=filename,**kwargs)

    @decorate_with(ROISignificance)
    def plot_significance(self,filename="significance.png",**kwargs):
        if not self.quiet: print 'Plotting significance'
        ROISignificance(self,**kwargs).show(filename=filename)

    @decorate_with(ROIAnalysis.TS)
    @select_quiet
    def TS(self,*args,**kwargs):
        return super(LandeROI,self).TS(*args,**kwargs)

    @decorate_with(ROIAnalysis.plot_sed)
    @select_quiet
    def plot_sed(self,filename='sed.png',galmap=False,**kwargs):
        return super(LandeROI,self).plot_sed(filename=filename,galmap=galmap,**kwargs)

    @staticmethod
    def load(*args,**kwargs):
        """ like ROIAnalysis.load, but makes sure to return a LandeROI object. """
        return LandeROI(ROIAnalysis.load(*args,**kwargs))

    def save(self,filename='roi.dat',*args,**kwargs):
        super(LandeROI,self).save(filename,*args,**kwargs)

class VerboseROI(LandeROI):
    def __init__(self,*args,**kwargs):
        super(VerboseROI,self).__init__(*args,**kwargs)

        if not self.quiet:
            print 'irf = ',self.sa.irf
            print 'Phase Factor = ',self.phase_factor
            print 'Bin Edges: %s' % [int(round(i)) for i in self.bin_edges.tolist()]
            print '-> Total (not negative) logLikelihood = ',-self.logLikelihood(self.parameters())

    @select_quiet
    def plot_sed(self,**kwargs):
        if not self.quiet: print 'Plotting SED'
        super(VerboseROI,self).plot_sed(**kwargs)

    @select_quiet
    def TS(self,which,**kwargs):
        ts = super(VerboseROI,self).TS(which=which,**kwargs)
        name=self.get_source(which).name
        if not self.quiet: print 'Source %s has TS = %s' % (name,ts)
        return ts

    @select_quiet
    def TS_ext(self,which,**kwargs):
        ts_ext = super(VerboseROI,self).TS_ext(which=which,**kwargs)
        name=self.get_source(which).name
        if not self.quiet: print 'Source %s has TS_ext = %s' % (name,ts_ext)
        return ts_ext

    def print_summary(self,*args,**kwargs):
        """ This function needs to shut up TS. """
        old_quiet=self.quiet
        self.quiet=True
        super(VerboseROI,self).print_summary(*args,**kwargs)
        print '-> Total (not negative) logLikelihood = ',-self.logLikelihood(self.parameters())
        self.quiet=old_quiet


    @select_quiet
    def localize(self,which,*args,**kwargs):
        """ My version of localize just prints a bunch of 
            stuff before & after the localization. """

        source = self.get_source(which)
        initdir=source.skydir
        name=source.name

        if not self.quiet:
            print;print
            print 'Beginning Localization of source %s' % name
            print '  Initial Position is ',galstr(initdir)

        retval=super(VerboseROI,self).localize(which,*args,**kwargs)

        if retval is not None and not self.quiet:
            print 'Done Localizing source %s' % name
            print '  Initial Position is:',galstr(initdir)
            print '  Final Position is:',galstr(retval[0])
            print '  Distance (deg):',retval[2]
            print '  lsigma (deg):',self.lsigma
            print '  Distance/lsigma:',retval[2]/self.lsigma
            print '  TS:',self.TS(which=source,quick=True,quiet=True)
            print '  Total (not negative) logLikelihood = ',-self.logLikelihood(self.parameters())
            print;print

        return retval

    @select_quiet
    def fit_extension(self,which,*args,**kwargs):

        source = self.get_source(which)
        name=source.name
        sm = source.spatial_model

        if not self.quiet:
            print;print
            print 'Beginning Localization of source %s' % name
            print '  Initial Spatial is ',sm.pretty_string()

            retval=super(VerboseROI,self).fit_extension(which,*args,**kwargs)
            print 'Done Localizing source %s' % name
            print '  Final Spatial is:',sm.pretty_string()
            print '  Final TS (quick=True) is:',self.TS(which=which,quick=True,quiet=True)
            if sm.can_shrink():
                print '  Final TS_ext (refit=False) is:',self.TS_ext(which=which,refit=False,quiet=True)
            print '  Total (not negative) logLikelihood = ',-self.logLikelihood(self.parameters())

        return retval

    @decorate_with(ROIAnalysis.fit)
    def fit(self,*args,**kwargs):
        if not self.quiet: print '\n\nFitting ROI'
        retval=super(VerboseROI,self).fit(*args,**kwargs)
        if not self.quiet: 
            print '  Total (not negative) logLikelihood = ',-self.logLikelihood(self.parameters())
            print '\n\n'
        return retval

    def __str__(self):
        """ This function prints out the ROI in a format that josh likes
            
            * it is compatable yaml
            * merges display of point + extended sources,
            * printing the direction of point sources
            * prints the (quick=True) TS of all components. """
        r=[]
        r.append('POINT+EXTENDED SOURCE FITS:')

        def format_ps(source):
            return '\n'.join(['  %s:'%source.name,
                              '    %s:'%(source.model.full_name()),
                              '      '+source.model.__str__(indent='      '),
                              '    dir: %s' % galstr(source.skydir),
                              '    dist: %s' % np.degrees(self.roi_dir.difference(source.skydir)),
                              '    TS: %.3f ' % (self.TS(which=source,quick=True,quiet=True))])

        def format_es(source):
            return '\n'.join(['  %s:'%source.name,
                              '    %s:'%(source.model.full_name()),
                              '      '+source.model.__str__(indent='      '),
                              '    %s:'%(source.spatial_model.full_name()),
                              '      '+source.spatial_model.__str__(indent='      '),
                              '    dist: %s' % np.degrees(self.roi_dir.difference(source.spatial_model.center)),
                              '    TS: %.3f ' % (self.TS(which=source,quick=True,quiet=True))])

        for source in self.get_sources():
            if np.any(source.model.free):
                r.append(format_ps(source) \
                         if isinstance(source,PointSource) 
                         else format_es(source))

        r.append('BACKGROUND SOURCE FITS:')

        def format_bg(source):
            return '\n'.join(['  %s:'%source.name,
                              '    %s:'%(source.smodel.full_name()),
                              '      '+source.smodel.__str__(indent='      '),
                              '    TS: %.3f ' % (self.TS(which=source,quick=True,quiet=True))])

        for source in self.bgm.diffuse_sources: 
            if not isinstance(source,ExtendedSource):
                r.append(format_bg(source))

        if (self.logl is not None) and (self.prev_logl is not None):
            r.append('Log likelihood change: %.2f'%(self.logl - self.prev_logl))

        return '\n\n'+'\n\n'.join(r)+'\n\n'

    @staticmethod
    def load(filename,*args,**kwargs):
        """ like ROIAnalysis.load, but makes sure to return a LandeROI object. """
        print 'Loading ROI from %s' % filename
        roi=VerboseROI(ROIAnalysis.load(filename,*args,**kwargs))
        print 'Loaded ROI:'
        roi.print_summary()
        return roi

    def save(self,filename='roi.dat',*args,**kwargs):
        print 'Saving ROI to %s' % filename
        self.print_summary()
        print 'Saved ROI:'
        super(VerboseROI,self).save(filename,*args,**kwargs)

load=VerboseROI.load




