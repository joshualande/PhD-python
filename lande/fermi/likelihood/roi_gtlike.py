""" Module for converting a pointlike analysis into a pyLikelihood analysis.

    Usage:

    # roi is a pointlike ROIAnalysis object
    roi=sa.roi_from_xml(...)

    from roi_gtlike import Gtlike

    # like is a pyLikelihood object
    gtlike = Gtlike(roi) 

    # the most relevant parameters are 
    gtlike = Gtlike(roi,
                  binsz       = 1/8.,  # spatial bins in gtlike
                  bigger_roi  = False, # bigger_roi=True => put pointlike circular ROI inside gtlike square ROI
    )

    # like is a BinnedAnalysis pyLikelihood object 
    # This can e used for pyLikelihood analysis
    like = gtlike.like

    # spectral fit with pyLikelihood
    like.fit()

    Author: Joshua Lande <joshualande@gmail.com>
"""
import os
from os.path import join
import shutil
import math
import numpy as np
from tempfile import mkdtemp

# pointlike
from uw.like.SpatialModels import RadiallySymmetricModel
from skymaps import SkyDir
from uw.utilities import keyword_options

# pylikelihood
from GtApp import GtApp
from BinnedAnalysis import BinnedObs,BinnedAnalysis
from UnbinnedAnalysis import UnbinnedObs, UnbinnedAnalysis
from pyLikelihood import ParameterVector

from lande.fermi.data import livetime

class Gtlike(object):

    common_defaults = (
        ("savedir",                 None, "Directory to put output files into. Default is to use a temporary file and delete it when done."),
        ("savedir_prefix",   '/scratch/', "Directory to put tempdir in."),
        ("optimizer",           "MINUIT", "Optimizer to use when fitting."),
        ("chatter",                    2, "Passed into the ScienceTools."),
    )

    defaults = common_defaults + (
            ('binsz',                   1/8., "bin size"),
            ("bigger_roi",             False, """If False, enscribe gtlike ROI in pointlike ROI. 
                                                 If True, enscribe pointlike ROI in gtlike ROI"""),
            ("proj",                   "ZEA", "Projection"),
            ("galactic",                True, "Coordinate system"),
            ("emin",                    None, "Minimum energy. Default is to get from ROI."),
            ("emax",                    None, "Maximum energy. Default is to get from ROI."),
            ("enumbins",                None, "Number of bins. Defualt is to get from ROI."),
            ("enable_edisp",           False, """ Enable energy dispersion. 
                                                  See https://confluence.slac.stanford.edu/display/ST/Energy+Dispersion+in+Binned+Likelihood"""),
            ("fix_pointlike_ltcube",   False, "Fix header of pointlike livetime cube so that it can be used by gtlike."),
            ("rfactor",                    2, "Goes into gtsrcmaps."),
            ("resample",               "yes", "Goes into gtsrcmaps."),
            ("minbinsz",                 0.1, "Goes into gtsrcmaps."),
            ("extended_dir_name",      None, "place to save converted extended sources."),
    )


    @staticmethod
    def make_evfile(roi,savedir):
        ft1files=roi.sa.pixeldata.ft1files
        if isinstance(ft1files,str):
            evfile=ft1files
        elif len(ft1files) == 1:
            evfile=ft1files[0]
        else:
            ft1list = join(savedir,'ft1files.lst')
            if not os.path.exists(ft1list):
                temp=open(ft1list,'w')
                temp.write('\n'.join(ft1files))
                temp.close()
            evfile='@'+ft1list
        return evfile


    @staticmethod
    def get_gtlike_irfs(roi):
        irfs=roi.sa.irf
        ct = roi.sa.pixeldata.conv_type
        if ct == 0 and 'FRONT' not in irfs: irfs += '::FRONT'
        if ct == 1 and 'BACK' not in irfs: irfs += '::BACK'
        return irfs

    @staticmethod
    def save_xml(roi, input_srcmdl_file, extended_dir_name):

        # shirnk all disk sources which are smaller than 0.05 degrees.
        # Otherwise, gtlike will crash. Anyway, no reason
        # to have extended source smaller than the binsize.
        shrink_list = [ source for source in roi.get_sources() if \
                       hasattr(source,'spatial_model') and \
                       isinstance(source.spatial_model,RadiallySymmetricModel) and \
                       source.spatial_model['sigma'] < 0.05  ]
        for src in shrink_list: src.spatial_model.shrink(size=0.05)

        roi.toXML(input_srcmdl_file,convert_extended=True,extended_dir_name=extended_dir_name,expand_env_vars=True)

        for src in shrink_list: src.spatial_model.unshrink()

    @staticmethod
    def get_ft2(roi):
        """ for now, only one ft2 file. """
        ft2files=roi.sa.pixeldata.ft2files
        if len(ft2files) > 1: raise Exception("Only 1 ft2 file at a time, for now")
        return ft2files[0]


    @keyword_options.decorate(defaults)
    def __init__(self, roi, **kwargs):
        """ Build a gtlike pyLikelihood object
            which is consistent with a pointlike roi. """
        keyword_options.process(self, kwargs)
        self.roi = roi

        if not roi.quiet: print 'Running a gtlike followup'

        self.old_dir=os.getcwd()
        if self.savedir is not None:
            self.savedata = True
            if not os.path.exists(self.savedir):
                os.makedirs(self.savedir)
        else:
            self.savedata = False
            self.savedir=mkdtemp(prefix=self.savedir_prefix)

        # put pfiles into savedir
        os.environ['PFILES']=self.savedir+';'+os.environ['PFILES'].split(';')[-1]

        if not roi.quiet: print 'Saving files to ',self.savedir

        if self.emin==None and self.emax==None and self.enumbins==None:
            self.emin,self.emax=roi.bin_edges[0],roi.bin_edges[-1]
            self.enumbins=len(roi.bin_edges)-1
        elif self.emin is not None and \
                self.emax is not None and \
                self.enumbins is not None:
            # all set
            pass
        else:
            raise Exception("emin, emax, and enumbins must all be set.")

        # Note that this formulation makes the gtlike slightly smaller than
        # the pointlike ROI (so the gtlike ROI is inside the pointlike ROI)
        roi_radius=np.degrees(max(_.radius_in_rad for _ in roi.bands))
        if self.bigger_roi:
            npix=int(math.ceil(2.0*roi_radius/self.binsz))
        else:
            npix=int(math.ceil(np.sqrt(2.0)*roi_radius/self.binsz))

        ct = roi.sa.pixeldata.conv_type

        cmap_file=join(self.savedir,'ccube.fits')
        srcmap_file=join(self.savedir,'srcmap.fits')
        bexpmap_file=join(self.savedir,'bexpmap.fits')
        input_srcmdl_file=join(self.savedir,'srcmdl.xml')
        cut_ft1=join(self.savedir,"ft1_cut.fits")


        ft2=Gtlike.get_ft2(roi)
        ltcube=roi.sa.pixeldata.ltcube 

        if self.fix_pointlike_ltcube:
            print 'Fixing pointlike ltcube %s' % ltcube
            livetime.fix_pointlike_ltcube(ltcube)
        
        irfs=Gtlike.get_gtlike_irfs(roi)

        if self.galactic:
            x,y,coordsys_str=roi.roi_dir.l(),roi.roi_dir.b(),'GAL'
        else:
            x,y,coordsys_str=roi.roi_dir.ra(),roi.roi_dir.dec(),'CEL'

        Gtlike.save_xml(roi, input_srcmdl_file, extended_dir_name=self.extended_dir_name)

        evfile=Gtlike.make_evfile(roi,self.savedir)

        if not os.path.exists(cut_ft1):
            if not roi.quiet: print 'Running gtselect'
            gtselect=GtApp('gtselect','dataSubselector')
            gtselect.run(infile=evfile,
                         outfile=cut_ft1,
                         ra=0, dec=0, rad=180,
                         tmin=0, tmax=0,
                         emin=self.emin, emax=self.emax,
                         zmax=180, convtype=ct,
                         chatter=self.chatter)
        else:
            if not roi.quiet: print '... Skiping gtselect'

        if not os.path.exists(cmap_file):
            if not roi.quiet: print 'Running gtbin (ccube)'
            gtbin=GtApp('gtbin','evtbin')
            gtbin.run(algorithm='ccube',
                      nxpix=npix, nypix=npix, binsz=self.binsz,
                      evfile=cut_ft1,
                      outfile=cmap_file,
                      scfile=ft2,
                      xref=x, yref=y, axisrot=0, proj=self.proj,
                      ebinalg='LOG', emin=self.emin, emax=self.emax, enumbins=self.enumbins,
                      coordsys=coordsys_str,
                      chatter=self.chatter)
        else:
            if not roi.quiet: print '... Skiping gtbin (ccube)'

        if not os.path.exists(bexpmap_file):
            # Use the default binning all sky, 1deg/pixel
            if not roi.quiet: print 'Running gtexpcube'
            gtexpcube=GtApp('gtexpcube2','Likelihood')
            gtexpcube.run(infile=ltcube,
                          cmap='none',
                          ebinalg='LOG', emin=self.emin, emax=self.emax, enumbins=self.enumbins,
                          outfile=bexpmap_file, proj='CAR',
                          nxpix=360, nypix=180, binsz=1,
                          irfs=irfs,
                          coordsys=coordsys_str,
                          chatter=self.chatter)
        else:
            if not roi.quiet: print '... Skiping gtexpcube'

        if not os.path.exists(srcmap_file):
            if not roi.quiet: print 'Running gtsrcmaps'
            gtsrcmaps=GtApp('gtsrcmaps','Likelihood')
            gtsrcmaps.run(scfile=ft2,
                          expcube=ltcube,
                          cmap=cmap_file,
                          srcmdl=input_srcmdl_file,
                          bexpmap=bexpmap_file,
                          outfile=srcmap_file,
                          irfs=irfs,
                          rfactor=self.rfactor,
                          resample=self.resample,
                          minbinsz=self.minbinsz,
                          chatter=self.chatter)
        else:
            if not roi.quiet: print '... Skiping gtsrcmaps'

        if not roi.quiet: print 'Creating Binned LIKE'
        obs=BinnedObs(srcMaps=srcmap_file,expCube=ltcube,binnedExpMap=bexpmap_file,irfs=irfs)

        self.like = BinnedAnalysis(binnedData=obs,srcModel=input_srcmdl_file,optimizer=self.optimizer)

        if self.enable_edisp:
            if not roi.quiet: print 'Enabeling energy dispersion'
            self.like.logLike.set_edisp_flag(True)


        if not roi.quiet: print 'Binned LIKE Created!'

    def __del__(self):
        if not self.savedata:
            if not self.roi.quiet: print 'Removing savedir',self.savedir
            shutil.rmtree(self.savedir)


class UnbinnedGtlike(object):

    defaults = Gtlike.common_defaults

    @keyword_options.decorate(defaults)
    def __init__(self, roi, **kwargs):
        """ Build a gtlike pyLikelihood object
            which is consistent with a pointlike roi. """
        keyword_options.process(self, kwargs)
        self.roi = roi

        if not roi.quiet: print 'Running a gtlike followup'

        self.old_dir=os.getcwd()
        if self.savedir is not None:
            self.savedata = True
            if not os.path.exists(self.savedir):
                os.makedirs(self.savedir)
        else:
            self.savedata = False
            self.savedir=mkdtemp(prefix=self.savedir_prefix)

        # put pfiles into savedir
        os.environ['PFILES']=self.savedir+';'+os.environ['PFILES'].split(';')[-1]

        if not roi.quiet: print 'Saving files to ',self.savedir

        cut_ft1=join(self.savedir,"ft1_cut.fits")
        input_srcmdl_file=join(self.savedir,'srcmdl.xml')
        expmap = join(self.savedir,"expmap.fits")

        ltcube=roi.sa.pixeldata.ltcube 
        ft2=Gtlike.get_ft2(roi)
        irfs=Gtlike.get_gtlike_irfs(roi)

        ct = roi.sa.pixeldata.conv_type
        radius = roi.sa.maxROI
        ra = roi.roi_dir.ra()
        dec = roi.roi_dir.dec()
        emin,emax=roi.bin_edges[0],roi.bin_edges[-1]

        Gtlike.save_xml(roi, input_srcmdl_file)

        evfile=Gtlike.make_evfile(roi,self.savedir)

        if not os.path.exists(cut_ft1):
            if not roi.quiet: print 'Running gtselect'
            gtselect=GtApp('gtselect','dataSubselector')
            gtselect.run(infile=evfile,
                         outfile=cut_ft1,
                         ra=ra, dec=dec, rad=radius,
                         tmin=0, tmax=0,
                         emin=emin, emax=emax,
                         zmax=180, convtype=ct,
                         chatter=self.chatter)
        else:
            if not roi.quiet: print '... Skiping gtselect'

        if not os.path.exists(expmap):
            # Run gtexpmap following suggestions from tutorial
            # pad 10deg on radius to account for nearby sources,
            # nlat has half degree pixels
            if not roi.quiet: print 'Running gtexpmap'
            gtexpmap=GtApp('gtexpmap')
            gtexpmap.run(evfile=cut_ft1,
                         scfile=ft2,
                         expcube=ltcube,
                         outfile=expmap,
                         irfs=irfs,
                         srcrad=radius+10,
                         nlong=int(np.ceil(0.5*(radius+10)*2)),
                         nlat=int(np.ceil(0.5*(radius+10)*2)),
                         nenergies=int(np.ceil(np.log10(emax)-np.log10(emin)))*4,
                         chatter=self.chatter,
                )
        else:
            if not roi.quiet: print '... Skiping gtexpmap'

        gtdiffrsp = GtApp('gtdiffrsp')
        if not roi.quiet: print 'Running gtdiffrsp'
        gtdiffrsp.run(evfile=cut_ft1,
                      scfile=ft2,
                      srcmdl=input_srcmdl_file,
                      irfs=irfs,
                      chatter=self.chatter,
                     )

        if not roi.quiet: print 'Creating Unbinned LIKE'
        obs = UnbinnedObs(eventFile=cut_ft1, scFile=ft2, expMap=expmap, expCube=ltcube, irfs=irfs)

        self.like = UnbinnedAnalysis(observation=obs,srcModel=input_srcmdl_file,optimizer=self.optimizer)

        if not roi.quiet: print 'Unbinned LIKE Created!'

    def __del__(self):
        if not self.savedata:
            if not self.roi.quiet: print 'Removing savedir',self.savedir
            shutil.rmtree(self.savedir)

