""" This is mostly Romain's code. """
import os
import readXml
# pointlike
from skymaps import SkyDir
from uw.utilities import keyword_options

# pylikelihood
from GtApp import GtApp

from datetime import date

IRF="P7SOURCE_V6"

def cmap(roi,like,ft1,ft2,cmap=None,ccube=None,npix=100,binsz=0.1,proj="ZEA"):
    if like.coordsystem==SkyDir.GALACTIC:
        x,y,coordsys_str=roi.roi_dir.l(),roi.roi_dir.b(),'GAL'
    else:
        x,y,coordsys_str=roi.roi_dir.ra(),roi.roi_dir.dec(),'CEL'
    if cmap is not None : 
        gtbin=GtApp('gtbin','evtbin')
        gtbin.run(algorithm='cmap',
                  nxpix=npix, nypix=npix, binsz=binsz,
                  evfile=ft1,
                  outfile=cmap,
                  scfile=ft2,
                  xref=x, yref=y, axisrot=0, proj=proj,
                  coordsys=coordsys_str)
        
    if ccube is not None:
        gtbin=GtApp('gtbin','evtbin')
        gtbin.run(algorithm='ccube',
                  nxpix=npix, nypix=npix, binsz=like.binsz,
                  evfile=ft1,
                  outfile=ccube,
                  scfile=ft2,
                  xref=x, yref=y, axisrot=0, proj=proj,
                  ebinalg='LOG', emin=like.emin, emax=like.emax, enumbins=like.enumbins,
                  coordsys=coordsys_str)

def modelmap(like,ft2,cmap,outfile):
    """Fonction to generate a model map.
    Use :
    like : a pylike object (BinnedAnalysis)
    ft2 : a spacecraft file
    cmap : a 3D cmap
    outfile : the name of the output"""

    #taking the needed data
    expcube=like.binnedData.expCube
    bexpmap=like.binnedData.binnedExpMap

    like.srcModel = "srcmodel_%s.xml"%date.today()
    like.writeXml()
    srcModel = readXml.SourceModel("srcmodel_%s.xml"%date.today())


    #creating a srcmap
    gtsrcmaps = GtApp('gtsrcmaps')
    gtsrcmaps['scfile'] = ft2
    gtsrcmaps['expcube'] = expcube
    gtsrcmaps['bexpmap'] = bexpmap
    gtsrcmaps['irfs'] =IRF
    gtsrcmaps['cmap'] = cmap
    gtsrcmaps['srcmdl'] = "srcmodel_%s.xml"%date.today()
    gtsrcmaps['outfile'] = "srcm_%s.fits"%date.today()
    gtsrcmaps['ptsrc'] = 'no'
    gtsrcmaps['debug'] = 'yes'
    
    gtsrcmaps.run()                                                                                         

    #creating the model map
    gtmodel= GtApp('gtmodel')
    gtmodel['srcmaps'] = "srcm_%s.fits"%date.today()
    gtmodel['srcmdl'] = "srcmodel_%s.xml"%date.today()
    gtmodel['outfile']=outfile
    gtmodel['irfs']=IRF
    gtmodel['expcube']=expcube
    gtmodel['bexpmap']=bexpmap
    gtmodel.run()

    os.system("rm -rf srcmodel_%s.xml"%date.today())
    os.system("rm -rf srcm_%s.fits"%date.today())


def excess_and_residual(roi,glike,ft1,ft2,outdir,name):
    #coordsys=glike.coordsystem
    #like=glike.like

    
    cmap3D="%s/CMAP_3D_%s.fits"%(outdir,name)
    #generating counts map
    cmap(roi,glike,ft1[0],ft2[0],cmap="%s/CMAP_2D_%s.fits"%(outdir,name),ccube=cmap3D,npix=200,binsz=0.1,proj="ZEA")
    like=glike.like
    
    #saving old data 
    like.srcModel = "srcmodel_save_%s.xml"%date.today()
    like.writeXml()
    
    
    #Create a model counts map including the source
    modelmap(like,ft2,cmap3D,"%s/model_countsmap_with_%s.fits"%(outdir,name))

    #Delete the source
    like.deleteSource(name)
    modelmap(like,ft2,cmap3D,"%s/model_countsmap_without_%s.fits"%(outdir,name))

    #Reload initial roi

    srcModel = readXml.SourceModel("srcmodel_%s.xml"%date.today())

    
    os.system("rm -rf srcmodel_%s.xml"%date.today())
