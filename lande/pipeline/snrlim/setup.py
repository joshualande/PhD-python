import numpy as np
import math

from skymaps import SkyDir
from uw.like.pointspec import DataSpecification, SpectralAnalysis
from uw.like.pointspec_helpers import get_default_diffuse
from uw.like.roi_catalogs import Catalog2FGL
from uw.like.roi_extended import ExtendedSource
from uw.like.SpatialModels import Disk
from uw.like.Models import PowerLaw
 
from lande.utilities.save import loaddict

from . import data

def build_roi(name, snrdata, latdata):

    snrdata=loaddict(snrdata)
    latdata=loaddict(latdata)

    roi_dir = SkyDir(*snrdata[name]['cel'])
    snrsize = snrdata[name]['size']

    if isinstance(snrsize,list) and len(snrsize) == 2:
        snrradius = math.sqrt(snrsize[0]*snrsize[1])/2.0
    else:
        snrradius = snrsize/2.0

    ds = DataSpecification(**latdata['data'])

    sa = SpectralAnalysis(ds,
                          binsperdec = 4,
                          emin       = 1e4,
                          emax       = 10**5.5,
                          irf        = "P7SOURCE_V6",
                          roi_dir    = roi_dir,
                          maxROI     = 10,
                          minROI     = 10,
                          event_class= 0)


    diffuse_sources = get_default_diffuse(**latdata['diffuse'])

    catalog = Catalog2FGL(**latdata['catalog'])

    roi=sa.roi(point_sources=[],
               diffuse_sources=diffuse_sources,
               catalogs=catalog)

    print 'bins',roi.bin_edges

    for source in roi.get_sources():
        if np.degrees(source.skydir.difference(roi_dir)) < snrradius + 0.5:
            roi.del_source(source)

    snr = ExtendedSource(
        name = name,
        model = PowerLaw(),
        spatial_model = Disk(sigma=snrradius, center=roi_dir)
    )

    roi.add_source(snr)

    return roi
