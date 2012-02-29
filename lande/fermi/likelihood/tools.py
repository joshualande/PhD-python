""" This file contains various function which I have found useful. """

from uw.like.roi_analysis import ROIAnalysis
from uw.like.Models import Model

from SED import SED


def gtlike_or_pointlike(f_gtlike, f_pointlike, like_or_roi, *args, **kwargs):
    """ Note, like_or_roi can be either an ROI or a spectral model. """

    from pyLikelihood import Function
    from BinnedAnalysis import BinnedAnalysis

    if isinstance(like_or_roi, BinnedAnalysis) or \
       isinstance(like_or_roi,Function):
        f=f_gtlike
    elif isinstance(like_or_roi, ROIAnalysis) or \
       isinstance(like_or_roi,Model):
        f=f_pointlike
    else:
        raise Exception("like_or_roi must be of type BinnedAnalysis or ROIAnalysis")
    return f(like_or_roi, *args, **kwargs)


def force_gradient(use_gradient):
    """ A kludge to force use_gradient everywhere! """
    from uw.like.roi_analysis import ROIAnalysis
    from lande.utilities.decorators import modify_defaults
    ROIAnalysis.fit=modify_defaults(use_gradient=use_gradient)(ROIAnalysis.fit)

def galstr(skydir):
    return 'SkyDir(%.3f,%.3f,SkyDir.GALACTIC)' % (skydir.l(),skydir.b())


