
from lande.utilities.save import savedict
from lande.fermi.likelihood.tools import force_gradient

from . setup import build_roi
from . help import pointlike_analysis, gtlike_analysis, pointlike_plots

def run(name, snrdata, latdata):

    force_gradient(use_gradient=False)

    roi = build_roi(name, snrdata, latdata)

    results = dict()

    kwargs = dict(plotdir='plotdir')

    results['pointlike']=pointlike_analysis(roi,name, **kwargs)
    pointlike_plots(roi)
    results['gtlike']=gtlike_analysis(roi,name, **kwargs)

    savedict(results, 'results_%s.yaml' % name)
