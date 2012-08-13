
from lande.fermi.likelihood.save import savedict

from snrhelp import pointlike_analysis, gtlike_analysis, pointlike_plots

def run(name, snrdata):

    roi = build_roi()

    results = dict()

    results['pointlike']=pointlike_analysis(roi,name)
    pointlike_plots(roi)
    results['gtlike']=gtlike_analysis(roi,name)

    savedict(results, 'results_%s.yaml' % name)
