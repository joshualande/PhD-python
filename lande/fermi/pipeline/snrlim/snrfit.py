from snrhelp import pointlike_analysis, gtlike_analysis, pointlike_plots

from lande.fermi.likelihood.save import savedict

parser = ArgumentParser()
parser.add_argument("--name", required=True, help="Name of the pulsar")
parser.add_argument("--snrdata", required=True)
args=parser.parse_args()

roi = build_roi()

results = dict()

results['pointlike']=pointlike_analysis(roi,name)
pointlike_plots(roi)
results['gtlike']=gtlike_analysis(roi,name)

savedict(results, 'results_%s.yaml' % name)
