import yaml
from os.path import expandvars, join, exists

from lande.utilities.tools import merge_dict
from lande.utilities.save import loaddict

class PWNResultsLoader(object):
    """ Class to load in results from analysis. """

    def __init__(self, pwndata, fitdir):
        self.pwndata = expandvars(pwndata)
        self.fitdir = expandvars(fitdir)

    def get_pwnlist(self):
        return sorted(yaml.load(open(self.pwndata)).keys())

    def get_results(self, pwn, require_all_exists=True, load_seds=False):
        all_results = ['results_%s_pointlike.yaml' % pwn, 
                       'results_%s_gtlike_at_pulsar.yaml' % pwn,
                       'results_%s_gtlike_point.yaml' % pwn,
                       'results_%s_gtlike_extended.yaml' % pwn,
                       'results_%s_variability_at_pulsar.yaml' % pwn,
                       'results_%s_variability_point.yaml' % pwn]
        all_results = [join(self.fitdir,pwn,i) for i in all_results]

        if not exists(all_results[0]):
            print '%s does not exist' % all_results[0]
            return None

        if require_all_exists:
            for i in all_results: 
                if not exists(i):
                    print '%s does not exist' % i
                    return None
        else:
            g = [loaddict(i) for i in all_results if exists(i)]
            results=merge_dict(*g)

        if load_seds:
            for h in ['at_pulsar', 'point', 'extended']:
                for binning in ['1bpd','2bpd','4bpd']:
                    sed = join(self.fitdir,pwn,'seds','sed_gtlike_%s_%s_%s.yaml' % (binning,h,pwn))
                    if require_all_exists and not exists(sed):
                        raise Exception("%s does not exist" % sed)
                    results[h]['gtlike'][binning] = loaddict(sed)
        return results


    def get_sed(self,pwn,binning,hypothesis):
        filename=join(self.fitdir,pwn,'seds','sed_gtlike_%s_%s_%s.yaml' % (binning, hypothesis, pwn))
        return yaml.load(open(filename))

