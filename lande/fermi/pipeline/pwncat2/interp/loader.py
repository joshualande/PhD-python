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

    def get_results(self, pwn, require_all_exists=True, get_seds=True, get_bandfits=True, get_variability=True):
        all_results = ['results_%s_pointlike.yaml' % pwn, 
                       'results_%s_gtlike_at_pulsar.yaml' % pwn,
                       'results_%s_gtlike_point.yaml' % pwn,
                       'results_%s_gtlike_extended.yaml' % pwn]
        if get_variability:
            all_results += ['results_%s_variability_at_pulsar.yaml' % pwn,
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

        g = [loaddict(i) for i in all_results if exists(i)]
        results=merge_dict(*g)

        for hypothesis in ['at_pulsar', 'point', 'extended']:

            if get_seds:
                for code,all_binning in [['gtlike',['1bpd','2bpd','4bpd']], ['pointlike',['4bpd']]]:
                    for binning in all_binning:
                        sed = join(self.fitdir,pwn,'seds','sed_%s_%s_%s_%s.yaml' % (code,binning,hypothesis,pwn))
                        if require_all_exists and not exists(sed):
                            raise Exception("%s does not exist" % sed)
                        if exists(sed):
                            if not results[hypothesis][code].has_key('seds'):
                                results[hypothesis][code]['seds']=dict()
                            results[hypothesis][code]['seds'][binning] = loaddict(sed)
            
            if get_bandfits:
                for code in ['gtlike']:
                    bandfit = join(self.fitdir,pwn,'data','bandfit_%s_%s_%s.yaml' % (code, hypothesis, pwn))
                    if require_all_exists and not exists(bandfit):
                        raise Exception("%s does not exist" % bandfit)
                    if exists(bandfit):
                        results[hypothesis][code]['bandfit'] = loaddict(bandfit)

        return results


    def get_sed(self,pwn,binning,hypothesis):
        sed=join(self.fitdir,pwn,'seds','sed_gtlike_%s_%s_%s.yaml' % (binning, hypothesis, pwn))
        return loaddict(sed)

