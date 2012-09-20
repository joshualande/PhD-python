import yaml
from os.path import expandvars, join, exists

from lande.utilities.tools import merge_dict
from lande.utilities.save import loaddict

class PWNResultsLoader(object):
    """ Class to load in results from analysis. """

    all_hypotheses = ['at_pulsar', 'point', 'extended']

    def __init__(self, pwndata, fitdir):
        self.pwndata = expandvars(pwndata)
        self.fitdir = expandvars(fitdir)

    def get_pwnlist(self):
        return sorted(yaml.load(open(self.pwndata)).keys())

    def get_results(self, pwn, require_all_exists=True, get_seds=True, get_bandfits=True, get_variability=True):
        filename = join(self.fitdir,pwn,'results_%s_general.yaml' % pwn)
        print filename,exists(filename)
        if not exists(filename): return None

        results = loaddict(filename)
        for hypothesis in self.all_hypotheses:
            results[hypothesis] = dict()

        for code in ['gtlike','pointlike']:
            for hypothesis in self.all_hypotheses:
                filename=join(self.fitdir,pwn,'results_%s_%s_%s.yaml' % (pwn,code,hypothesis))
                if exists(filename):
                    results[hypothesis][code] = loaddict(filename)
                else:
                    if require_all_exists:
                        raise Exception("%s does not exist" % filename)

        if get_variability:
            for hypothesis in ['at_pulsar','point']:
                filename =join(self.fitdir,pwn,'results_%s_variability_%s.yaml' % (pwn,hypothesis))
                if exists(filename):
                    results[hypothesis]['variability'] = loaddict(filename)
                else:
                    if require_all_exists:
                        raise Exception('%s does not exist' % filename)

        if get_seds:
            for hypothesis in self.all_hypotheses:
                for code,all_binning in [['gtlike',['1bpd','2bpd','4bpd']], ['pointlike',['4bpd']]]:
                    results[hypothesis][code]['seds'] = dict()
                    for binning in all_binning:
                        filename = join(self.fitdir,pwn,'seds','sed_%s_%s_%s_%s.yaml' % (code,binning,hypothesis,pwn))
                        if exists(filename):
                            results[hypothesis][code]['seds'][binning] = loaddict(filename)
                        else:
                            if require_all_exists:
                                raise Exception("%s does not exist" % filename)
            
        if get_bandfits:
            for hypothesis in self.all_hypotheses:
                    bandfit = join(self.fitdir,pwn,'data','bandfit_%s_%s_%s.yaml' % ('gtlike', hypothesis, pwn))
                    if exists(bandfit):
                        results[hypothesis][code]['bandfit'] = loaddict(bandfit)
                    else:
                        if require_all_exists:
                            raise Exception("%s does not exist" % bandfit)
        return results


    def get_sed(self,pwn,binning,hypothesis):
        sed=join(self.fitdir,pwn,'seds','sed_gtlike_%s_%s_%s.yaml' % (binning, hypothesis, pwn))
        return loaddict(sed)

