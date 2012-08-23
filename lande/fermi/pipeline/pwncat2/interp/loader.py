import yaml
from os.path import expandvars, join, exists

class ResultsLoader(object):
    """ Class to load in results from analysis. """

    def __init__(self, pwndata,fitdir):
        self.pwndata = expandvars(pwndata)
        self.fitdir = expandvars(fitdir)

    def get_pwnlist(self):
        return yaml.load(open(self.pwndata)).keys()

    def get_results(self, pwn):
        all_results = ['results_%s_pointlike.yaml' % pwn, 
                       'results_%s_gtlike_at_pulsar.yaml' % pwn,
                       'results_%s_gtlike_point.yaml' % pwn,
                       'results_%s_gtlike_extended.yaml' % pwn]
        all_results = [join(self.fitdir,pwn,i) for i in all_results]

        for i in all_results: 
            if not exists(i):
                print '%s does not exist' % i
                return None
        g = [yaml.load(open(i)) for i in f]
        return merge_dict(*g)

    def get_sed(self,pwn,binning,hypothesis):
        filename=join(self.fitdir,pwn,'seds','sed_gtlike_%s_%s_%s.yaml' % (binning, hypothesis, pwn))
        return yaml.load(open(filename))

