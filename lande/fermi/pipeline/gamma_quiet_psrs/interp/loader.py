from os.path import expandvars, join, exists
from lande.utilities.save import loaddict

class ResultsException(Exception): pass

class ResultsLoader(object):
    """ Class to load in results from analysis. """

    all_hypotheses = ['at_pulsar', 'point', 'extended']

    def __init__(self, gamma_quiet_psrs_data, analysisdir, verbosity=True):
        self.gamma_quiet_psrs_data = loaddict(expandvars(gamma_quiet_psrs_data))
        self.analysisdir = expandvars(analysisdir)
        self.verbosity = verbosity

    def get_psr_list(self):
        return self.gamma_quiet_psrs_data.keys()


    def get_results(self, psr, require_all_exists=True, get_seds=True, get_variability=True, verbosity=None):
        def load_me_maybe(filename):
            if require_all_exists and not exists(filename): 
                raise ResultsException("%s does not exist" % filename)
            elif not exists(filename):
                return dict()
            return loaddict(filename)

        results = dict(
            params=load_me_maybe(join(self.analysisdir,psr,'results_%s_general.yaml' % psr))
        )
        for hypothesis in self.all_hypotheses:
            results[hypothesis]=dict()
            for code in ['gtlike','pointlike']:
                results[hypothesis][code]=load_me_maybe(join(self.analysisdir,psr,'results_%s_%s_%s.yaml' % (psr,code,hypothesis)))
        return results
