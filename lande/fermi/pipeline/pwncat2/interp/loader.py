import yaml
from os.path import expandvars, join, exists

from uw.pulsar.phase_range import PhaseRange

from lande.utilities.tools import merge_dict
from lande.utilities.save import loaddict

class PWNResultsException(Exception): pass

class PWNResultsLoader(object):
    """ Class to load in results from analysis. """

    all_hypotheses = ['at_pulsar', 'point', 'extended']

    def __init__(self, pwndata, fitdir, phase_shift=None, verbosity=True):
        self.pwndata = expandvars(pwndata)
        self.fitdir = expandvars(fitdir)
        if phase_shift is not None:
            self.phase_shift = expandvars(phase_shift)
        else:
            print 'WARNING: PHASE SHIFT IS NONE'
            self.phase_shift = None

        self.verbosity = verbosity

    def get_pwnlist(self):
        return sorted(yaml.load(open(self.pwndata)).keys())

    def all_exists(self,pwn, get_seds=True, get_variability=True):
        try:
            if self.verbosity:
                print 'Checking if results for %s exist' % pwn
            self.get_results(pwn, require_all_exists=True, get_seds=get_seds, get_variability=get_variability, verbosity=False)
            if self.verbosity:
                print ' * Results for %s exists!' % pwn
            return True
        except PWNResultsException, ex:
            if self.verbosity:
                print ' * Results for %s do not exist:' % pwn,ex
            return False

    def get_results(self, pwn, require_all_exists=True, get_seds=True, get_variability=True, verbosity=None):
        filename = join(self.fitdir,pwn,'results_%s_general.yaml' % pwn)
        if verbosity or (verbosity is None and self.verbosity):
            print 'Getting results for %s' % pwn

        if not exists(filename): 
            raise PWNResultsException("%s does not exist" % filename)

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
                        raise PWNResultsException("%s does not exist" % filename)

        if get_variability:
            for hypothesis in ['at_pulsar','point']:
                filename =join(self.fitdir,pwn,'results_%s_variability_%s.yaml' % (pwn,hypothesis))
                if exists(filename):
                    results[hypothesis]['variability'] = loaddict(filename)
                else:
                    if require_all_exists:
                        raise PWNResultsException('%s does not exist' % filename)

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
                                raise PWNResultsException("%s does not exist" % filename)

        
        results['raw_phase'] = PhaseRange(results['phase'])

        del results['phase'] # to avoid ambiguity

        if self.phase_shift is not None:
            import pickle
            ps=pickle.load(open(self.phase_shift))
            shift = float(ps[pwn.replace('PSRJ','J')]['shift'])
            """ Convention taken from private correspondence with Matthew Kerr:
                    me: ok
                        so i am trying to implement phase shifts in my analysis
                        did you decide the right way to do it
                    Matthew: it seems to be +
                    me: so correct_phase = phase_in_file + shift
                    Matthew: yes """
            results['shifted_phase'] = results['raw_phase'].offset(shift)

        return results


    def get_sed(self,pwn,binning,hypothesis):
        sed=join(self.fitdir,pwn,'seds','sed_gtlike_%s_%s_%s.yaml' % (binning, hypothesis, pwn))
        return loaddict(sed)

