""" Code to clasify regions. """
import copy
import yaml
from os.path import expandvars

from . loader import PWNResultsLoader


def compare_classifications(pwndata,fitdir,pwn_classification):

    loader = PWNResultsLoader(pwndata, fitdir, verbosity=False)

    auto=PWNAutomaticClassifier(loader)
    manual=PWNManualClassification(pwn_classification)

    pwnlist = loader.get_pwnlist()

    for pwn in pwnlist:
        manual_results = manual.get_classification(pwn)

        if manual_results['spatial_model'] is None or \
               manual_results['spectral_model'] is None or \
               manual_results['source_class'] is None:
            print "%s - no manual classification" % pwn
        elif not loader.all_exists(pwn, get_variability=False):
            print "%s - results don't exist" % pwn
        else:
            auto_results = auto.get_classification(pwn)

            if auto_results['spatial_model'] == manual_results['spatial_model'] and \
                    auto_results['spectral_model'] == manual_results['spectral_model'] and \
                    auto_results['source_class'] == manual_results['source_class']:
                print '%s - classifications agree' % pwn
            else:
                print '%s - classifications do not agree' % pwn
                print '                    %20s %20s' % ('automatic','manual')
                print '  - spectral model: %20s %20s' % (auto_results['spectral_model'],manual_results['spectral_model'])
                print '  - spatial model:  %20s %20s' % (auto_results['spatial_model'],manual_results['spatial_model'])
                print '  - source class:   %20s %20s' % (auto_results['source_class'],manual_results['source_class'])





class PWNClassifier(object):

    allowed_source_class = ['Pulsar', 'PWN', 'Confused', 'Upper_Limit']
    allowed_spatial_models = ['At_Pulsar','Point','Extended']
    allowed_spectral_models = ['FileFunction','PowerLaw','PLSuperExpCutoff']


    def get_results(self, pwn):

        classifier = self.get_classification(pwn)

        spatial_model=classifier['spatial_model']
        spectral_model=classifier['spectral_model']
        source_class = classifier['source_class']

        if spatial_model is None or spectral_model is None or source_class is None:
            print '%s has not been classified yet, skipping' % pwn
            return None

        assert spatial_model in PWNClassifier.allowed_spatial_models
        assert spectral_model in PWNClassifier.allowed_spectral_models
        assert source_class in PWNClassifier.allowed_source_class

        results = self.loader.get_results(pwn, require_all_exists=True, get_variability=False)

        if results is None:
            print 'Results for %s is not done yet, skipping' % pwn
            return None

        point_gtlike = results['point']['gtlike']
        extended_gtlike = results['extended']['gtlike']

        gtlike = results[spatial_model.lower()]['gtlike']
        pointlike = results[spatial_model.lower()]['pointlike']

        point_cutoff=results['point']['gtlike']['test_cutoff']

        d = copy.copy(classifier)

        # likelihood stuff

        d['ts_point'] = max(point_gtlike['TS']['reoptimize'],0)

        if source_class in ['Confused', 'Pulsar', 'PWN']:
            d['ts_ext'] = max(extended_gtlike['TS']['reoptimize']-point_gtlike['TS']['reoptimize'],0)
            d['ts_cutoff'] = max(point_cutoff['TS_cutoff'],0)

        d['ts_var'] = None

        # spectral stuff


        if source_class != 'Upper_Limit':

            if spectral_model in ['PowerLaw','FileFunction']:
                d['flux'] = gtlike['flux']['flux']
                d['flux_err'] = gtlike['flux']['flux_err']

                if spectral_model == 'PowerLaw':
                    d['prefactor'] = gtlike['spectrum']['Prefactor']
                    d['prefactor_err'] = gtlike['spectrum']['Prefactor_err']

                    d['index'] = gtlike['spectrum']['Index']
                    d['index_err'] = gtlike['spectrum']['Index_err']

                    d['model_scale'] = gtlike['spectrum']['Scale']

                elif spectral_model == 'FileFunction':

                    d['prefactor'] = gtlike['spectrum']['Normalization']
                    d['prefactor_err'] = gtlike['spectrum']['Normalization']

                    d['index'] = None
                    d['index_err'] = None

                d['cutoff'] = None
                d['cutoff_err'] = None

            elif spectral_model == 'PLSuperExpCutoff':
                d['flux'] = gtlike['test_cutoff']['hypothesis_1']['flux']['flux']
                d['flux_err'] = gtlike['test_cutoff']['hypothesis_1']['flux']['flux_err']
            
                d['prefactor'] = gtlike['test_cutoff']['hypothesis_1']['spectrum']['Prefactor']
                d['prefactor_err'] = gtlike['test_cutoff']['hypothesis_1']['spectrum']['Prefactor']

                d['index'] = gtlike['test_cutoff']['hypothesis_1']['spectrum']['Index1']
                d['index_err'] = gtlike['test_cutoff']['hypothesis_1']['spectrum']['Index1_err']

                d['model_scale'] = gtlike['test_cutoff']['hypothesis_1']['spectrum']['Scale']

                d['cutoff'] = gtlike['test_cutoff']['hypothesis_1']['spectrum']['Cutoff']
                d['cutoff_err'] = gtlike['test_cutoff']['hypothesis_1']['spectrum']['Cutoff_err']

        # spatial stuff

        d['ra'] = pointlike['position']['equ'][0]
        d['dec'] = pointlike['position']['equ'][0]

        d['glon'] = pointlike['position']['gal'][0]
        d['glat'] = pointlike['position']['gal'][0]

        if spatial_model in [ 'Point', 'Extended' ]: 

            ellipse = pointlike['spatial_model']['ellipse']
            if ellipse.has_key('lsigma'):
                d['poserr'] = ellipse['lsigma']
            else:
                print 'WARNING: localization failed for %s' % pwn
                d['poserr'] = None

        if spatial_model == 'Extended':
            d['extension'] = pointlike['spatial_model']['Sigma']
            d['extension_err'] = pointlike['spatial_model']['Sigma_err']

        return d


class PWNAutomaticClassifier(PWNClassifier):

    def __init__(self, loader):
        self.loader = loader

    def get_classification(self, pwn):
        results = self.loader.get_results(pwn, require_all_exists=True, get_variability=False)

        point_gtlike = results['point']['gtlike']
        extended_gtlike = results['extended']['gtlike']
        cutoff=point_gtlike['test_cutoff']

        ts_point = max(point_gtlike['TS']['reoptimize'],0)
        ts_ext = max(extended_gtlike['TS']['reoptimize']-point_gtlike['TS']['reoptimize'],0)

        assert point_gtlike['spectrum']['name'] == extended_gtlike['spectrum']['name']
        spectral_name = point_gtlike['spectrum']['name']

        try:
            ts_cutoff = max(cutoff['TS_cutoff'],0)
        except:
            ts_cutoff = None

        if ts_point > 25:

            if ts_cutoff > 16:
                source_class = 'Pulsar'
                spectral_model = 'PLSuperExpCutoff'
            else:
                source_class = 'Confused'
                spectral_model = spectral_name

            if ts_ext > 16:
                spatial_model = 'Extended'
            else:
                spatial_model = 'Point'

        else:
            source_class = 'Upper_Limit'
            spatial_model = 'At_Pulsar'
            spectral_model = spectral_name

        return dict(source_class=source_class, 
                    spatial_model=spatial_model,
                    spectral_model=spectral_model)

    @staticmethod
    def get_automatic_classify(pwndata, fitdir):
        loader = PWNResultsLoader(pwndata, fitdir)
        pwnlist = loader.get_pwnlist()

        classifier = PWNAutomaticClassifier(loader=loader)

        return {pwn:classifier.get_classification(pwn) for pwn in pwnlist if loader.all_exists(pwn, get_variability=False)}


class PWNManualClassification(object):
    """ Classify a PWN """


    def __init__(self, pwn_classification):
        self.pwn_classifications = yaml.load(open(expandvars(pwn_classification)))

    def get_classification(self, pwn):
        return self.pwn_classifications[pwn]

    @staticmethod
    def get_manual_classify(pwndata, fitdir):
        """ Create an empty dict which can be used for manual classification. """
        loader = PWNResultsLoader(pwndata, fitdir)
        pwnlist = loader.get_pwnlist()
        return {pwn:dict(source_class=None, 
                        spatial_model=None,
                        spectral_model=None) for pwn in pwnlist}
