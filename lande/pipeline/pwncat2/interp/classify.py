""" Code to clasify regions. """

from . loader import PWNResultsLoader

def auto_classify(pwndata, fitdir):
    loader = PWNResultsLoader(pwndata, fitdir)
    pwnlist = loader.get_pwnlist()

    classifier = PWNClassifier(loader=loader)

    return {pwn:classifier.automatic_clasification(pwn) for pwn in pwnlist}


class PWNClassifier(object):
    """ Classify a PWN """

    def __init__(self, loader, pwn_classifications=None):
        self.loader = loader
        self.pwn_classifications = pwn_classifications

    def get_fit_parameters(pwn):

        classifier = self.manual_clasification(pwn)
        spatial_model=classifier['spatial_model']
        spectral_model=classifier['spectral_model']
        results = dict()
        return results


    def manual_clasification(self, pwn):
        assert self.pwn_classifications is not None
        return yaml.load(self.pwn_classifications)[pwn]

    def automatic_clasification(self, pwn):
        results = self.loader.get_results(pwn)

        point_gtlike = results['point']['gtlike']
        extended_gtlike = results['extended']['gtlike']
        cutoff=point_gtlike['test_cutoff']

        ts_point = max(point_gtlike['TS']['reoptimize'],0)
        ts_ext = max(extended_gtlike['TS']['reoptimize']-point_gtlike['TS']['reoptimize'],0)

        assert point_gtlike['model']['name'] == extended_gtlike['model']['name']
        spectral_name = point_gtlike['model']['name']

        try:
            ts_cutoff = max(cutoff['TS_cutoff'],0)
        except:
            ts_cutoff = None

        if ts_point > 25:

            if ts_cutoff > 16:
                source_class = 'pulsar'

                spectral_model = 'PLSuperExpCutoff'
                spatial_model = 'point'
            else:
                if ts_ext > 16:
                    source_class = 'confused'
                    spatial_model = 'extended'
                    spectral_model = spectral_name
                else:
                    source_class = 'confused'
                    spatial_model = 'point'
                    spectral_model = spectral_name

        else:
            source_class = 'upper_limit'
            spatial_model = 'at_pulsar'
            spectral_model = spectral_name

        return dict(source_class=source_class, 
                    spatial_model=spatial_model,
                    spectral_model=spectral_model)





class BestHypothesis(object):
    """ For legacy """
    def __init__(self, results):
        self.results = results

        at_pulsar_gtlike = results['at_pulsar']['gtlike']
        at_pulsar_pointlike = results['at_pulsar']['pointlike']
        

        point_gtlike = results['point']['gtlike']
        point_pointlike = results['point']['pointlike']
        
        extended_gtlike = results['extended']['gtlike']
        extended_pointlike = results['extended']['pointlike']

        self.cutoff=point_gtlike['test_cutoff']

        self.ts_point = max(point_gtlike['TS'],0)
        #self.ts_ext = max(extended_gtlike['ts_ext'],0)
        self.ts_ext = max(extended_gtlike['TS']-point_gtlike['TS'],0)

        try:
            self.ts_cutoff = max(self.cutoff['TS_cutoff'],0)
        except:
            self.ts_cutoff = None

        if self.ts_point > 25:
            if self.ts_ext > 16 and self.ts_cutoff > 16:
                self.type = 'confused'
            if self.ts_ext > 16:
                self.gtlike = extended_gtlike
                self.pointlike = extended_pointlike
                self.type = 'extended'
            elif self.ts_cutoff > 16:
                self.gtlike = point_gtlike
                self.pointlike = point_pointlike
                self.type = 'cutoff'
            else:
                self.gtlike = point_gtlike
                self.pointlike = point_pointlike
                self.type = 'point'
        else:
            self.type = 'upper_limit'
            self.gtlike = at_pulsar_gtlike
            self.pointlike = at_pulsar_pointlike

