""" Code to clasify regions. """
import copy
from os.path import expandvars
import math
import itertools

import numpy as np
import yaml

from lande.pysed import units

from . loader import PWNResultsLoader, PWNResultsException


isnan = lambda x: isinstance(x,float) and math.isnan(x)
equal = lambda x,y: (x==y) or (isnan(x) and isnan(y))

class PWNClassifierException(Exception):
    pass

def compare_classifications(pwndata,fitdir,pwn_classification):

    loader = PWNResultsLoader(pwndata, fitdir, verbosity=False)

    auto=PWNAutomaticClassifier(loader)
    manual=PWNManualClassifier(loader,pwn_classification)

    pwnlist = loader.get_pwnlist()

    for pwn in pwnlist:
        try:
            manual_results = manual.get_classification(pwn)

            try:
                auto_results = auto.get_classification(pwn)


                if equal(auto_results['spatial_model'],manual_results['spatial_model']) and \
                   equal(auto_results['spectral_model'],manual_results['spectral_model']) and \
                   auto_results['source_class'] == manual_results['source_class']:
                    print '%s - classifications agree' % pwn
                else:
                    print '%s - classifications do not agree' % pwn
                    print '                    %20s %20s' % ('automatic','manual')
                    print '  - spectral model: %20s %20s' % (auto_results['spectral_model'],manual_results['spectral_model'])
                    print '  - spatial model:  %20s %20s' % (auto_results['spatial_model'],manual_results['spatial_model'])
                    print '  - source class:   %20s %20s' % (auto_results['source_class'],manual_results['source_class'])
            except PWNClassifierException:
                print "%s - no automatic classification" % pwn

        except PWNClassifierException:
            print "%s - no manual classification" % pwn


class PWNClassifier(object):

    abbreviated_source_class_mapper = dict(
        Pulsar='M', # Emission "Magnetospheric"
        Pulsar_Confused='M*',
        PWN='W', # From Pulsar "Wind"
        Confused='U',
        Upper_Limit='L')

    """
    expanded_source_class_mapper = dict(
        Pulsar='Magnetospheric', # Emission "Magnetospheric"
        Pulsar_Confused='Magnetospheric*',
        PWN='PWN', # From Pulsar "Wind"
        Confused='Unidentified',
        Upper_Limit='Upper_Limit')
    """

    allowed_source_class = ['Pulsar', 'Pulsar_Confused', 'PWN', 'Confused', 'Upper_Limit']
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

        assert source_class in PWNClassifier.allowed_source_class
        assert source_class == 'Upper_Limit' or spatial_model in PWNClassifier.allowed_spatial_models
        assert source_class == 'Upper_Limit' or spectral_model in PWNClassifier.allowed_spectral_models 

        results = self.loader.get_results(pwn, require_all_exists=True, get_variability=True)

        if results is None:
            print 'Results for %s is not done yet, skipping' % pwn
            return None

        point_gtlike = results['point']['gtlike']
        extended_gtlike = results['extended']['gtlike']

        if isnan(spatial_model): # upper limits
            gtlike = results['at_pulsar']['gtlike']
            pointlike = results['at_pulsar']['pointlike']
        else:
            gtlike = results[spatial_model.lower()]['gtlike']
            pointlike = results[spatial_model.lower()]['pointlike']

        at_pulsar_cutoff=results['at_pulsar']['gtlike']['test_cutoff']

        d = copy.copy(classifier)

        d['raw_phase'] = results['raw_phase']

        if results.has_key('shifted_phase'):
            d['shifted_phase'] = results['shifted_phase']

        # likelihood stuff

        d['ts_point'] = max(point_gtlike['TS']['reoptimize'],0)

        d['abbreviated_source_class'] = self.abbreviated_source_class_mapper[source_class]

        if source_class in ['Confused', 'Pulsar', 'Pulsar_Confused', 'PWN']:
            d['ts_ext'] = max(extended_gtlike['TS']['reoptimize']-point_gtlike['TS']['reoptimize'],0)
            d['ts_cutoff'] = max(at_pulsar_cutoff['hypothesis_1']['TS']['reoptimize']-at_pulsar_cutoff['hypothesis_0']['TS']['reoptimize'],0)

            alt_models = [point_gtlike['altdiff'][dist,halo,TS] for dist,halo,TS in itertools.product(['SNR','Lorimer'],[4,10],[150,100000])]
            if np.any([i is None for i in alt_models]):
                d['ts_altdiff'] = None
                print 'BAD = '
            else:
                all_TS = [ i['TS']['reoptimize'] if i is not None else None for i in alt_models]
                d['ts_altdiff'] = max(min(all_TS),0)
                print d['ts_point'], all_TS, d['ts_altdiff'] 

        elif source_class == 'Upper_Limit':
            pass
        else:
            raise Exception("...")

        if spatial_model == 'Extended':
            # note, does not make sense to use extended spatial model for variability
            d['ts_var'] = results['point']['variability']['TS_var']['gtlike']
        elif isinstance(spatial_model,float) and np.isnan(spatial_model):
            # upper limit
            d['ts_var'] = results['at_pulsar']['variability']['TS_var']['gtlike']
        else:
            d['ts_var'] = results[spatial_model.lower()]['variability']['TS_var']['gtlike']

        # spectral stuff

        d['flux'] = None
        d['flux_err'] = None

        d['energy_flux'] = None
        d['energy_flux_err'] = None

        d['prefactor'] = None
        d['prefactor_err'] = None

        d['normalization'] = None
        d['normalization_err'] = None

        d['index'] = None
        d['index_err'] = None

        d['model_scale'] = None

        d['cutoff'] = None
        d['cutoff_err'] = None

        convert_prefactor = lambda x: units.convert(x, 'ph/cm^2/s/MeV', 'ph/cm^2/s/erg')

        if source_class != 'Upper_Limit':

            if spectral_model in ['PowerLaw','FileFunction']:
                d['flux'] = gtlike['flux']['flux']
                d['flux_err'] = gtlike['flux']['flux_err']

                assert gtlike['flux']['flux_units'] == 'ph/cm^2/s'
                d['energy_flux'] = gtlike['flux']['eflux']
                d['energy_flux_err'] = gtlike['flux']['eflux_err']
                assert gtlike['flux']['eflux_units'] == 'erg/cm^2/s'

                d['spectrum'] = gtlike['spectrum']

                if spectral_model == 'PowerLaw':

                    # Note, prefactor is 
                    d['prefactor'] = convert_prefactor(gtlike['spectrum']['Prefactor'])
                    d['prefactor_err'] = convert_prefactor(gtlike['spectrum']['Prefactor_err'])

                    d['index'] = -1*gtlike['spectrum']['Index']
                    d['index_err'] = np.abs(gtlike['spectrum']['Index_err'])

                    d['model_scale'] = gtlike['spectrum']['Scale']

                elif spectral_model == 'FileFunction':

                    d['normalization'] = gtlike['spectrum']['Normalization']
                    d['normalization_err'] = gtlike['spectrum']['Normalization_err']

            elif spectral_model == 'PLSuperExpCutoff':
                h1 = gtlike['test_cutoff']['hypothesis_1']

                d['spectrum'] = h1['spectrum']

                d['flux'] = h1['flux']['flux']
                d['flux_err'] = h1['flux']['flux_err']

                d['energy_flux'] = h1['flux']['eflux']
                d['energy_flux_err'] = h1['flux']['eflux_err']

                assert h1['flux']['flux_units'] == 'ph/cm^2/s'
                assert h1['flux']['eflux_units'] == 'erg/cm^2/s'
            
                d['prefactor'] = convert_prefactor(h1['spectrum']['Prefactor'])
                d['prefactor_err'] = convert_prefactor(h1['spectrum']['Prefactor_err'])

                d['index'] = -1*h1['spectrum']['Index1']
                d['index_err'] = np.abs(h1['spectrum']['Index1_err'])

                d['model_scale'] = h1['spectrum']['Scale']

                d['cutoff'] = h1['spectrum']['Cutoff']
                d['cutoff_err'] = h1['spectrum']['Cutoff_err']
        else:
            # add in upper limits
            pass

        # spatial stuff

        d['ra'] = pointlike['position']['equ'][0]
        d['dec'] = pointlike['position']['equ'][1]

        d['glon'] = pointlike['position']['gal'][0]
        d['glat'] = pointlike['position']['gal'][1]

        d['poserr'] = None

        if spatial_model in [ 'Point', 'Extended' ]: 

            ellipse = pointlike['spatial_model']['ellipse']
            if ellipse.has_key('lsigma'):
                d['poserr'] = ellipse['lsigma']
            else:
                print 'WARNING: localization failed for %s' % pwn
                d['poserr'] = None

        d['extension'] = None
        d['extension_err'] = None
        if spatial_model == 'Extended':
            d['extension'] = pointlike['spatial_model']['Sigma']
            d['extension_err'] = pointlike['spatial_model']['Sigma_err']

        d['powerlaw_flux_upper_limit'] = None
        d['powerlaw_energy_flux_upper_limit'] = None

        d['cutoff_flux_upper_limit'] = None
        d['cutoff_energy_flux_upper_limit'] = None

        if source_class in ['Confused', 'Pulsar', 'Pulsar_Confused', 'PWN']:
            pass
        elif source_class == 'Upper_Limit':
            d['powerlaw_flux_upper_limit'] = gtlike['powerlaw_upper_limit']['flux']
            d['powerlaw_energy_flux_upper_limit'] = gtlike['powerlaw_upper_limit']['eflux']

            assert gtlike['powerlaw_upper_limit']['flux_units'] == 'ph/cm^2/s'
            assert gtlike['powerlaw_upper_limit']['eflux_units'] == 'erg/cm^2/s'

            d['cutoff_flux_upper_limit'] = gtlike['cutoff_upper_limit']['flux']
            d['cutoff_energy_flux_upper_limit'] = gtlike['cutoff_upper_limit']['eflux']
            
            assert gtlike['cutoff_upper_limit']['flux_units'] == 'ph/cm^2/s'
            assert gtlike['cutoff_upper_limit']['eflux_units'] == 'erg/cm^2/s'
        else:
            raise Exception("...")

        # add in 4BPD SED
        sed = gtlike['seds']['4bpd']
        sed_ts = np.asarray(sed['Test_Statistic'])
        sed_ts = np.where(sed_ts>0,sed_ts,0)
        sed_lower_energy = np.asarray(sed['Energy']['Lower'])
        sed_upper_energy = np.asarray(sed['Energy']['Upper'])
        sed_middle_energy = np.asarray(sed['Energy']['Value'])
        sed_prefactor = np.asarray(sed['dNdE']['Value'])
        sed_prefactor_lower_err = np.asarray(sed['dNdE']['Lower_Error'])
        sed_prefactor_upper_err = np.asarray(sed['dNdE']['Upper_Error'])
        sed_prefactor_upper_limit = np.asarray(sed['dNdE']['Upper_Limit'])

        assert sed['Energy']['Units'] == 'MeV'
        assert sed['dNdE']['Units'] == 'ph/cm^2/s/erg'

        d['sed_4bpd'] = sed

        # Note, when overall source is not significant, do not include upper limits
        d['sed_ts'] = sed_ts
        d['sed_lower_energy'] = sed_lower_energy
        d['sed_upper_energy'] = sed_upper_energy
        d['sed_middle_energy'] = sed_middle_energy
        if source_class in ['Confused', 'Pulsar', 'Pulsar_Confused', 'PWN']:
            significant = (sed_ts >= 4)
            d['sed_prefactor'] = np.where(significant, sed_prefactor, np.nan)
            d['sed_prefactor_lower_err'] = np.where(significant, sed_prefactor_lower_err, np.nan)
            d['sed_prefactor_upper_err'] = np.where(significant, sed_prefactor_upper_err, np.nan)
            d['sed_prefactor_upper_limit'] = np.where(~significant, sed_prefactor_upper_limit, np.nan)
        elif source_class == 'Upper_Limit':
            array_from = lambda x: np.asarray([x]*len(sed_ts))
            d['sed_prefactor'] = array_from(np.nan)
            d['sed_prefactor_lower_err'] = array_from(np.nan)
            d['sed_prefactor_upper_err'] = array_from(np.nan)
            d['sed_prefactor_upper_limit'] = sed_prefactor_upper_limit
        else:
            raise Exception("...")

        return d


class PWNAutomaticClassifier(PWNClassifier):

    def __init__(self, loader):
        self.loader = loader

    def get_classification(self, pwn):

        try:
            results = self.loader.get_results(pwn, require_all_exists=True, get_variability=False)
        except PWNResultsException:
            raise PWNClassifierException("No results for %s" % pwn)

        at_pulsar_gtlike = results['at_pulsar']['gtlike']
        point_gtlike = results['point']['gtlike']
        extended_gtlike = results['extended']['gtlike']
        at_pulsar_cutoff=at_pulsar_gtlike['test_cutoff']

        ts_point = max(point_gtlike['TS']['reoptimize'],0)
        ts_ext = max(extended_gtlike['TS']['reoptimize']-point_gtlike['TS']['reoptimize'],0)

        assert point_gtlike['spectrum']['name'] == extended_gtlike['spectrum']['name']
        spectral_name = point_gtlike['spectrum']['name']

        try:
            ts_cutoff = max(at_pulsar_cutoff['hypothesis_1']['TS']['reoptimize']-at_pulsar_cutoff['hypothesis_0']['TS']['reoptimize'],0)

        except:
            ts_cutoff = None

        if ts_point >= 25:

            if ts_cutoff >= 9:
                spatial_model = 'At_Pulsar'
                spectral_model = 'PLSuperExpCutoff'
                if ts_ext >= 16:
                    source_class = 'Pulsar_Confused'
                else:
                    source_class = 'Pulsar'
            else:
                if ts_ext >= 16:
                    spatial_model = 'Extended'
                    source_class = 'Confused'
                    spectral_model = spectral_name

                else:
                    spatial_model = 'Point'
                    source_class = 'Confused'
                    spectral_model = spectral_name

        else:
            source_class = 'Upper_Limit'
            spatial_model = np.nan
            spectral_model = np.nan

        return dict(source_class=source_class, 
                    spatial_model=spatial_model,
                    spectral_model=spectral_model)

    @staticmethod
    def get_automatic_classify(pwndata, fitdir):
        loader = PWNResultsLoader(pwndata, fitdir)
        pwnlist = loader.get_pwnlist()

        classifier = PWNAutomaticClassifier(loader=loader)

        return {pwn:classifier.get_classification(pwn) for pwn in pwnlist if loader.all_exists(pwn, get_variability=False)}


class PWNManualClassifier(PWNClassifier):
    """ Classify a PWN """


    def __init__(self, loader, pwn_classification):
        self.loader = loader
        self.pwn_classifications = yaml.load(open(expandvars(pwn_classification)))

    def get_classification(self, pwn):
        if pwn in self.pwn_classifications:
            return self.pwn_classifications[pwn]
        else:
            raise PWNClassifierException("No classification for %s" % pwn)

    @staticmethod
    def get_manual_classify(pwndata, fitdir):
        """ Create an empty dict which can be used for manual classification. 
            The acutal manual classification has to be created in by hand (duh). """
        loader = PWNResultsLoader(pwndata, fitdir)
        pwnlist = loader.get_pwnlist()
        return {pwn:dict(source_class=None, 
                        spatial_model=None,
                        spectral_model=None) for pwn in pwnlist}
