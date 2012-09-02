from os.path import expandvars

import yaml
import numpy as np
import atpy

from skymaps import SkyDir

from uw.pulsar.phase_range import PhaseRange

from lande.utilities.tools import OrderedDefaultDict

from . writer import TableWriter
from lande.pipeline.pwncat2.interp.classify import PWNClassifier
from lande.pipeline.pwncat2.interp.loader import PWNResultsLoader


def auxiliary_table(pwndata, fitdir, filename, pwn_classification):

    results_loader = PWNResultsLoader(
        pwndata=pwndata,
        fitdir=fitdir)

    classifier = PWNClassifier(results_loader, pwn_classification)

    pwnlist = results_loader.get_pwnlist()
    pwnlist = pwnlist[0:15]

    table=atpy.Table(name='Off_Peak')


    def add_float(name, *args, **kwargs):
        table.add_empty_column(name, np.float, *args, **kwargs); return name
        table[name][:]=np.nan
        return name

    def add_string(name, width, *args, **kwargs):
        table.add_empty_column(name, np.dtype((str, width)), *args, **kwargs)
        return name

    maxwidth=max([len(i) for i in pwnlist])
    psr_name=add_string('PSR', maxwidth, shape=len(pwnlist))

    len_class = max(map(len,PWNClassifier.allowed_source_class))
    source_class_name=add_string('Source_Class', len_class)


    # Phase Stuff

    off_peak_min_name=add_float('Off_Peak_Min')
    off_peak_max_name=add_float('Off_Peak_Max')
    second_off_peak_min_name=add_float('Second_Off_Peak_Min')
    second_off_peak_max_name=add_float('Second_Off_Peak_Max')

    # Significance stuff
    ts_point_name=add_float('TS_point')
    ts_ext_name=add_float('TS_ext')
    ts_cutoff_name=add_float('TS_cutoff')

    # Spectral Stuff

    len_spectral = max(map(len,PWNClassifier.allowed_spectral_models))
    spectral_model_name=add_string('Spectral_Model', len_spectral)

    flux_name = add_float('Flux_(0.1-316)', unit='ph/cm^2/s')
    flux_err_name = add_float('Flux_err_(0.1-316)', unit='ph/cm^2/s')

    prefactor_name = add_float('Prefactor', unit='ph/cm^2/s/MeV')
    prefactor_err_name = add_float('Prefactor_err', unit='ph/cm^2/s/MeV')

    scale_name = add_float('Scale', unit='MeV')

    index_name = add_float('Gamma')
    index_err_name = add_float('Gamma_err')

    cutoff_name = add_float('Energy_cutoff', unit='MeV')
    cutoff_err_name = add_float('Energy_err_cutoff', unit='MeV')

    data = yaml.load(open(expandvars('$pwncode/data/pwncat2_phase_lande.yaml')))

    # Spatial Stuff

    len_spatial = max(map(len,PWNClassifier.allowed_spatial_models))
    spatial_model_name=add_string('Spatial_Model',len_spatial)

    ra_name = add_float('RAJ2000', unit='deg')
    dec_name = add_float('DECJ2000', unit='deg')

    glon_name = add_float('GLON', unit='deg')
    glat_name = add_float('GLAT', unit='deg')

    poserr_name = add_float('Position_Error', unit='deg')

    extension_name = add_float('Extension', unit='deg')
    extension_err_name = add_float('Extension_Err', unit='deg')

    for i,pwn in enumerate(pwnlist):
        print pwn

        phase=PhaseRange(data[pwn]['phase'])

        r = classifier.get_results(pwn)

        if r is None:
            continue

        phase=PhaseRange(data[pwn]['phase'])

        table[psr_name][i]=pwn

        source_class = r['source_class']
        table[source_class_name][i] = source_class

        if phase.is_continuous():
            a,b = phase.tolist(dense=True)
            table[off_peak_min_name][i]=a
            table[off_peak_max_name][i]=b
        else:
            ranges = phase.split_ranges()
            assert len(ranges) == 2
            a,b = ranges[0].tolist(dense=True)
            table[off_peak_min_name][i]=a
            table[off_peak_max_name][i]=b
            a,b = ranges[1].tolist(dense=True)
            table[second_off_peak_min_name][i]=a
            table[second_off_peak_max_name][i]=b

        # likelihood stuff

        table[ts_point_name][i]=r['ts_point']

        if source_class in ['Confused', 'Pulsar', 'PWN']:
            table[ts_ext_name][i]=r['ts_ext']
            table[ts_cutoff_name][i]=r['ts_cutoff']

        # spectral stuff

        if source_class in ['Confused', 'Pulsar', 'PWN']:

            table[flux_name][i]=r['flux']
            table[flux_err_name][i]=r['flux_err']

            table[spectral_model_name][i]=r['spectral_model']

            if r['spectral_model'] in ['PowerLaw','PLSuperExpCutoff']:
                # Don't include prefactor for FileFunction model, since the units are wrong

                table[prefactor_name][i]=r['prefactor']
                table[prefactor_err_name][i]=r['prefactor_err']

                table[index_name][i]=r['index']
                table[index_err_name][i]=r['index_err']

                table[scale_name][i]=r['model_scale']

            if r['spectral_model'] == 'PLSuperExpCutoff':

                table[cutoff_name][i]=r['cutoff']
                table[cutoff_err_name][i]=r['cutoff_err']

        # spatial stuff

        spatial_model = r['spatial_model']
        table[spatial_model_name][i] = spatial_model

        table[ra_name][i] = r['ra']
        table[dec_name][i] = r['dec']

        table[glon_name][i] = r['glon']
        table[glat_name][i] = r['glat']

        if spatial_model in [ 'Point', 'Extended' ]:
            table[poserr_name][i] = r['poserr']

        """
        if spatial_model == 'Extended':
            table[extension_name][i] = r['extension']
            table[extension_err_name][i] = r['extension_err']
        """

    table.write(expandvars(filename), overwrite=True)
        
