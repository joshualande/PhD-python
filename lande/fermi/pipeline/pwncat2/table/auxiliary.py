from os.path import expandvars

import yaml
import numpy as np
import atpy

from skymaps import SkyDir

from lande.utilities.tools import OrderedDefaultDict

from . writer import TableWriter
from lande.fermi.pipeline.pwncat2.interp.classify import PWNClassifier,PWNManualClassifier, PWNClassifierException
from lande.fermi.pipeline.pwncat2.interp.loader import PWNResultsLoader


def auxiliary_table(pwndata, 
                    phase_shift, 
                    fitdir, filename, pwn_classification):

    loader = PWNResultsLoader(
        pwndata=pwndata,
        fitdir=fitdir,
        phase_shift=phase_shift
        )

    classifier = PWNManualClassifier(loader=loader, pwn_classification=pwn_classification)

    pwnlist = loader.get_pwnlist()

    npwn = len(pwnlist)

    table=atpy.Table(name='Off_Peak')


    def add_float(name, **kwargs):
        table.add_empty_column(name, np.dtype('float32'), shape=npwn, **kwargs)
        table[name][:]=np.nan
        return name

    def add_int(name, **kwargs):
        table.add_empty_column(name, np.dtype('uint32'), shape=npwn, **kwargs)
        table[name][:]=np.nan
        return name

    def add_vector_float(name, size, *args, **kwargs):
        table.add_empty_column(name, np.dtype('float32'), shape=(npwn, size), **kwargs)
        table[name][:]=np.nan
        return name

    def add_vector_int(name, size, *args, **kwargs):
        table.add_empty_column(name, np.dtype('uint32'), shape=(npwn, size), **kwargs)
        table[name][:]=np.nan
        return name

    def add_string(name, width, *args, **kwargs):
        table.add_empty_column(name, np.dtype((str, width)), shape=npwn, *args, **kwargs)
        return name

    maxwidth=max([len(i) for i in pwnlist])
    psr_name=add_string('PSR', maxwidth)

    len_class = max(map(len,PWNClassifier.abbreviated_source_class_mapper.values()))
    classification_name=add_string('Classification_OP', len_class)


    # Phase Stuff
    off_peak_min_name=add_float('Min_Phase_OP')
    off_peak_max_name=add_float('Max_Phase_OP')
    second_off_peak_min_name=add_float('Min_2_Phase_OP')
    second_off_peak_max_name=add_float('Max_2_Phase_OP')

    # Significance stuff
    ts_point_name=add_float('TS_point_OP')
    ts_ext_name=add_float('TS_ext_OP')
    ts_cutoff_name=add_float('TS_cutoff_OP')
    ts_altdiff_name = add_float(r'TS_altdiff_OP')
    ts_var_name=add_float('TS_var_OP')

    # Spectral Stuff

    len_spectral = max(map(len,PWNClassifier.allowed_spectral_models))
    spectral_model_name=add_string('Spectral_Model_OP', len_spectral)

    energy_units = 'MeV'
    energy_flux_units = 'erg/cm^2/s'
    flux_units = 'ph/cm^2/s'
    prefactor_units = 'ph/cm^2/s/erg'

    flux_name = add_float('Flux_OP', unit=flux_units)
    flux_err_name = add_float('Unc_Flux_OP', unit=flux_units)

    energy_flux_name = add_float('EFlux_OP', unit=energy_flux_units)
    energy_flux_err_name = add_float('Unc_EFlux_OP', unit=energy_flux_units)


    prefactor_name = add_float('Prefactor_OP', unit=prefactor_units)
    prefactor_err_name = add_float('Unc_Prefactor_OP', unit=prefactor_units)

    normalization_name = add_float('Normalization_OP')
    normalization_err_name = add_float('Unc_Normalization_OP')


    scale_name = add_float('Scale_OP', unit=energy_units)

    index_name = add_float('Index_OP')
    index_err_name = add_float('Unc_Index_OP')

    cutoff_name = add_float('Energy_Cutoff_OP', unit=energy_units)
    cutoff_err_name = add_float('Unc_Energy_Cutoff_OP', unit=energy_units)

    # Spatial Stuff

    len_spatial = max(map(len,PWNClassifier.allowed_spatial_models))
    spatial_model_name=add_string('Spatial_Model_OP',len_spatial)

    ra_name = add_float('RAJ2000_OP', unit='deg')
    dec_name = add_float('DECJ2000_OP', unit='deg')

    glon_name = add_float('GLON_OP', unit='deg')
    glat_name = add_float('GLAT_OP', unit='deg')

    poserr_name = add_float('Unc_Position_OP', unit='deg')

    extension_name = add_float('Extension_OP', unit='deg')
    extension_err_name = add_float('Unc_Extension_OP', unit='deg')

    powerlaw_flux_upper_limit_name = add_float('PowerLaw_Flux_UL_OP')
    powerlaw_energy_flux_upper_limit_name =add_float('PowerLaw_EFlux_UL_OP')

    cutoff_flux_upper_limit_name = add_float('Cutoff_Flux_UL_OP')
    cutoff_energy_flux_upper_limit_name =add_float('Cutoff_EFlux_UL_OP')

    sed_size=14
    sed_lower_energy_name = add_vector_float('SED_Lower_Energy_OP', size=sed_size, unit=energy_units)
    sed_upper_energy_name = add_vector_float('SED_Upper_Energy_OP', size=sed_size, unit=energy_units)
    sed_middle_energy_name = add_vector_float('SED_Center_Energy_OP', size=sed_size, unit=energy_units)

    sed_ts_name = add_vector_float('SED_TS_OP', size=sed_size)

    sed_prefactor_name = add_vector_float('SED_Prefactor_OP', size=sed_size, unit='ph*cm**-2*s**-1*erg**-1')
    sed_prefactor_lower_err_name = add_vector_float('SED_Neg_Unc_Prefactor_OP', size=sed_size, unit='ph*cm**-2*s**-1*erg**-1')
    sed_prefactor_upper_err_name = add_vector_float('SED_Pos_Unc_Prefactor_OP', size=sed_size, unit='ph*cm**-2*s**-1*erg**-1')
    sed_prefactor_upper_limit_name = add_vector_float('SED_Prefactor_UL_OP', size=sed_size, unit='ph*cm**-2*s**-1*erg**-1')

    for i,pwn in enumerate(pwnlist):
        print pwn

        try:
            r = classifier.get_results(pwn)
        except PWNClassifierException:
            print 'Skipping %s' % pwn
            continue

        phase=r['shifted_phase']

        table[psr_name][i]=pwn.replace('PSRJ','J')

        source_class = r['source_class']
        table[classification_name][i] = r['abbreviated_source_class']

        assert source_class in PWNClassifier.allowed_source_class

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
        ts_point = r['ts_point']
        table[ts_point_name][i]=ts_point

        if source_class in ['Confused', 'Pulsar', 'Pulsar_Confused', 'PWN']:
            table[ts_ext_name][i]=r['ts_ext']
            table[ts_cutoff_name][i]=r['ts_cutoff']
        elif source_class == 'Upper_Limit':
            pass
        else:
            raise Exception("...")

        if source_class in ['Pulsar','Pulsar_Confused']:
            table[ts_altdiff_name][i]=r['ts_altdiff']
        elif source_class in ['Confused', 'PWN', 'Upper_Limit']:
            pass
        else:
            raise Exception("...")

        table[ts_var_name][i]=r['ts_var']

        # spectral stuff

        spectral_model = r['spectral_model']
        if not type(spectral_model) == str and np.isnan(spectral_model):
            spectral_model = 'NULL'

        table[spectral_model_name][i]= spectral_model

        table[energy_flux_name][i] = r['energy_flux']
        table[energy_flux_err_name][i] = r['energy_flux_err']

        table[flux_name][i]=r['flux']
        table[flux_err_name][i]=r['flux_err']

        table[prefactor_name][i]=r['prefactor']
        table[prefactor_err_name][i]=r['prefactor_err']

        table[normalization_name][i]=r['normalization']
        table[normalization_err_name][i]=r['normalization_err']

        table[index_name][i]=r['index']
        table[index_err_name][i]=r['index_err']

        table[scale_name][i]=r['model_scale']

        table[cutoff_name][i]=r['cutoff']
        table[cutoff_err_name][i]=r['cutoff_err']

        # spatial stuff

        spatial_model = r['spatial_model']
        if not type(spatial_model) == str and np.isnan(spatial_model):
            spatial_model = 'NULL'
        table[spatial_model_name][i] = spatial_model

        table[ra_name][i] = r['ra']
        table[dec_name][i] = r['dec']

        table[glon_name][i] = r['glon']
        table[glat_name][i] = r['glat']

        table[poserr_name][i] = r['poserr']

        table[extension_name][i] = r['extension']
        table[extension_err_name][i] = r['extension_err']

        # Add powerlaw upper limit
        table[powerlaw_flux_upper_limit_name][i] = r['powerlaw_flux_upper_limit']
        table[powerlaw_energy_flux_upper_limit_name][i] = r['powerlaw_energy_flux_upper_limit']

        # Add cutoff upper limit
        table[cutoff_flux_upper_limit_name][i] = r['cutoff_flux_upper_limit']
        table[cutoff_energy_flux_upper_limit_name][i] = r['cutoff_energy_flux_upper_limit']

        # Add SED results


        table[sed_lower_energy_name] = r['sed_lower_energy'] 
        table[sed_upper_energy_name] = r['sed_upper_energy'] 
        table[sed_middle_energy_name] = r['sed_middle_energy'] 

        table[sed_ts_name][i] = r['sed_ts'] 

        table[sed_prefactor_name][i] = r['sed_prefactor'] 
        table[sed_prefactor_lower_err_name][i] = r['sed_prefactor_lower_err'] 
        table[sed_prefactor_upper_err_name][i] = r['sed_prefactor_upper_err'] 
        table[sed_prefactor_upper_limit_name][i] = r['sed_prefactor_upper_limit'] 

    table.write(expandvars(filename), overwrite=True)
        
