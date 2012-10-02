from os.path import expandvars

import yaml
import numpy as np
import atpy

from skymaps import SkyDir

from lande.utilities.tools import OrderedDefaultDict

from . writer import TableWriter
from lande.fermi.pipeline.pwncat2.interp.classify import PWNClassifier,PWNManualClassifier
from lande.fermi.pipeline.pwncat2.interp.loader import PWNResultsLoader


def auxiliary_table(pwndata, phase_shift, fitdir, filename, pwn_classification):

    loader = PWNResultsLoader(
        pwndata=pwndata,
        fitdir=fitdir,
        phase_shift=phase_shift)

    classifier = PWNManualClassifier(loader=loader, pwn_classification=pwn_classification)

    pwnlist = loader.get_pwnlist()
    #pwnlist = pwnlist[10:20]

    npwn = len(pwnlist)

    table=atpy.Table(name='Off_Peak')


    def add_float(name, **kwargs):
        table.add_empty_column(name, np.float, shape=npwn, **kwargs)
        table[name][:]=np.nan
        return name

    def add_vector_float(name, size, *args, **kwargs):
        table.add_empty_column(name, np.float, shape=(npwn, size), **kwargs)
        table[name][:]=np.nan
        return name

    def add_string(name, width, *args, **kwargs):
        table.add_empty_column(name, np.dtype((str, width)), shape=npwn, *args, **kwargs)
        return name

    maxwidth=max([len(i) for i in pwnlist])
    psr_name=add_string('PSR', maxwidth)

    len_class = max(map(len,PWNClassifier.allowed_source_class))
    classification_name=add_string('Classification', len_class)


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

    energy_units = 'MeV'
    energy_flux_units = 'erg/cm^2/s'
    flux_units = 'ph/cm^2/s'
    prefactor_units = 'ph/cm^2/s/erg'

    energy_flux_name = add_float('Energy_Flux', unit=energy_flux_units)
    energy_flux_err_name = add_float('Energy_Flux_Error', unit=energy_flux_units)

    flux_name = add_float('Flux', unit=flux_units)
    flux_err_name = add_float('Flux_Error', unit=flux_units)

    prefactor_name = add_float('Prefactor', unit=prefactor_units)
    prefactor_err_name = add_float('prefactor_Error', unit=prefactor_units)

    normalization_name = add_float('Normalization')
    normalization_err_name = add_float('Normalization_Error')


    scale_name = add_float('Scale', unit=energy_units)

    index_name = add_float('Gamma')
    index_err_name = add_float('Gamma_Error')

    cutoff_name = add_float('Energy_Cutoff', unit=energy_units)
    cutoff_err_name = add_float('Energy_Cutoff_Error', unit=energy_units)

    # Spatial Stuff

    len_spatial = max(map(len,PWNClassifier.allowed_spatial_models))
    spatial_model_name=add_string('Spatial_Model',len_spatial)

    ra_name = add_float('RA_J2000', unit='deg')
    dec_name = add_float('DEC_J2000', unit='deg')

    glon_name = add_float('GLON', unit='deg')
    glat_name = add_float('GLAT', unit='deg')

    poserr_name = add_float('Position_Error', unit='deg')

    extension_name = add_float('Extension', unit='deg')
    extension_err_name = add_float('Extension_Error', unit='deg')

    powerlaw_flux_upper_limit_name = add_float('PowerLaw_Flux_Upper_Limit')
    powerlaw_energy_flux_upper_limit_name =add_float('PowerLaw_Energy_Flux_Upper_Limit')

    cutoff_flux_upper_limit_name = add_float('Cutoff_Flux_Upper_Limit')
    cutoff_energy_flux_upper_limit_name =add_float('Cutoff_Energy_Flux_Upper_Limit')

    sed_size=14
    sed_ts_name = add_vector_float('SED_TS', size=sed_size)
    sed_lower_energy_name = add_vector_float('SED_Lower_Energy', size=sed_size, unit=energy_units)
    sed_upper_energy_name = add_vector_float('SED_Upper_Energy', size=sed_size, unit=energy_units)
    sed_middle_energy_name = add_vector_float('SED_Middle_Energy', size=sed_size, unit=energy_units)
    sed_prefactor_name = add_vector_float('SED_Prefactor', size=sed_size, unit='ph/cm^2/s/MeV')
    sed_prefactor_lower_err_name = add_vector_float('SED_Prefactor_Lower_Error', size=sed_size, unit='ph/cm^2/s/erg')
    sed_prefactor_upper_err_name = add_vector_float('SED_Prefactor_Upper_Error', size=sed_size, unit='ph/cm^2/s/erg')
    sed_prefactor_upper_limit_name = add_vector_float('SED_Prefactor_Upper_Limit', size=sed_size, unit='ph/cm^2/s/erg')

    bandfits_size = 3
    bandfits_ts_name = add_vector_float('Bandfits_TS', size=bandfits_size)
    bandfits_lower_energy_name = add_vector_float('Bandfits_Lower_Energy', size=bandfits_size, unit=energy_units)
    bandfits_upper_energy_name = add_vector_float('Bandfits_Upper_Energy', size=bandfits_size, unit=energy_units)
    bandfits_middle_energy_name = add_vector_float('Bandfits_Middle_Energy', size=bandfits_size, unit=energy_units)

    bandfits_flux_name = add_vector_float('Bandfits_Flux', size=bandfits_size, unit=flux_units)
    bandfits_flux_err_name = add_vector_float('Bandfits_Flux_Error', size=bandfits_size, unit=flux_units)
    bandfits_flux_upper_limit_name = add_vector_float('Bandfits_Flux_Upper_Limit', size=bandfits_size, unit=flux_units)

    bandfits_energy_flux_name = add_vector_float('Bandfits_energy_flux', size=bandfits_size, unit=energy_flux_units)
    bandfits_energy_flux_err_name = add_vector_float('Bandfits_energy_flux_Error', size=bandfits_size, unit=energy_flux_units)
    bandfits_energy_flux_upper_limit_name = add_vector_float('Bandfits_energy_flux_upper_limit', size=bandfits_size, unit=energy_flux_units)

    bandfits_prefactor_name = add_vector_float('Bandfits_Prefactor', size=bandfits_size, unit=prefactor_units)
    bandfits_prefactor_err_name = add_vector_float('Bandfits_Prefactor_Error', size=bandfits_size, unit=prefactor_units)
    
    bandfits_index_name = add_vector_float('Bandfits_Index', size=bandfits_size)
    bandfits_index_err_name = add_vector_float('Bandfits_Index_Error', size=bandfits_size)

    for i,pwn in enumerate(pwnlist):
        print pwn

        r = classifier.get_results(pwn)

        if r is None: continue

        phase=r['shifted_phase']

        table[psr_name][i]=pwn

        source_class = r['source_class']
        table[classification_name][i] = source_class

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

        if source_class in ['Confused', 'Pulsar', 'PWN']:
            table[ts_ext_name][i]=r['ts_ext']
            table[ts_cutoff_name][i]=r['ts_cutoff']

        # spectral stuff

        table[spectral_model_name][i]=r['spectral_model']

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
        table[spatial_model_name][i] = spatial_model

        table[ra_name][i] = r['ra']
        table[dec_name][i] = r['dec']

        table[glon_name][i] = r['glon']
        table[glat_name][i] = r['glat']

        if spatial_model in [ 'Point', 'Extended' ]:
            table[poserr_name][i] = r['poserr']

        if spatial_model == 'Extended':
            table[extension_name][i] = r['extension']
            table[extension_err_name][i] = r['extension_err']

        if source_class == 'Upper_Limit':
            # Add powerlaw upper limit
            table[powerlaw_flux_upper_limit_name][i] = r['powerlaw_flux_upper_limit']
            table[powerlaw_energy_flux_upper_limit_name][i] = r['powerlaw_energy_flux_upper_limit']

            # Add cutoff upper limit
            table[cutoff_flux_upper_limit_name][i] = r['cutoff_flux_upper_limit']
            table[cutoff_energy_flux_upper_limit_name][i] = r['cutoff_energy_flux_upper_limit']

        # Add SED results

        table[sed_ts_name][i] = r['sed_ts'] 

        table[sed_lower_energy_name] = r['sed_lower_energy'] 
        table[sed_upper_energy_name] = r['sed_upper_energy'] 
        table[sed_middle_energy_name] = r['sed_middle_energy'] 

        table[sed_prefactor_name][i] = r['sed_prefactor'] 
        table[sed_prefactor_lower_err_name][i] = r['sed_prefactor_lower_err'] 
        table[sed_prefactor_upper_err_name][i] = r['sed_prefactor_upper_err'] 
        table[sed_prefactor_upper_limit_name][i] = r['sed_prefactor_upper_limit'] 

        # Add bandfit results
        # ...
        table[bandfits_ts_name] = r['bandfits_ts']

        table[bandfits_lower_energy_name] = r['bandfits_lower_energy']
        table[bandfits_upper_energy_name] = r['bandfits_upper_energy']
        table[bandfits_middle_energy_name] = r['bandfits_middle_energy']

        table[bandfits_flux_name] = r['bandfits_flux']
        table[bandfits_flux_err_name] = r['bandfits_flux_err']
        table[bandfits_flux_upper_limit_name] = r['bandfits_flux_upper_limit']

        table[bandfits_energy_flux_name] = r['bandfits_energy_flux']
        table[bandfits_energy_flux_err_name] = r['bandfits_energy_flux_err']
        table[bandfits_energy_flux_upper_limit_name] = r['bandfits_energy_flux_upper_limit']

        table[bandfits_prefactor_name] = r['bandfits_prefactor']
        table[bandfits_prefactor_err_name] = r['bandfits_prefactor_err']
        
        table[bandfits_index_name] = r['bandfits_index']
        table[bandfits_index_err_name] = r['bandfits_index_err']

    table.write(expandvars(filename), overwrite=True)
        
