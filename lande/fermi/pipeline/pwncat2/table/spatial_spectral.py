from os.path import expandvars

import yaml
import numpy as np

from skymaps import SkyDir

from lande.utilities.tools import OrderedDefaultDict

from . writer import TableWriter
from . format import PWNFormatter
from lande.fermi.pipeline.pwncat2.interp.classify import PWNManualClassifier, PWNClassifierException
from lande.fermi.pipeline.pwncat2.interp.loader import PWNResultsLoader

from lande.fermi.pipeline.pwncat2.interp.bigfile import PulsarCatalogLoader

def spatial_spectral_table(pwndata, 
                           phase_shift, 
                           fitdir, savedir, pwn_classification, filebase, table_type,
                           bigfile_filename):
    assert table_type == 'latex'

    format=PWNFormatter(table_type=table_type, precision=2)

    loader = PWNResultsLoader(
        pwndata=pwndata,
        fitdir=fitdir,
        phase_shift=phase_shift
        )

    classifier = PWNManualClassifier(loader=loader, pwn_classification=pwn_classification)

    table = OrderedDefaultDict(list)

    psr_name='PSR'
    classification_name = 'Type'
    ts_point_name=r'$\tspoint$'
    ts_ext_name=r'$\tsext$'
    ts_cutoff_name = r'$\tscutoff$'
    ts_altdiff_name = r'$\tsaltdiff$'
    eflux_name = r'Energy Flux'
    index_name = r'$\Gamma$'
    cutoff_name = r'$\Ecutoff$'

    pwnlist = loader.get_pwnlist()
    #pwnlist = pwnlist[10:20]

    pcl = PulsarCatalogLoader(bigfile_filename=bigfile_filename)

    young = [ psr for psr in pwnlist if 'm' not in pcl.get_pulsar_classification(psr.replace('PSRJ','J'))]
    msps = [ psr for psr in pwnlist if 'm' in pcl.get_pulsar_classification(psr.replace('PSRJ','J'))]
    print 'young',young
    print 'msps',msps

    sorted_pwnlist = young + msps

    first_msp_index = None

    for pwn in sorted_pwnlist:
        print pwn

        try:
            r = classifier.get_results(pwn)

            if r['source_class'] == 'Upper_Limit': 
                continue

            if first_msp_index is None and pwn in msps:
                first_msp_index = len(table[psr_name])

            table[psr_name].append(format.pwn(pwn))
            table[classification_name].append(r['abbreviated_source_class'])

            def david_format_ts(x):
                if x >= 100:
                    return format.value(x,precision=0) + '.'
                else:
                    return format.value(x,precision=1)

            table[ts_point_name].append(david_format_ts(r['ts_point']))
            table[ts_ext_name].append(david_format_ts(r['ts_ext']))
            table[ts_cutoff_name].append(david_format_ts(r['ts_cutoff']))
            table[ts_altdiff_name].append(david_format_ts(r['ts_altdiff']) if r['ts_altdiff'] is not None else format.nodata)

            def david_format_flux(x, y):
                if x >= 10:
                    return format.error(x,y, precision=1)
                else:
                    return  format.error(x,y, precision=2)

            table[eflux_name].append(david_format_flux(r['energy_flux']/1e-11,r['energy_flux_err']/1e-11))
            if r['spectral_model'] in ['PowerLaw','PLSuperExpCutoff']:
                table[index_name].append(format.error(r['index'],r['index_err']))
            elif pwn == 'PSRJ0534+2200':
                table[index_name].append(r'\tablenotemark{a}')
            elif pwn == 'PSRJ0835-4510':
                table[index_name].append(r'\tablenotemark{b}')
            else:
                table[index_name].append(format.nodata)

            if r['spectral_model'] == 'PLSuperExpCutoff':
                table[cutoff_name].append(format.error(r['cutoff']/1e3,r['cutoff_err']/1e3, precision=2))
            else:
                table[cutoff_name].append(format.nodata)

        except PWNClassifierException, ex:
            print 'Skipping %s: %s' % (pwn,ex)
            table[psr_name].append(format.pwn(pwn))
            table[classification_name].append('None')
            table[ts_point_name].append('None')
            table[ts_ext_name].append('None')
            table[ts_cutoff_name].append('None')
            table[eflux_name].append('None')
            table[index_name].append('None')
            table[cutoff_name].append('None')

    table[psr_name][0] = '\multicolumn{8}{c}{Young Pulsars} \\\\\n\hline\n' + table[psr_name][0]
    table[psr_name][first_msp_index] = '\cutinhead{Millisecond Pulsars}\n' + table[psr_name][first_msp_index]
               
    writer=TableWriter(table, savedir, filebase)
    writer.write_latex(
                preamble=r'\tabletypesize{\tiny}',
                units={
                    eflux_name:r'($10^{-11}\,\efluxunits$)',
                    cutoff_name:r'(GeV)',
                },
               )
