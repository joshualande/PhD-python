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
    assert table_type in ['latex', 'confluence']

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
    if table_type == 'confluence':
        ts_point_name='TS_point'
        ts_ext_name='TS_ext'
        ts_cutoff_name = 'TS_cutoff'
        flux_name = 'Flux'
        index_name = 'Gamma'
        cutoff_name = 'E_cutoff'
    elif table_type == 'latex':
        ts_point_name=r'$\tspoint$'
        ts_ext_name=r'$\tsext$'
        ts_cutoff_name = r'$\tscutoff$'
        flux_name = r'Flux'
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

            table[ts_point_name].append(format.value(r['ts_point'],precision=1))
            table[ts_ext_name].append(format.value(r['ts_ext'],precision=1))
            table[ts_cutoff_name].append(format.value(r['ts_cutoff'],precision=1))

            table[flux_name].append(format.error(r['flux']/1e-9,r['flux_err']/1e-9))
            if r['spectral_model'] in ['PowerLaw','PLSuperExpCutoff']:
                table[index_name].append(format.error(r['index'],r['index_err']))
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
            table[flux_name].append('None')
            table[index_name].append('None')
            table[cutoff_name].append('None')

    table[psr_name][0] = '\multicolumn{8}{c}{Young Pulsars} \\\\\n\hline\n' + table[psr_name][0]
    table[psr_name][first_msp_index] = '\cutinhead{Millisecond Pulsars}\n' + table[psr_name][first_msp_index]
               
    writer=TableWriter(table, savedir, filebase)
    if table_type == 'confluence':
        writer.write_confluence(
                         units={
                             flux_name:r'(10^-9 ph cm^-2 s^-1)',
                             cutoff_name:r'(GeV)',
                         })
    elif table_type == 'latex':
        writer.write_latex(
                    preamble=r'\tabletypesize{\tiny}',
                    units={
                        flux_name:r'($10^{-9}\,\fluxunits$)',
                        cutoff_name:r'(GeV)',
                    },
                   )
