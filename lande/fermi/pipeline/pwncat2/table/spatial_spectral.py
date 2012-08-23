from os.path import expandvars

import yaml
import numpy as np

from skymaps import SkyDir

from uw.pulsar.phase_range import PhaseRange

from lande.utilities.tools import OrderedDefaultDict

from . writer import TableWriter
from . format import PWNFormatter
from lande.fermi.pipeline.pwncat2.interp.classify import BestHypothesis
from lande.fermi.pipeline.pwncat2.interp.loader import ResultsLoader


def spatial_spectral_table(pwndata, fitdir, savedir, filebase, table_type):
    assert table_type in ['latex', 'confluence']

    format=PWNFormatter(table_type=table_type, precision=2)

    resloader = ResultsLoader(
        pwndata=pwndata,
        fitdir=fitdir)

    table = OrderedDefaultDict(list)

    psr_name='PSR'
    phase_name='Off Peak'
    if table_type == 'confluence':
        ts_point_name='TS_point'
        ts_ext_name='TS_ext'
        ts_cutoff_name = 'TS_cutoff'
        flux_name = 'F_(0.1-316)'
        index_name = 'Gamma'
        cutoff_name = 'E_cutoff'
    elif table_type == 'latex':
        ts_point_name=r'$\ts_\text{point}$'
        ts_ext_name=r'\tsext'
        ts_cutoff_name = r'$\ts_\text{cutoff}$'
        flux_name = r'$F_{0.1-316}$'
        index_name = r'$\Gamma$'
        cutoff_name = r'$E_\text{cutoff}$'

    data = yaml.load(open(expandvars('$pwncode/data/pwncat2_phase_lande.yaml')))

    pwnlist = resloader.get_pwnlist()

    for pwn in pwnlist:
        print pwn

        phase=PhaseRange(data[pwn]['phase'])


        results = resloader.get_results(pwn)

        if results is None:
            print 'Skipping %s' % pwn
            # job crashed/not finished
            table[psr_name].append(format.pwn(pwn))
            table[phase_name].append(phase.pretty_format())
            table[ts_point_name].append('None')
            table[ts_ext_name].append('None')
            table[ts_cutoff_name].append('None')
            table[flux_name].append('None')
            table[index_name].append('None')
            table[cutoff_name].append('None')
        else:

            point_gtlike = results['point']['gtlike']
            point_pointlike = results['point']['pointlike']

            extended_pointlike = results['extended']['pointlike']
            extended_gtlike = results['extended']['gtlike']
            

            b = BestHypothesis(results)
            gtlike = b.gtlike
            pointlike = b.pointlike
            type = b.type
            cutoff = b.cutoff

            ts_point = b.ts_point
            ts_ext = b.ts_ext
            ts_cutoff = b.ts_cutoff

            if type == 'upper_limit': 
                continue

            phase=PhaseRange(data[pwn]['phase'])

            table[psr_name].append(format.pwn(pwn))
            table[phase_name].append(phase.pretty_format())

            table[ts_point_name].append(format.value(ts_point,precision=1))
            table[ts_ext_name].append(format.value(ts_ext,precision=1))
            table[ts_cutoff_name].append(format.value(ts_cutoff,precision=1))

            if pwn == 'PSRJ0534+2200':
                table[flux_name].append(format.error(gtlike['flux']['flux']/1e-9,gtlike['flux']['flux_err']/1e-9))
                table[index_name].append(format.nodata)
                table[cutoff_name].append(format.nodata)
            elif type == 'cutoff':
                table[flux_name].append(format.error(cutoff['flux_1']['flux']/1e-9,cutoff['flux_1']['flux_err']/1e-9))
                table[index_name].append(format.error(cutoff['model_1']['Index1'],cutoff['model_1']['Index1_err']))
                table[cutoff_name].append(format.error(cutoff['model_1']['Cutoff'],cutoff['model_1']['Cutoff_err']))
            elif type in ['extended','point']:
                table[flux_name].append(format.error(gtlike['flux']['flux']/1e-9,gtlike['flux']['flux_err']/1e-9))
                table[index_name].append(format.error(-1*gtlike['model']['Index'],-1*gtlike['model']['Index_err']))
                table[cutoff_name].append(format.nodata)
            elif type == 'confused':
                table[flux_name].append('???')
                table[index_name].append('???')
                table[cutoff_name].append('???')


                
    writer=TableWriter(table, savedir, filebase)
    if table_type == 'confluence':
        writer.write_confluence(
                         units={
                             flux_name:r'(10^-9 erg cm^-2 s^-1)',
                             cutoff_name:r'(MeV)',
                         })
    elif table_type == 'latex':
        writer.write_latex(
                    preamble=r'\tabletypesize{\tiny}',
                    units={
                        flux_name:r'($10^{-9}$\ erg\,cm$^{-2}$\,s$^{-1}$)',
                        cutoff_name:r'(MeV)',
                    },
                   )
