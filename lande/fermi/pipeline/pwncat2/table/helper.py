from os.path import join as j
from os.path import expandvars
import StringIO
from textwrap import dedent
import shutil
import os.path

import yaml

from lande.utilities.tools import merge_dict

from lande.utilities.table import latex_table,confluence_table,TableFormatter

base='$pwndata/spectral/v29/'

fitdir=expandvars(j(base,'analysis/'))
savedir=expandvars(j(base,'tables'))

if not os.path.exists(savedir): os.makedirs(savedir)


class PWNFormatter(TableFormatter):
    def pwn(self, pwn):
        pwn = pwn.replace('PSR','')
        if not self.confluence:
            pwn = pwn.replace('-','$-$')
        return pwn


def write_confluence(table, filebase, **kwargs):

    t = confluence_table(table, **kwargs)
    os.chdir(savedir)
    open('%s.confluence' % filebase,'w').write(t)


def write_latex(table, filebase, preamble='',**kwargs):

    t = latex_table(table, **kwargs)

    lines = t.split('\n')
    if lines[-1] == '': 
        lines=lines[0:-1]

    header = lines[0]
    footer= lines[-1]

    t = '\n'.join(lines[1:-1])

    os.chdir(savedir)

    open('%s.tex' % filebase,'w').write(t)

    open('temp.tex','w').write(dedent(r"""
        \documentclass{aastex}
        \usepackage{amsmath}

        \input{$pwnpaper/style/style.tex}

        \begin{document}
        %s
        %s
        \input{%s}
        %s
        \end{document}""" % (header,preamble,filebase,footer)))

    os.system('pdflatex temp.tex')
    shutil.move('temp.pdf','%s.pdf' % filebase)
    for i in ['temp.tex','temp.aux','temp.log']:
        os.remove(i)

def get_results(pwn):
    f = [j(fitdir,pwn,i) for i in ['results_%s_pointlike.yaml' % pwn, 
#                                   'results_%s_extul_point.yaml' % pwn,
                                   'results_%s_gtlike_at_pulsar.yaml' % pwn,
                                   'results_%s_gtlike_point.yaml' % pwn,
                                   'results_%s_gtlike_extended.yaml' % pwn]]
    for i in f: 
        if not os.path.exists(i):
            print '%s does not exist' % i
            return None
    g = [yaml.load(open(i)) for i in f]
    return merge_dict(*g)

def get_sed(pwn,binning,hypothesis):
    filename=j(fitdir,pwn,'seds','sed_gtlike_%s_%s_%s.yaml' % (binning, hypothesis, pwn))
    return yaml.load(open(filename))

def get_pwnlist():
    pwnlist=sorted(yaml.load(open(expandvars('$pwncode/data/pwncat2_data_lande.yaml'))).keys())
    return pwnlist


class BestHypothesis(object):
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

