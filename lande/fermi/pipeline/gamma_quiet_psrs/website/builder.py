from os.path import join, expandvars, relpath
import traceback

import yaml

from lande.utilities.website import t2t
from lande.utilities.tools import OrderedDefaultDict
from lande.utilities.table import t2t_table
from lande.fermi.pipeline.gamma_quiet_psrs.interp.loader import ResultsLoader


class WebsiteBuilder(object):

    def __init__(self, gamma_quiet_psrs_data, analysisdir, webdir):
        self.loader=ResultsLoader(gamma_quiet_psrs_data=gamma_quiet_psrs_data,
                                  analysisdir=analysisdir)
        self.webdir=expandvars(webdir)
        self.relpath=relpath(self.loader.analysisdir, self.webdir)
        
    def summary_table(self, psr_list):
        table = OrderedDefaultDict(lambda:['None']*len(psr_list))
        for i,psr in enumerate(psr_list):
            table['names'][i]='[%s %s.html]' % (psr,psr)

            results = self.loader.get_results(psr, require_all_exists=False)

            try:
                table['TS_at_pulsar (p)'][i] = '%.1f' % results['at_pulsar']['pointlike']['TS']['noquick']
            except:
                print traceback.format_exc()

            try:
                table['TS_point (p)'][i] = '%.1f' % results['point']['pointlike']['TS']['noquick']
            except:
                print traceback.format_exc()

            try:
                table['TS_extended (p)'][i] = '%.1f' % results['extended']['pointlike']['TS']['noquick']
            except:
                print traceback.format_exc()

            try:
                table['TS_at_pulsar (g)'][i] = '%.1f' % results['at_pulsar']['gtlike']['TS']['reoptimize']
            except:
                print traceback.format_exc()

            try:
                table['TS_point (g)'][i] = '%.1f' % results['point']['gtlike']['TS']['reoptimize']
            except:
                print traceback.format_exc()

            try:
                table['TS_extended (g)'][i] = '%.1f' % results['extended']['gtlike']['TS']['reoptimize']
            except:
                print traceback.format_exc()

        print 'table',table

        return t2t_table(table)

    def build_main(self):
        psr_list = self.loader.get_psr_list()


        filename=join(self.webdir,'index.t2t')

        website=[
            'Gamma Quiet PSRs',
            '',
            '',
            self.summary_table(psr_list)
        ]

        t2t(website,filename)

    def build_each_psr(self,psr):
        filename=join(self.webdir,'%s.t2t' % psr)

        get_sed_table = lambda *args: '|| ' + ' | '.join(['[%s/fits/%s/seds/%s]' % (self.relpath,psr,i) for i in args]) + ' |\n\n'

        get_plot_table = lambda *args: '|| ' + ' | '.join(['[%s/fits/%s/plots/%s]' % (self.relpath,psr,i) for i in args]) + ' |\n\n'

        website = [
            '%s+' % psr,
            '',
            '',
            self.summary_table([psr]),
            '',
            'SEDs',
            get_sed_table(*['sed_pointlike_4bpd_%s_%s.png' % (i,psr) for i in self.loader.all_hypotheses]),
            get_sed_table(*['sed_gtlike_2bpd_%s_%s.png' % (i,psr) for i in self.loader.all_hypotheses]),
            get_sed_table(*['sed_gtlike_4bpd_%s_%s.png' % (i,psr) for i in self.loader.all_hypotheses]),
            'TS Maps',
            get_plot_table(*['tsmap_source_%s_%s_5.0deg.png' % (i,psr) for i in self.loader.all_hypotheses]),
            get_plot_table(*['tsmap_source_%s_%s_10.0deg.png' % (i,psr) for i in self.loader.all_hypotheses]),
            get_plot_table(*['tsmap_residual_%s_%s_5.0deg.png' % (i,psr) for i in self.loader.all_hypotheses]),
            get_plot_table(*['tsmap_residual_%s_%s_10.0deg.png' % (i,psr) for i in self.loader.all_hypotheses]),
            'Plots',
            get_plot_table(*['source_%s_%s_5.0deg.png' % (i,psr) for i in self.loader.all_hypotheses]),
            get_plot_table(*['source_%s_%s_10.0deg.png' % (i,psr) for i in self.loader.all_hypotheses]),
            get_plot_table(*['sources_%s_%s_5.0deg.png' % (i,psr) for i in self.loader.all_hypotheses]),
            get_plot_table(*['sources_%s_%s_10.0deg.png' % (i,psr) for i in self.loader.all_hypotheses]),
            '```',
            yaml.dump(self.loader.get_results(psr,require_all_exists=False)),
            '```',
        ]

        t2t(website,filename)

    def build_all_psrs(self):
        psr_list = self.loader.get_psr_list()
        for psr in psr_list:
            self.build_each_psr(psr)
