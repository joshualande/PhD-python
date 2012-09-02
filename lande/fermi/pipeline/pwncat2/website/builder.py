from skymaps import SkyDir

import numpy as np
import yaml
from StringIO import StringIO
import os
from os.path import join,exists,expandvars
from collections import OrderedDict

import asciitable

from lande.utilities.website import t2t

from lande.pipeline.pwncat2.interp.loader import PWNResultsLoader




class TableFormatter(object):

    def __init__(self, loader):
        self.loader = loader

    def format(self, pwnlist):

        flux_name=r'F(0.1-100)'
        gamma_name=r'Gamma'

        table = OrderedDict()
        for k in ['PSR', 
                  'TS_at_pulsar_ptlike', 'TS_loc_ptlike', 'TS_ext_ptlike', 'TS_cutoff_ptlike', 'TS_var_ptlike', 'disp',
                  'TS_at_pulsar_gtlike', 'TS_loc_gtlike', 'TS_ext_gtlike', 'TS_cutoff_gtlike', 'TS_var_gtlike',
                 ]:
            table[k] = ['None']*len(pwnlist)

        for i,pwn in enumerate(pwnlist):

            table['PSR'][i]='[%s %s.html]' % (pwn,pwn)

            results = self.loader.get_results(pwn, require_all_exists=False)
            if results is None: continue

            pt_at_pulsar=results['at_pulsar']['pointlike']
            pt_point=results['point']['pointlike']
            pt_extended=results['extended']['pointlike']


            bold = lambda text, doit=True: '**%s**' % text if doit else text

            ts_at_pulsar=pt_at_pulsar['TS']
            ts_point = pt_point['TS']
            ts_gauss = pt_extended['TS']
            ts_ext = ts_gauss - ts_point

            table['TS_at_pulsar_ptlike'][i] = bold('%.1f' % ts_at_pulsar, ts_at_pulsar>25)


            ts_loc = ts_point - ts_at_pulsar
            table['TS_loc_ptlike'][i] = bold('%.1f' % (ts_loc), ts_point>25)
                
            table['TS_ext_ptlike'][i] = bold('%.1f' % ts_ext, ts_point > 25 and ts_ext > 16)



            displacement = np.degrees(SkyDir(*pt_point['position']['equ']).difference(SkyDir(*pt_at_pulsar['position']['equ'])))
            table['disp'][i] = '%.2f' % displacement

            try:
                ts_cutoff = pt_point['test_cutoff']['TS_cutoff']
                table['TS_cutoff_ptlike'][i] = bold('%.1f' % ts_cutoff, ts_cutoff > 16)
            except:
                pass

            try:
                ts_var_ptlike = results['point']['variability']['TS_var']['pointlike']
                table['TS_var_ptlike'][i] = '%.1f' % max(ts_var_ptlike,0)
            except:
                pass

            if results['at_pulsar'].has_key('gtlike'):
                gt_at_pulsar=results['at_pulsar']['gtlike']

                ts_at_pulsar=gt_at_pulsar['TS']['reoptimize']

                table['TS_at_pulsar_gtlike'][i] = bold('%.1f' % ts_at_pulsar, ts_at_pulsar>25)


                if results['point'].has_key('gtlike'):
                    gt_point=results['point']['gtlike']

                    ts_point = gt_point['TS']['reoptimize']
                    ts_loc = ts_point - ts_at_pulsar

                    table['TS_loc_gtlike'][i] = bold('%.1f' % (ts_loc), ts_point>25)

                    if results['extended'].has_key('gtlike'):
                        gt_extended=results['extended']['gtlike']

                        ts_gauss = gt_extended['TS']['reoptimize']
                        ts_ext = ts_gauss - ts_point
                        table['TS_ext_gtlike'][i] = bold('%.1f' % ts_ext, ts_point > 25 and ts_ext > 16)

                    try:
                        ts_cutoff = gt_point['test_cutoff']['TS_cutoff']
                        table['TS_cutoff_gtlike'][i] = bold('%.1f' % ts_cutoff, ts_cutoff > 16)
                    except:
                        pass

            try:
                ts_var_gtlike = results['point']['variability']['TS_var']['gtlike']
                table['TS_var_gtlike'][i] = '%.1f' % max(ts_var_gtlike,0)
            except:
                pass

        return self.get_t2t_table(table)

    @staticmethod
    def get_t2t_table(table, **kwargs):

        outtable=StringIO()

        asciitable.write(table, outtable, 
                         Writer=asciitable.FixedWidth,
                         names=table.keys(),
                         **kwargs
                        )
        t=outtable.getvalue()

        # this is required by t2t for tables
        # see for example: http://txt2tags.org/markup.html
        t='||' + t[2:]
        return t



class WebsiteBuilder(object):
    def __init__(self,
                 pwndata,
                 fitdir,
                 webdir):

        self.fitdir=expandvars(fitdir)
        self.webdir=expandvars(webdir)

        self.relpath=os.path.relpath(self.fitdir, self.webdir)

        if not os.path.exists(self.webdir): os.makedirs(self.webdir)

        self.loader = PWNResultsLoader(pwndata, self.fitdir)

        self.pwnlist=self.loader.get_pwnlist()

        self.formatter = TableFormatter(self.loader)

    def build(self):
        self.build_main_website()
        self.build_all_pages()

    def build_all_pages(self):
        for pwn in self.pwnlist: 
            self.build_each_page(pwn)

    def build_main_website(self):

        index_t2t = []
        index_t2t.append('PWNCatalog+\n\n')
        t=self.formatter.format(self.pwnlist)
        index_t2t.append(str(t))
        t2t(index_t2t, join(self.webdir,'index.t2t'))

    def build_each_page(self,pwn):
        index_t2t = []
        index_t2t.append(pwn+'\n\n')
        index_t2t.append('([back index.html])')
        t=self.formatter.format([pwn])
        index_t2t.append(str(t))
        index_t2t.append('')
        index_t2t.append('[Analysis Folder %s/%s]\n' % (self.relpath,pwn))
        index_t2t.append('[log (pointlike) %s/%s/log_run_%s.txt]\n' % (self.relpath,pwn,pwn))

        index_t2t.append('[results (pointlike) %s/%s/results_%s_pointlike.yaml]\n' % (self.relpath,pwn,pwn))

        get_img_table = lambda *args: index_t2t.append('|| ' + ' | '.join(['[%s/%s/%s]' % (self.relpath,pwn,i) for i in args]) + ' |\n\n')
        get_sed_table = lambda *args: index_t2t.append('|| ' + ' | '.join(['[%s/%s/%s]' % (self.relpath,pwn,i) for i in args]) + ' |\n\n')

        title = lambda i: index_t2t.append('\n\n== %s ==' % i)

        title('Phase Info')
        get_img_table('plots/phaseogram_%s.png' % (pwn),'plots/phase_vs_time_%s.png' % (pwn))

        all = ['at_pulsar', 'point', 'extended']

        title('Big Residual TS map')
        get_img_table(*['plots/tsmap_residual_%s_%s_10deg.png' % (i,pwn) for i in all])

        index_t2t.append('[tsmap_residual_%s_%s_10deg.fits %s/%s/data/tsmap_residual_%s_%s_10deg.fits]' % (pwn,'at_pulsar',self.relpath,pwn,'at_pulsar',pwn))

        title('SED')
        get_sed_table(*['seds/sed_gtlike_4bpd_%s_%s.png' % (i,pwn) for i in all])

        title('gtlike Cutoff test')
        get_sed_table(*['plots/test_cutoff_%s_%s.png' % (i,pwn) for i in ['at_pulsar', 'point']])

                
        title('Source TS Maps')
        get_img_table(*['plots/tsmap_source_%s_%s_5deg.png' % (i,pwn) for i in all])
        get_img_table(*['plots/band_tsmap_source_%s_%s_5deg.png' % (i,pwn) for i in all])

        title('Residual TS Maps')
        get_img_table(*['plots/tsmap_residual_%s_%s_5deg.png' % (i,pwn) for i in all])
        get_img_table(*['plots/band_tsmap_residual_%s_%s_5deg.png' % (i,pwn) for i in all])

        title('New Source TS Maps')
        get_img_table(*['plots/tsmap_newsrc_%s_%s_5deg.png' % (i,pwn) for i in all])
        get_img_table(*['plots/band_tsmap_newsrc_%s_%s_5deg.png' % (i,pwn) for i in all])

        title('at_pulsar Smoothed Counts Diffuse Subtracted (0.1)')
        get_img_table(*['plots/sources_%s_%s_5deg_0.1deg.png' % (i,pwn) for i in all])
        get_img_table(*['plots/band_sources_%s_%s_5deg_0.1deg.png' % (i,pwn) for i in all])

        title('Smoothed Counts BG Source Subtracted (0.1)')
        get_img_table(*['plots/source_%s_%s_5deg_0.1deg.png' % (i,pwn) for i in all])



        title('Band Smoothed Counts BG Source Subtracted (0.1)')
        get_img_table(*['plots/band_source_%s_%s_5deg_0.1deg.png' % (i,pwn) for i in all])

        title('gtlike SED')
        get_sed_table(*['seds/sed_gtlike_1bpd_%s_%s.png' % (i,pwn) for i in all])
        get_sed_table(*['seds/sed_gtlike_2bpd_%s_%s.png' % (i,pwn) for i in all])



        title('Pointlike SEDs')
        get_sed_table(*['seds/sed_pointlike_%s_%s.png' % (i,pwn) for i in all])

        title('Extra: Source TS Maps (10 deg)')
        get_img_table(*['plots/tsmap_source_%s_%s_10deg.png' % (i,pwn) for i in all])
        get_img_table(*['plots/band_tsmap_source_%s_%s_10deg.png' % (i,pwn) for i in all])

        title('Extra: Residual TS Maps')
        get_img_table(*['plots/tsmap_residual_%s_%s_10deg.png' % (i,pwn) for i in all])
        get_img_table(*['plots/band_tsmap_residual_%s_%s_10deg.png' % (i,pwn) for i in all])

        title('Extra: New Source TS Maps (10 deg)')
        get_img_table(*['plots/tsmap_newsrc_%s_%s_10deg.png' % (i,pwn) for i in all])
        get_img_table(*['plots/band_tsmap_newsrc_%s_%s_10deg.png' % (i,pwn) for i in all])


        title('Extra: Smoothed Counts (0.25deg)')
        get_img_table(*['plots/source_0.25_%s_%s_5deg.png' % (i,pwn) for i in all])

        title('Extra: Smoothed Counts (0.25deg)')
        get_img_table(*['plots/sources_0.25_%s_%s_5deg.png' % (i,pwn) for i in all])


        title('Extra: Band Smoothed Counts (0.25)')
        get_img_table(*['plots/band_source_%s_%s_5deg_0.25deg.png' % (i,pwn) for i in all])
        get_img_table(*['plots/band_sources_%s_%s_5deg_0.25deg.png' % (i,pwn) for i in all])

        title('Counts (0.1)')
        get_img_table(*['plots/counts_residual_%s_%s_5deg_0.1deg.png' % (i,pwn) for i in all])

        get_img_table(*['plots/counts_source_%s_%s_5deg_0.1deg.png' % (i,pwn) for i in all])


        title('Extra: Counts (0.25)')
        get_img_table(*['plots/counts_source_%s_%s_5deg_0.25deg.png' % (i,pwn) for i in all])
        get_img_table(*['plots/counts_residual_%s_%s_5deg_0.25deg.png' % (i,pwn) for i in all])


        t2t(index_t2t, join(self.webdir,'%s.t2t' % pwn))
