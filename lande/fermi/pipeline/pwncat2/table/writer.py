from os.path import expandvars
from textwrap import dedent
import shutil
from lande.utilities.table import latex_table,confluence_table

class TableWriter(object):
    def __init__(self, table, savedir, filebase):
        self.savedir = expandvars(savedir)
        self.table = table
        self.filebase= filebase

    def write_confluence(self, **kwargs):

        t = confluence_table(self.table, **kwargs)
        os.chdir(self.savedir)
        open('%s.confluence' % self.filebase,'w').write(t)

    def write_latex(self, preamble='',**kwargs):

        t = latex_table(self.table, **kwargs)

        lines = t.split('\n')
        if lines[-1] == '': 
            lines=lines[0:-1]

        header = lines[0]
        footer= lines[-1]

        t = '\n'.join(lines[1:-1])

        os.chdir(self.savedir)

        open('%s.tex' % self.filebase,'w').write(t)

        open('temp.tex','w').write(dedent(r"""
            \documentclass{aastex}
            \usepackage{amsmath}

            \input{$pwnpaper/style/style.tex}

            \begin{document}
            %s
            %s
            \input{%s}
            %s
            \end{document}""" % (header,preamble,self.filebase,footer)))

        os.system('pdflatex temp.tex')
        shutil.move('temp.pdf','%s.pdf' % self.filebase)
        for i in ['temp.tex','temp.aux','temp.log']:
            os.remove(i)


