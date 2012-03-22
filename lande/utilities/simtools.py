#!/usr/bin/env python
from os.path import join, expandvars, exists
from os import makedirs
import collections

from itertools import product

from uw.utilities import keyword_options 
from . lists import flatten

class SimBuilder(object):

    defaults = (
        ('params', None, 'Extra params to pass into the script'),
        ('extra', 'sim', 'Extra params to pass into the script'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, savedir, code, num, **kwargs):
        keyword_options.process(self, kwargs)

        self.savedir = expandvars(savedir)
        self.code = code
        self.num = num

    def build(self):
        if self.params is not None:

            keys = flatten(self.params.keys())
            values = self.params.values()

            for perm in product(*values):
                perm = flatten(perm)

                f = lambda x: x if isinstance(x,str) else '%g' % x

                base = self.extra + '_' + '_'.join('%s_%s' % (f(k),f(v)) for k,v in zip(keys,perm))

                args= ' '.join('--%s=%s' % (f(k),f(v)) for k,v in zip(keys,perm))

                subdir = join(self.savedir, base)

                for i in range(self.num):
                    istr='%0*d' % (5,i)

                    jobdir = join(subdir,istr)
                    if not exists(jobdir): makedirs(jobdir)

                    run = join(jobdir,'run.sh')
                    open(run,'w').write("python %s %g %s" % (self.code, i, args))

            submit_all = join(self.savedir,'submit_all.sh')
            open(submit_all,'w').write("submit_all */*/run.sh $@")

        else:
            for i in range(self.num):
                istr='%0*d' % (5,i)

                jobdir = join(self.savedir,istr)
                if not exists(jobdir): makedirs(jobdir)

                run = join(jobdir,'run.sh')
                open(run,'w').write("python %s %g" % (self.code,i))

            submit_all = join(self.savedir,'submit_all.sh')
            open(submit_all,'w').write("submit_all */run.sh $@")
