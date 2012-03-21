#!/usr/bin/env python
from os.path import join, expandvars, exists
from os import makedirs

from uw.utilities import keyword_options 
from itertools import product

class SimBuilder(object):

    defaults = (
        ('params', None, 'Extra params to pass into the script'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, savedir, code, num, **kwargs):
        keyword_options.process(self, kwargs)

        self.savedir = expandvars(savedir)
        self.code = code
        self.num = num

    def build(self):

        if params is None:

            keys = params.keys()
            values = params.values()

            for perm in product(values):

                base = '_'.join(p if not isinstance(p, tuple) else '_'.join(p) for p in perm)
                extra = ' '.join('--%s' for k,v in zip(keys,values))

                print 'base',base
                print 'extra',extra

                for i in range(self.num):
                    istr='%0*d' % (5,i)

                    jobdir = join(self.savedir,istr)
                    if not exists(jobdir): makedirs(jobdir)

                    run = join(jobdir,'run.sh')
                    open(run,'w').write("python %s %g" % (self.code,i))

            submit_all = join(self.savedir,'submit_all.sh')
            open(submit_all,'w').write("submit_all */run.sh $@")

        else:
            for i in range(self.num):
                istr='%0*d' % (5,i)

                jobdir = join(self.savedir,istr)
                if not exists(jobdir): makedirs(jobdir)

                run = join(jobdir,'run.sh')
                open(run,'w').write("python %s %g" % (self.code,i))

            submit_all = join(self.savedir,'submit_all.sh')
            open(submit_all,'w').write("submit_all */run.sh $@")
