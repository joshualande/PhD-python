#!/usr/bin/env python
from os.path import join, expandvars, exists
from os import makedirs
from collections import defaultdict
from itertools import product

import yaml

from uw.utilities import keyword_options 

from . lists import flatten
from . files import locate
from . save import savedict

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

            for k,v in self.params.items():
                if not isinstance(v,list) and not isinstance(v,tuple):
                    self.params[k] = [v]

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


class SimMerger(object):

    def __init__(self, savedir, keys):
        self.savedir = expandvars(savedir)
        self.keys = keys

        self._merge()

    @staticmethod
    def traverse(data, keys):
        """ Pull out of a complicated data strucutre
            a value where keys is a list of the 
            keys to follow.

            Example:

                >>> d = dict(a=dict(b=dict(c=[0,dict(a='treasure')])))
                >>> print SimMerger.traverse(d, ['a', 'b', 'c', 1, 'a'])
                treasure
        """
        if len(keys) == 1:
            return data[keys[0]]
        else:
            return SimMerger.traverse(data[keys[0]], keys[1:])

    def _merge(self):
        self.results = defaultdict(list)

        all_results=locate('results_*.yaml', self.savedir)

        for i,r in enumerate(all_results):
            if i % 10==0: print '%s' % i

            x = yaml.load(open(r))

            if x is None: continue

            for k,v in self.keys.items():

                # if each reults file contains a straight array, assume
                # that each of the entries in the results file needs to 
                # be parsed identically. I hope this is a reasonable assumption.
                if isinstance(x,list):
                    for i in x: self.results[k].append(SimMerger.traverse(i,v))
                else:
                    self.results[k].append(SimMerger.traverse(x,v))
                    

    def save(self, filename):
        savedict(filename, self.results)
        

if __name__ == "__main__":
    import doctest
    doctest.testmod()
