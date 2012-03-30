#!/usr/bin/env python
from os.path import join, expandvars, exists, splitext
from os import makedirs
from collections import defaultdict
from itertools import product

import numpy as np
import yaml

from uw.utilities import keyword_options 

from . lists import flatten
from . files import locate
from . save import savedict, loaddict

class SimBuilder(object):

    defaults = (
        ('params', {}, 'Extra params to pass into the script'),
        ('front', 'sim'),
        ('extra', ''),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, savedir, code, num, **kwargs):
        keyword_options.process(self, kwargs)

        self.savedir = expandvars(savedir)
        self.code = code
        self.num = num

    def build(self):

        ext=splitext(self.code)[-1]
        if ext == '.py':
            program='python'
        elif ext == '.sh':
            program='sh'
        else:
            raise Exception("...")

        for k,v in self.params.items():
            if not isinstance(v,list) and not isinstance(v,tuple):
                self.params[k] = [v]

        keys = flatten(self.params.keys())
        values = self.params.values()

        no_multiples = not np.any([len(v)>1 for v in values])

        for perm in product(*values):
            perm = flatten(perm)

            f = lambda x: x if isinstance(x,str) else '%g' % x

            if no_multiples:
                base = '.'
            else:
                base = self.front + '_' + '_'.join('%s_%s' % (f(k),f(v)) for k,v in zip(keys,perm) if len(self.params[k])>1)

            args= ' '.join('--%s=%s' % (f(k),f(v)) for k,v in zip(keys,perm))

            subdir = join(self.savedir, base)

            for i in range(self.num):
                istr='%0*d' % (5,i)

                jobdir = join(subdir,istr)
                if not exists(jobdir): makedirs(jobdir)

                run = join(jobdir,'run.sh')

                open(run,'w').write("%s %s %g %s %s" % (program,self.code, i, args, self.extra))

        submit_all = join(self.savedir,'submit_all.sh')
                                      
        if no_multiples:
            run='*/run.sh'
        else:
            run='*/*/run.sh'
        open(submit_all,'w').write("submit_all %s $@" % run)

class SimMerger(object):

    defaults = (
        ('filename', '*.yaml', 'Name of files to merge. Can allow wildcard matches.'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, savedir, keys, **kwargs):
        keyword_options.process(self, kwargs)

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
                >>> print SimMerger.traverse(d, ['non_existing'])
                Traceback (most recent call last):
                    ...
                KeyError: 'Unable to traverse structure with key non_existing'

            If one of the keys is a list of keys, this function will traverse the list
            in whatever way it can:

                >>> d = dict(a=dict(b1=dict(d=0)))
                >>> print SimMerger.traverse(d, [ 'a', ['b1', 'b2'], 'd' ])
                0
                >>> print SimMerger.traverse(d, [ 'a', ['b2', 'b3'], 'd' ])
                Traceback (most recent call last):
                    ...
                KeyError: "Unable to traverse structure with any of the keys ['b2', 'b3']"

        """
        k = keys[0]
        if isinstance(k,list):
            found=False
            for i in k:
                try:
                    sub_data=data[i]
                    found=True
                    break
                except:
                    pass
            if not found: 
                raise KeyError("Unable to traverse structure with any of the keys %s" % k)
        else:
            try:
                sub_data=data[k]
            except:
                raise KeyError("Unable to traverse structure with key %s" % k)

        if len(keys) == 1:
            return sub_data
        else:
            return SimMerger.traverse(sub_data, keys[1:])

    def _merge(self):
        self.results = defaultdict(list)

        all_results=locate(self.filename, self.savedir)

        for i,r in enumerate(all_results):
            if i % 10==0: print '%s' % i

            x = loaddict(r)

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
