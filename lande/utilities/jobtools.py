#!/usr/bin/env python
from os.path import join, expandvars, exists, splitext
from os import makedirs
from collections import defaultdict
from itertools import product

import numpy as np
import yaml

from uw.utilities import keyword_options 

from . lists import flatten, islist, duplicates
from . files import locate
from . save import savedict, loaddict

class JobBuilder(object):
    """ Example usage:

            # version 27 - more options
            params=dict(edisp=[True,False], simbg=[True,False], emin=[1e2,1e3], emax=1e5, flux=1e-5,
                        index=[2,2.66], cuts=[True,False], zenithcut=[100,180], savedata=True)
            params['time','phibins']=[['2fgl',0], ['my2fgl',0], ['my2fgl',9], ['2years',0], ['2years',9]]
            b = JobBuilder(
                savedir='$w44simdata/v27',
                code='$w44simcode/simspec.py',
                num=1,
                params=params,
                )
            b.build()
    """
    defaults = (
        ('params', {}, 'Extra params to pass into the script'),
        ('front', 'job'),
        ('extra', ''),
        ('num', None),
        ('short_folder_names', False, """ By default, folder names are long: emin_1e2_emax_1e5. 
                                         short_folder_names will make folder names shorter: 1e2_1e5. """)
    )

    @keyword_options.decorate(defaults)
    def __init__(self, savedir, code, **kwargs):
        keyword_options.process(self, kwargs)

        self.savedir = expandvars(savedir)
        self.code = code

        print 'Putting jobs in %s' % self.savedir

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
        if len(duplicates(keys)) > 0:
            raise Exception("Duplicate keys for params: %s" % (', '.join(duplicates(keys))))
        values = self.params.values()
        num = flatten([len(v) if not islist(k) else [len(v)]*len(k) for k,v in self.params.items()])

        no_multiples = not np.any([n>1 for n in num])

        for perm in product(*values):
            perm = flatten(perm)

            def f(x):
                if isinstance(x,str):
                    return x
                elif isinstance(x,bool):
                    return str(x)
                else:
                    return '%g' % x

            if no_multiples:
                base = '.'
            else:
                if self.short_folder_names:
                    base = '_'.join(f(v) for k,n,v in zip(keys,num,perm) if n>1)
                else:
                    base = self.front + '_' + '_'.join('%s_%s' % (f(k),f(v)) for k,n,v in zip(keys,num,perm) if n>1)

            args = []
            for k,v in zip(keys,perm):
                if v is True:
                    args.append('\\\n    --%s' % f(k))
                elif v is False:
                    pass # no flag
                else:
                    args.append('\\\n    --%s=%s' % (f(k),f(v)))
            args = ' '.join(args)

            subdir = join(self.savedir, base)
            if not exists(subdir): makedirs(subdir)

            if self.num == None:

                run = join(subdir,'run.sh')
                open(run,'w').write("%s %s %s %s" % (program, self.code, args, self.extra))

            else:
                for i in range(self.num):
                    istr='%0*d' % (5,i)

                    if self.num > 1:
                        jobdir = join(subdir,istr)
                    else:
                        jobdir = subdir
                    if not exists(jobdir): makedirs(jobdir)

                    run = join(jobdir,'run.sh')

                    open(run,'w').write("%s %s %g %s %s" % (program, self.code, i, args, self.extra))

        submit_all = join(self.savedir,'submit_all.sh')
                                      
        if no_multiples or self.num in [1,None]:
            run='*/run.sh'
        else:
            run='*/*/run.sh'
        open(submit_all,'w').write("submit_all %s $@" % run)

class JobMerger(object):

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
                >>> print JobMerger.traverse(d, ['a', 'b', 'c', 1, 'a'])
                treasure
                >>> print JobMerger.traverse(d, ['non_existing'])
                Traceback (most recent call last):
                    ...
                KeyError: 'Unable to traverse structure with key non_existing'

            If one of the keys is a list of keys, this function will traverse the list
            in whatever way it can:

                >>> d = dict(a=dict(b1=dict(d=0)))
                >>> print JobMerger.traverse(d, [ 'a', ['b1', 'b2'], 'd' ])
                0
                >>> print JobMerger.traverse(d, [ 'a', ['b2', 'b3'], 'd' ])
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
            return JobMerger.traverse(sub_data, keys[1:])

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
                    for i in x: self.results[k].append(JobMerger.traverse(i,v))
                else:
                    self.results[k].append(JobMerger.traverse(x,v))
                    

    def save(self, filename):
        savedict(filename, self.results)
        

if __name__ == "__main__":
    import doctest
    doctest.testmod()
