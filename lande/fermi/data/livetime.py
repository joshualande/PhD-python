import os
import math
from os.path import exists, join
from textwrap import dedent
from tempfile import NamedTemporaryFile

import numpy as np
import pyfits

from GtApp import GtApp
import skymaps

from uw.like.roi_monte_carlo import MonteCarlo

from lande.utilities.jobtools import JobBuilder
from lande.utilities.shell import format_command

def fix_pointlike_ltcube(ltcube):
    """ Modify a pointlike ltcube so that it can be used in gtlike. """

    f = pyfits.open(ltcube, mode='update')

    for exposure in ['EXPOSURE', 'WEIGHTED_EXPOSURE']:
        if exposure in [i.name for i in f]:
            h=f[exposure].header

            if 'NDSKEYS' not in h.keys():
                h.update('DSTYP1', 'TIME')
                h.update('DSUNI1', 's')
                h.update('DSVAL1', 'TABLE')
                h.update('DSREF1', ':GTI')
                h.update('NDSKEYS', 1)
    f.flush()

def gtltcube(**kwargs):
    l=GtApp('gtltcube', 'Likelihood')
    return l.run(**kwargs)

def gtltsum(**kwargs):
    l=GtApp('gtltsum')
    return l.run(**kwargs)


def batch_gtltcube(evfile, scfile, outfile, savedir, njobs=100, **kwargs):
    """ Create a folder savedir where a livetime cube can be computd
        in parallel. """
    if not exists(savedir): os.makedirs(savedir)

    first,last=MonteCarlo.get_time_from_ft2(scfile)

    times = np.linspace(first, last, njobs + 1)

    l = max(len('%d' % t) for t in times)

    tmins = times[:-1]
    tmaxs = times[1:]

    for tmin,tmax in zip(tmins,tmaxs):
        tmin_fmt='%0*d' % (l,tmin)
        tmax_fmt='%0*d' % (l,tmax)
        subdir = join(savedir,'times_%s_%s' % (tmin_fmt,tmax_fmt))

        if not exists(subdir): os.makedirs(subdir)

        cut_evfile=os.path.basename(evfile).replace('.fits','_%s_%s.fits')
        gtselect_command = format_command('gtselect',
                                          infile=evfile,
                                          outfile=cut_evfile,
                                          ra=0, dec=0, rad=180,
                                          tmin=tmin, tmax=tmax,
                                          emin=1, emax=1e6,
                                          zmax=180)

        ltcube_kwargs = dict(
            evfile=cut_evfile,
            outfile='ltcube_%s_%s.fits' % (tmin_fmt,tmax_fmt),
            scfile=scfile,
            tmin=tmin,
            tmax=tmax,)
        ltcube_kwargs.update(kwargs)

        ltcube_command = format_command('gtltcube',**ltcube_kwargs)

        jobfile = join(subdir,'run.sh')
        open(jobfile,'w').write('\n\n'.join([gtselect_command,ltcube_command]))

    submit_all = join(savedir,'submit_all.sh')
    open(submit_all,'w').write('submit_all */run.sh $@')

    merge = join(savedir,'merge.py')
    open(merge,'w').write(dedent("""\
        from lande.fermi.data.livetime import recursive_gtltsum
        from glob import iglob
        recursive_gtltsum(iglob("*/ltcube*.fits"),"%s")""" % outfile))


def recursive_gtltsum(infiles, outfile):
    """ Run gtltsum recursively to sum all livetime cubes. """
    infiles = list(infiles) # in case generator
    if len(infiles) < 2: raise Exception('Must sum >= 2 livetime cubes')

    temp = NamedTemporaryFile(suffix='.fits',delete=False)
    tempfile = temp.name

    recursive_gtltsum.i=1

    def sum_ltcube(infile1,infile2):
        print 'Merging file %s and %s (%s/%s)' % (infile1, infile2,recursive_gtltsum.i,len(infiles))
        recursive_gtltsum.i+=1
        if infile1 == outfile:
            gtltsum(infile1=infile1, infile2=infile2, outfile=tempfile)
            return tempfile
        else:
            gtltsum(infile1=infile1, infile2=infile2, outfile=outfile)
            return outfile

    accum=reduce(sum_ltcube, infiles)
    if accum == tempfile:
        os.move(tempfile, outfile)
    else:
        os.remove(tempfile)


def pointlike_ltcube(evfile,scfile,outfile,dcostheta,binsz, zmax, cone_angle,dir,quiet=False):
    """ Run the pointlike ltcube code. """

    gti = skymaps.Gti(evfile)

    for weighted,clobber in [[False,True],[True,False]]:

        lt = skymaps.LivetimeCube(
            cone_angle = cone_angle,
            dir        = dir,
            zcut       = math.cos(math.radians(zmax)),
            pixelsize  = binsz,
            quiet      = quiet,
            weighted   = weighted)

        lt.load(scfile,gti)

        extension = 'WEIGHTED_EXPOSURE' if weighted else 'EXPOSURE'
        lt.write(outfile,extension,clobber)

    fix_pointlike_ltcube(outfile)
