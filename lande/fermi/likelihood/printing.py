import numpy as np

from pyLikelihood import StringVector

from . save import get_all_names, get_full_energy_range, get_roi_dir, get_sources, get_background, get_skydir, logLikelihood

from lande.utilities.tools import OrderedDefaultDict
from lande.utilities.table import fixed_width_table


def gtlike_summary(like, sdir=None, galactic=True, maxdist=5, sep='-'*90, indent='    '):
    """ Print out a summary of information in a gtlike ROI. """

    emin,emax=get_full_energy_range(like)

    if sdir is None:
        sdir = get_roi_dir(like)
    
    sources = np.asarray(get_sources(like))
    diffs = np.degrees([get_skydir(like, name).difference(sdir) for name in sources])

    # cut on distance
    sources, diffs = sources[diffs <= maxdist], diffs[diffs <= maxdist]

    # sort
    sources = sources[np.argsort(diffs)]

    background = get_background(like)
    
    def format(names, do_position, do_ts, do_flux, do_norm):
        if len(names) < 1:
            return ''

        dat = OrderedDefaultDict(list)
        for name in names:
            dat['name'].append(name)

            source = like.logLike.getSource(name)
            spectrum = source.spectrum()

            model = spectrum.genericName()
            dat['model'].append(model)

            if do_position:
                skydir = get_skydir(like, name)
                dist = np.degrees(skydir.difference(sdir))

                dat['dist'].append('%.1f' % dist)
                dat['l' if galactic else 'ra'].append('%.2f' % (skydir.l() if galactic else skydir.ra()))
                dat['b' if galactic else 'dec'].append('%.2f' % (skydir.b() if galactic else skydir.dec()))

            if do_ts:
                ts = like.Ts(name, reoptimize=False)
                dat['TS'].append('%.1f' % ts)

            if do_flux:
                flux = like.flux(name,emin,emax)
                dat['flux8'].append('%.2f' % (flux/1e-8))


            parNames = StringVector()
            spectrum.getParamNames(parNames)

            if do_norm:
                norm_name = parNames[0]
                norm = spectrum.getParam(norm_name).getTrueValue()
                dat['norm'].append('%.2f' % norm)

            if len(parNames) > 1:
                index_name = parNames[1]
                index = spectrum.getParam(index_name).getTrueValue()
                dat['index'].append('%.2f' % index)
            else:
                dat['index'].append('')

            free = StringVector()
            spectrum.getFreeParamNames(free)
            any_free = len(free) > 0

            dat['Free'].append('True' if any_free else 'False')

        return fixed_width_table(dat, indent=indent)

    ll=logLikelihood(like)
    return '\n'.join([
        sep,
        'Point + Extended Sources',
        sep,
        format(sources, do_position=True, do_ts=True, do_flux=True, do_norm=False),
        sep,
        'Background sources',
        sep,
        format(background, do_position=False, do_ts=False, do_flux=False, do_norm=True),
        sep,
        'logLikelihood=%s' % ll,
        sep,])

summary=gtlike_summary
