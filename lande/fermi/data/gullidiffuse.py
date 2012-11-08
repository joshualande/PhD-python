""" Code to load in Guli's alternate diffuse models, described here:

        https://confluence.slac.stanford.edu/display/SCIGRPS/Diffuse+Models+for+SNR+Cat+Systematics

    Author: Joshua Lande <joshualande@gmail.com> """
from os.path import basename, join, exists

from uw.like.pointspec_helpers import get_diffuse_source
from uw.like.Models import PowerLaw, Constant

def get_gulli_diffuse(dist, halo, TS, version, event_class,
                      diffdir="/afs/slac/g/glast/groups/diffuse/SNRCatalog",
                      fit_index = False,
                      fit_isotropic=False,
                      add_loop=True, fit_loop=False,
                      add_lobes=True, fit_lobes=False,
                      gal_components=None,
                      verbosity=False
               ):
    """ Code to load in Gulli's Galactic and istropic diffuse files for Pointlike analysis 
    
        fit_index: Fit index of diffuse sources (modulated by a power law). Default is to not fit index 
            (and module instead by a constant).
        fit_isotropic: Fit a scaling powerlaw times the isotropic.
        add_loop & add_lobes: Add the loop and lobe feature. This is recommented.
        fit_loop & fit_lobes: fit the loop and lobes feature. This is not recommended.
            See work by Jean Ballet for best practices https://confluence.slac.stanford.edu/x/jJs7Bw

        Note, to load in these diffuse files, you can use the pointlike code:
            
            >>> spectral_analysis = ...

            >>> roi = spectral_analysis.roi(
            ...    diffuse_sources = get_gulli_diffuse(dist='Lorimer', halo=10, TS=150, version='2', ...
            ... )
    """
    assert (dist in ['Lorimer', 'SNR']) and (halo in [4,10]) and (TS in [150, 100000]) and (version in [1,2]) and (event_class in ['clean', 'source'])

    diffuse = []

    def get_gal(component, free=True):
        if version==1:
            gfile_base = '%s_z%s_Ts%s_v%s_mapcube_fixed_' % (dist,halo,TS,version)
        elif version==2:
            gfile_base = '%s_z%s_Ts%s_v%s_mapcube_' % (dist,halo,TS,version)
        filename = join(diffdir, gfile_base + component + '.fits.gz')
        if verbosity: print ' * Loading file %s' % filename

        if not exists(filename): raise Exception("File %s does not exist." % filename)

        gmodel = PowerLaw(norm=1, index=0) if fit_index else Constant()
        ds = get_diffuse_source('MapCubeFunction',filename,gmodel,None,name=basename(filename))
        ds.smodel.free[:]=free
        return ds

    def get_iso(free=True):
        if version==1:
            raise Exception("No isotropic file for v1 of gulli diffuse.")
        elif version==2:
            filename = join(diffdir, '%s_z%s_Ts%s_v%s_%s_isotropic.txt' % (dist,halo,TS,version,event_class))

        if verbosity: print ' * Loading file %s' % filename

        if not exists(filename): raise Exception("File %s does not exist." % filename)

        ds=get_diffuse_source('ConstantValue',None,'FileFunction',filename,name=basename(filename))
        ds.smodel.free[:]=free
        return ds

    if gal_components is None:
        gal_components=["CO_1","CO_2","CO_3","CO_4","HI_1","HI_2","HI_3","HI_4","ICnorm"]

    diffuse = [get_gal(i) for i in gal_components]

    if add_loop:
        diffuse.append(get_gal('LoopInorm', free=fit_loop))
    if add_lobes:
        diffuse.append(get_gal('Lobesnorm', free=fit_lobes))
    
    diffuse.append(get_iso(free=fit_isotropic))

    return diffuse


