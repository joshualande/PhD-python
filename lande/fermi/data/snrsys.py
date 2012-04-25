""" Code to load in Guli's alternate diffuse models, described here:

        https://confluence.slac.stanford.edu/display/SCIGRPS/Diffuse+Models+for+SNR+Cat+Systematics

    Author: Joshua Lande <joshualande@gmail.com> """
from os.path import basename, join as j, exists

from uw.like.pointspec_helpers import get_diffuse_source
from uw.like.Models import PowerLaw, Constant

def get_multiple_diffuse(dist, halo, TS, version,
                         diffdir="/afs/slac/g/glast/groups/diffuse/SNRCatalog",
                         ifile = '/afs/slac/g/glast/groups/diffuse/rings/2year/isotrop_2year_P76_source_v0.txt',
                         fit_index = False,
                         fit_isotropic = False,
                         add_loop=True, add_lobes=True,
               ):
    """ Code to load in Gulli's diffuse files for Pointlike analysis 
        and also one isotropic file. 
    
        fit_index: Fit index of diffuse sources. Default is to not fit index.
        fit_isotropic: Fit isotropic component. Default is to not fit isotropic

        Note, to load in these diffuse files, you can use the pointlike code:
            
            >>> spectral_analysis = ...

            >>> roi = spectral_analysis.roi(
            ...    diffuse_sources = get_multiple_diffuse(dist='Lorimer', halo=10, TS=150, version='2', ...
            ... )
    """
    assert (dist in ['Lorimer', 'SNR']) and (halo in [4,10]) and (TS in [150, 100000]) and (version in [1,2])

    gfile_base = '%s_z%s_Ts%s_v%s_mapcube_fixed_' % (dist,halo,TS, version)

    diffuse = []

    gal_components=["CO_1","CO_2","CO_3","CO_4","HI_1","HI_2","HI_3","HI_4","IC"]
    if add_loop:
        gal_components.append('LoopInorm')
    if add_lobes:
        gal_components.append('Lobesnorm')

    galactic_files = [j(diffdir, gfile_base + i + '.fits.gz') for i in gal_components]

    for file in galactic_files:

        if not exists(file): raise Exception("ERROR: file %s does not exist." % file)

        gmodel = PowerLaw(p=[1,1],index_offset=1) if fit_index else Constant()
        diffuse.append(
            get_diffuse_source('MapCubeFunction',file,gmodel,None,name=basename(file))
        )

    imodel = Constant(free=[fit_isotropic])
    diffuse.append(
        get_diffuse_source('ConstantValue',None,imodel,ifile,name=basename(ifile))
    )

    return diffuse


