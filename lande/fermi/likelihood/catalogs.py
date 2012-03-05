from os.path import join

def get_2fgl():
    catalog_basedir = "/afs/slac/g/glast/groups/catalog/P7_V4_SOURCE"
    d = dict(ft2=join(catalog_basedir,"ft2_2years.fits"),
             ltcube=join(catalog_basedir,"ltcube_24m_pass7.4_source_z100_t90_cl0.fits")
    )
    return d
