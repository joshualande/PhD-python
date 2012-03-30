""" Store paths to important catalog data. 
"""
from glob import glob
dict2fgl = dict(ft2="/afs/slac/g/glast/groups/catalog/P7_V4_SOURCE/ft2_2years.fits",
                ltcube="/afs/slac/g/glast/groups/catalog/P7_V4_SOURCE/ltcube_24m_pass7.4_source_z100_t90_cl0.fits",
                ft1=glob("/afs/slac/g/glast/groups/catalog/P7_V4_SOURCE/pass7.3_pre_source_merit_*_pass7.4_source_z100_t90_cl0.fits"),
               )
