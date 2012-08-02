import copy
import numpy as np

def sort_sources(sources, skydir):
    """ Sort soruces based upon distance from a skydir. """
    sources = copy.copy(sources)
    sources.sort(key=lambda src:src.skydir.difference(skydir))
    return sources

def freeze_far_away(roi, skydir, max_free):
    """ Freeze sources far away from the ROI. keep at most max_free sources free. """
    frozen = dict()

    sources=roi.get_sources()
    sorted = sort_sources(sources, skydir)

    i=0
    for source in sorted:
        if np.any(source.model.free==True):
            i+=1
            if i > max_free:
                frozen[source.name]=source.model.free.copy()
                roi.modify(which=source, free=False)
    return frozen

def unfreeze_far_away(roi, frozen):
    for name,free in frozen.items():
        roi.modify(which=name, free=free)
