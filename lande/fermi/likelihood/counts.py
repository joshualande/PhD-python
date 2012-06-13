import numpy as np

def pointlike_observed_counts(roi):
    """ Usage:
            observed_counts = pointlike_observed_counts(roi)
    """
    observed_counts = sum(b.pix_counts.sum() for b in roi.bands)
    return observed_counts


def pointlike_model_counts(roi,which):
    """ Usage:
            model_counts observed_counts = pointlike_model_counts(roi,'source_name')
    """
    manager,index=roi.mapper(which)
    assert manager in [roi.psm,roi.dsm]
    if manager == roi.psm:
        model_counts=sum(b.ps_counts[index] for b in roi.bands)
    elif manager == roi.dsm:
        model_counts=sum(b.bg_counts[index] for b in roi.bands)
    return model_counts

def pointlike_total_model_counts(roi):
    model_counts=sum(b.bg_all_counts+b.ps_all_counts for b in roi.bands)
    return model_counts



# Kluge for now
# TODO, write gtlike equivalent code
observed_counts=pointlike_observed_counts
model_counts=pointlike_model_counts
total_model_counts=pointlike_total_model_counts
