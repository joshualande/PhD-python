from lande.fermi.likelihood.counts import model_counts, observed_counts

def is_significant(roi, name, allowed_fraction, verbosity=False):
    oc=observed_counts(roi)
    mc=model_counts(roi,name)
    fraction=float(mc)/oc
    if verbosity: 
        print ' .. Source %s predicts %0.2f%% of total counts' % (name,fraction)
    if fraction < allowed_fraction:
        return False
    else:
        return True

def get_insignificant_diffuse(roi, allowed_fraction, verbosity=False):
    insignificant_list = []
    for source in get_background(roi):
        if is_significant(roi, source, allowed_fraction, verbosity):
            if verbosity: print '... keep source.'
        else:
            insignificant_list.append(source)
    return insignificant_list


def delete_insignificant_diffuse(roi, *args, **kwargs):
    insignificant = get_insignificant_diffuse(roi, *args, **kwargs)
    for source in insignificant:
        roi.del_source(source)

def freeze_insignificant_diffuse(roi, *args, **kwargs):
    """ Freeze insignificant components of the galactic diffuse
        model following the algorithm proposed by Jean Ballet:

            https://confluence.slac.stanford.edu/display/SCIGRPS/gtlike+with+many+diffuse+components
    """
    insignificant = get_insignificant_diffuse(roi, *args, **kwargs)
    for source in insignificant:
        roi.modify(which=source, free=False)

