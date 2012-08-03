from uw.utilities.parmap import LimitMapper

def all_params_limited(roi):
    sources=roi.psm.point_sources.tolist()+roi.dsm.diffuse_sources.tolist()
    for source in sources:
        model = roi.get_model(source)
        for p in model.param_names:
            if not isinstance(model.get_mapper(p),LimitMapper):
                return False
    return True

def zero_prefactor_lower_limits(roi):
    """ Set to 0 the prefactor lower limits for all parameters in the ROI. """
    sources=roi.psm.point_sources.tolist()+roi.dsm.diffuse_sources.tolist()
    for source in sources:
        model=roi.get_model(source)
        mapper = model.get_mapper(0)
        assert isinstance(mapper,LimitMapper)
        new_mapper = mapper.copy()
        new_mapper.lower = 0
        model.set_mapper(0, new_mapper)
