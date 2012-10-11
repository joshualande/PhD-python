
def gtlike_modify(like, name, free=None, freeze_spectral_shape=None):
    """ Freeze a source in a gtlike ROI. 
    
        The method for modifying the ROI
        follows the code in SuperState.py 
        
        I am not sure why the modificaiton
        has to be done in this particular way. """
    from pyLikelihood import StringVector

    source = like.logLike.getSource(name)
    spectrum = like[name].src.spectrum()

    if free is not None:
        parNames = StringVector()
        spectrum.getParamNames(parNames)
        for parName in parNames:
            index = like.par_index(name, parName)
            par = like.params()[index]
            par.setFree(free)

    if freeze_spectral_shape:

        normpar = like.normPar(name)
        norm_free=normpar.isFree()

        gtlike_modify(like,name,free=False)

        normpar = like.normPar(name)
        normpar.setFree(norm_free)

    like.syncSrcParams(name)

modify=gtlike_modify
