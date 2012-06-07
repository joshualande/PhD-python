""" Classes to save the state of fit parameters.
    
    Works like LikelihoodState with the added advantage that
    
    (a) it will change back spectral models that have changed
        This is useful if you are trying out a new spectral model or something

    (b) Allows the parameter like to be passed into restore to
        restore into a new like object.

    (c) fixes a bounds error bug (see LK-73).

"""
import pyLikelihood


class _Parameter(object):
    """ Copy of pyLikelihood.LikelihoodState with my attempted fix
        to a bounds error bug (see LK-73). """
    def __init__(self, par):
        self.par = par
        self.value = par.getValue()
        self.minValue, self.maxValue = par.getBounds()
        self.free = par.isFree()
        self.scale = par.getScale()
        self.error = par.error()
        self.alwaysFixed = par.alwaysFixed()
    def setDataMembers(self, par=None):
        if par is None:
            par = self.par
        par.setBounds(-float('inf'), float('inf'))
        par.setValue(self.value)
        par.setBounds(self.minValue, self.maxValue)
        par.setFree(self.free)
        par.setScale(self.scale)
        par.setError(self.error)
        par.setAlwaysFixed(self.alwaysFixed)


class SuperState(object):
    def __init__(self, like):
        self.like = like
        self.sources = dict()

        all_names = like.sourceNames()

        for name in all_names:

            spectrum = like[name].src.spectrum()

            parameters=pyLikelihood.ParameterVector()
            spectrum.getParams(parameters)

            self.sources[name] = d = dict(parameters=dict(), spectrum=spectrum)
            for p in parameters:
                d['parameters'][p.getName()] = _Parameter(p)

    def restore(self, like=None):
        if like is None: like = self.like

        for sname,v in self.sources.items():

            spectrum = v['spectrum']
            parameters = v['parameters']

            # Set back to old spectrum
            like.setSpectrum(sname,spectrum)

            # Reset all parameters
            for pname,pcache in parameters.items():
                index = like.par_index(sname, pname)
                like_par = like.params()[index]
                pcache.setDataMembers(like_par)

        like.syncSrcParams()
