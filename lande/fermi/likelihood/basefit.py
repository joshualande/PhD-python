import yaml

from uw.utilities import keyword_options

from lande.utilities.tools import tolist
from lande.utilities.save import loaddict

class BaseFitter(object):
    """ BaseFitter is a base class for all of my
        gtlike/pointlike fitter objects.

        The intention of this code is to provide a uniform 
        interface to performing a given analysis with gtlike/pointlike:

            analysis = AnalysisObject(like, name, verbosity=True, other_params)

            # create a dictionary of the results
            d = analysis.todict()

            # save results to a YAML file
            analysis.save('results.yaml')
        
        Object must create a self.results dictionary.
    """

    defaults = (
        ('verbosity', False, 'Make lots of noise'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, results, **kwargs):
        keyword_options.process(self, kwargs)

        if isinstance(results,dict):
            self.results = results
        elif isinstance(results, str):
            self.results = loaddict(results)
        else:
            raise Exception("Unrecognized results %s" % results)

    def todict(self):
        """ Pacakge up the results of the SED fit into
            a nice dictionary. """
        return tolist(self.results)

    def __str__(self):
        results = self.todict()
        return yaml.dump(results)

    def save(self,filename,**kwargs):
        """ Save SED data points to a file. """
        if hasattr(filename,'write'):
            filename.write(self.__str__())
        else:
            f=open(filename,'w')
            f.write(self.__str__())
            f.close()
