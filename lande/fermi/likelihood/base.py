import yaml
from lande.utilities.tools import tolist
from lande.utilities.save import loaddict

class BaseFitter(object):
    """ Base object for gtlike/pointlike fitter objects.
        
        Object must create a self.results dictionary.
    """

    defaults = (
        ('verbosity', False, 'Make lots of noise'),
    )

    def __init__(self, results):
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
