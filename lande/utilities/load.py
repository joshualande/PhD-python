import imp

def import_module(filename):
    """ import a python module from a pathname. 

        Usage:
            module = import_module('/path/to/module.py')
            ...
    """
    return imp.load_source('module',filename)

