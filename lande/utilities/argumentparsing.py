import sys
import argparse

def parse_strip_known_args(parser):
    """ parser is an argparse object.
        Parse the args that have been spcified and return
        the Namespace object attributed from it. But
        also remove those arguments from the command line
        arguments so a future argument parser can
        work on the remaining arguments. Useful if 
        you want you function to just grab one or two
        arguments hand pass the others to some more 
        general parsing code. """
    args,extra=parser.parse_known_args()
    sys.argv=[sys.argv[0]]+extra # remove from argv the parsed flags
    return args

def argparse_to_kwargs(args):
    """ Takes in an argparse 'args' object and returns a dictionary of the
        parameters.  """
    assert isinstance(args,argparse.Namespace)
    return dict(args._get_kwargs())
