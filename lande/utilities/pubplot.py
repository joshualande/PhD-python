import sys
from argparse import ArgumentParser
from os.path import expandvars

import pylab as P
from matplotlib import rc

from . argumentparsing import parse_strip_known_args

def set_latex_defaults():
    rc('ps',usedistiller='xpdf')
    rc('text', usetex=True)
    rc('font', family='serif', serif="Computer Modern Roman")



def get_bw():
    parser = ArgumentParser()
    parser.add_argument("--bw", action="store_true", default=False)
    args=parse_strip_known_args(parser)
    global bw
    bw=args.bw
    return bw


def save(base):
    base=expandvars(base)
    if bw:
        P.savefig('%s_bw.pdf' % base)
        P.savefig('%s_bw.eps' % base)
    else:
        P.savefig('%s_color.pdf' % base)
        P.savefig('%s_color.eps' % base)

