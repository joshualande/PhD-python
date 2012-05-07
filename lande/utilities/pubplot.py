import sys
from argparse import ArgumentParser

import pylab as P
from matplotlib import rc

def set_latex_defaults():
    rc('ps',usedistiller='xpdf')
    rc('text', usetex=True)
    rc('font', family='serif', serif="Computer Modern Roman")



def get_bw():
    parser = ArgumentParser()
    parser.add_argument("--bw", action="store_true", default=False)
    args,extra=parser.parse_known_args()
    sys.argv=[sys.argv[0]]+extra # remove from argv the --bw flag

    global bw
    bw=args.bw
    return bw


def save(base):
    if bw:
        P.savefig('%s_bw.pdf' % base)
        P.savefig('%s_bw.eps' % base)
    else:
        P.savefig('%s_color.pdf' % base)
        P.savefig('%s_color.eps' % base)

