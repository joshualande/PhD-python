from math import *


from scipy.stats import chi2
import scipy.special as sp


def TS2sigma(TS,dof):
    # one-sided Chi^2 test
    pval_1 = chi2.cdf(TS, dof)
    sigma=sqrt(2)*sp.erfinv(pval_1)
    print "TS=%.2f\t->\t%.2f sigma"%(TS,sigma)
    return sigma
