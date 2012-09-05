#from . binning import close, all_energy, high_energy, higher_energy, one_bin_per_dec, two_bin_per_dec, four_bin_per_dec
import numpy as np

close=lambda x,y: np.allclose(x,y, rtol=0, atol=1)
all_energy=lambda emin,emax: close([emin,emax],[1e2,10**5.5]) or close([emin,emax],[1e2,1e5])
high_energy=lambda emin,emax: close([emin,emax],[10**4,10**5.5])
higher_energy=lambda emin,emax: close([emin,emax],[10**4.5,10**5.5])

def one_bin_per_dec(emin,emax):
    assert close(emin,1e2) and (close(emax,1e5) or close(emax,10**5.5))
    if close(emax,1e5):
        return np.logspace(2,5,4)
    elif close(emax,10**5.5):
        return [1e2,1e3,1e4,10**5.5]

def two_bin_per_dec(emin,emax):
    assert close(emin,1e2) and (close(emax,1e5) or close(emax,10**5.5))
    if close(emax,1e5):
        return np.logspace(2,5,7)
    elif close(emax,10**5.5):
        return np.logspace(2,5.5,8)

def four_bin_per_dec(emin,emax):
    assert close(emin,1e2) and (close(emax,1e5) or close(emax,10**5.5))
    if close(emax,1e5):
        return np.logspace(2,5,13)
    elif close(emax,10**5.5):
        return np.logspace(2,5.5,15)

    
