import numpy as np

def isclose(a,b, rtol=1e-05, atol=1e-08):
    return np.abs(a - b) <= (atol + rtol * np.abs(b))

