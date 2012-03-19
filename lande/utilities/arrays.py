import numpy as np

def isclose(a,b, rtol=1e-05, atol=1e-08):
    return np.abs(a - b) <= (atol + rtol * np.abs(b))

def nzip(a,b):
    """ Zip two arrays together in a vector fashion:
        
            >>> np.all(nzip([1,2,3],[4,5,6]) == [1,4,2,5,3,6])
            True
    """
    assert len(a) == len(b)
    c=np.empty(2*len(a))
    c[::2] = a
    c[1::2] = b
    return c

if __name__ == "__main__":
    import doctest
    doctest.testmod()
