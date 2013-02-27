import operator

import numpy as np

def recursive_allclose(a,b, *args, **wargs):
    """ 
        >>> np.allclose(1e-20, 1.1e-20)
        True
        >>> np.allclose(1e-20, 1e-5)
        False
        >>> recursive_allclose(dict(a=dict(b=1e-20)), dict(a=dict(b=1.1e-20)))
        True
        >>> recursive_allclose(dict(a=dict(b=1e-20)), dict(a=dict(b=1e-5)))
        False
        >>> recursive_allclose(dict(a=[1e-20]),dict(a=[1.1e-20]))
        True
        >>> recursive_allclose(dict(a=[1e-20]),dict(a=[1e-5]))
        False
        >>> recursive_allclose('cat','cat')
        True
        >>> recursive_allclose('cat','dog')
        False
    """
    assert False,"This function was never finished"
    if isinstance(a,dict):
        assert isinstance(b,dict) and a.keys() == b.keys(), '%s, %s' % (a.keys(),b.keys())
        return reduce(operator.and_,[recursive_allclose(a[i],b[i], *args, **wargs) for i in a.keys()])
    elif isinstance(a,str):
        assert isinstance(b,str), '%s, %s' % (type(a),type(b))
        return a==b
    else:
        return np.allclose(a,b, *args, **wargs)

def almost_equal(a,b, rtol=1e-05, atol=1e-08):
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
