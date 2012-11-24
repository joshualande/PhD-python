import numpy as np

def overlaps(x1,x2, y1,y2):
    """Given two numeric ranges, returns a flag True or False
       indicating whether they overlap.

       Assumes that the ranges are ordered (smallest,largest). If
       that assumption is wrong, incorrect results may occur.

        Code taken from:
            http://bytes.com/topic/python/answers/457949-determing-whether-two-ranges-overlap

        Partially overlapping:

            >>> overlaps(0,5,4,10)
            True
            >>> overlaps(4,10,0,5)
            True

        Not overlapping:

            >>> overlaps(0,4,5,10)
            False
            >>> overlaps(5,10,0,4)
            False
            >>> overlaps(10,15,0,5)
            False

        Fully Overlapping:
            >>> overlaps(0,10,4,5)
            True
            >>> overlaps(4,5,0,10)
            True
    """

    # Fully overlapping cases:
    # x1 <= y1 <= y2 <= x2
    # y1 <= x1 <= x2 <= y2
    # Partially overlapping cases:
    # x1 <= y1 <= x2 <= y2
    # y1 <= x1 <= y2 <= x2
    # Non-overlapping cases:
    # x1 <= x2 < y1 <= y2
    # y1 <= y2 < x1 <= x2
    assert x1<x2 and y1<y2
    return not (x2 < y1 or y2 < x1)


def log_binomial(n,m):
    """ Log of the binomial coefficient (n choose m).

        Implementation taken from 
            http://grocid.net/2012/07/02/quick-and-dirty-way-to-calculate-large-binomial-coefficients-in-python/
    """
    log_fac = lambda n: np.sum(np.log(np.arange(2,n+1)))
    return log_fac(n)-log_fac(m)-log_fac(n-m)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
