from collections import Counter

def islist(x):
    return hasattr(x, "__iter__") and not isinstance(x, basestring)

def flatten(x):
    """flatten(sequence) -> list

        Taken from http://kogs-www.informatik.uni-hamburg.de/~meine/python_tricks

    Returns a single, flat list which contains all elements retrieved
    from the sequence and all recursively contained sub-sequences
    (iterables).

    Examples:

        >>> [1, 2, [3,4], (5,6)]
        [1, 2, [3, 4], (5, 6)]
        >>> flatten([[[1,2,3], (42,None)], [4,5], [6], 7, (8,9,10)])
        [1, 2, 3, 42, None, 4, 5, 6, 7, 8, 9, 10]"""

    result = []
    for el in x:
        if islist(el):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result

def duplicates(list):
    """ Find all duplicate items in list.

        >>> duplicates([1,2,3])
        []
        >>> duplicates([1,2,2,2,3])
        [2]
    """
    c = Counter(list)
    return [k for k,v in c.items() if v>1]


if __name__ == "__main__":
    import doctest
    doctest.testmod()
