#!/usr/bin/env python
"""
argsort : argmax, argmin and argsort for array_like and Python iterables.

This module was written by Matthias Cuntz while at Department of
Computational Hydrosystems, Helmholtz Centre for Environmental
Research - UFZ, Leipzig, Germany, and continued while at Institut
National de Recherche pour l'Agriculture, l'Alimentation et
l'Environnement (INRAE), Nancy, France.

Copyright (c) 2014-2020 Matthias Cuntz - mc (at) macu (dot) de
Released under the MIT License; see LICENSE file for details.

* Written Dec 2014 by Matthias Cuntz (mc (at) macu (dot) de)
* Added argmin, argmax, Jul 2019, Matthias Cuntz
* Using numpy docstring format, extending examples from numpy docstrings, May 2020, Matthias Cuntz

.. moduleauthor:: Matthias Cuntz

The following functions are provided

.. autosummary::
   argmax
   argmin
   argsort
"""
from __future__ import division, absolute_import, print_function
import numpy as np


__all__ = ['argmax', 'argmin', 'argsort']


def argmax(a, *args, **kwargs):
    """
    Wrapper for numpy.argmax, numpy.ma.argmax, and using max for Python iterables.

    Passes all keywords directly to the individual routines, i.e.

        numpy.argmax(a, axis=None, out=None)

        numpy.ma.argmax(self, axis=None, fill_value=None, out=None)

    No keyword will be passed to max routine for Python iterables.

    Parameters
    ----------
    a : array_like
        input array, masked array, or Python iterable
    *args : optional
        all arguments of numpy.argmax or numpy.ma.argmax
    **kwargs : optional
        all keyword arguments of numpy.argmax or numpy.ma.argmax

    Returns
    -------
    index_array : ndarray, int
        Array of indices of the largest element in input array `a`.
        It has the same shape as `a.shape` with the dimension along `axis` removed.
        a[np.unravel_index(argmax(a), a.shape)] is the maximum value of a.

    Examples
    --------
    >>> import numpy as np

    # One-dimensional array
    >>> a = np.array([0,4,6,2,1,5,3,5])
    >>> ii = argmax(a)
    >>> print(ii)
    2
    >>> print(a[ii])
    6

    # One-dimensional masked array
    >>> a = np.ma.array([0,4,6,2,1,5,3,5], mask=[0,0,1,1,0,0,0,0])
    >>> ii = argmax(a)
    >>> print(ii)
    5
    >>> print(a[ii])
    5
    >>> ii = argmax(a, fill_value=6)
    >>> print(ii)
    2

    # List
    >>> a = [0,4,6,2,1,5,3,5]
    >>> ii = argmax(a)
    >>> print(ii)
    2
    >>> print(a[ii])
    6

    >>> # from numpy.argmax docstring
    >>> a = np.arange(6).reshape(2,3) + 10
    >>> a
    array([[10, 11, 12],
           [13, 14, 15]])
    >>> np.argmax(a)
    5
    >>> np.argmax(a, axis=0)
    array([1, 1, 1])
    >>> np.argmax(a, axis=1)
    array([2, 2])

    # Indexes of the maximal elements of a N-dimensional array:
    >>> ind = np.unravel_index(np.argmax(a, axis=None), a.shape)
    >>> ind
    (1, 2)
    >>> a[ind]
    15

    >>> b = np.arange(6)
    >>> b[1] = 5
    >>> b
    array([0, 5, 2, 3, 4, 5])
    >>> np.argmax(b)  # Only the first occurrence is returned.
    1

    Notes
    -----
    argmax for iterables was taken from
    https://stackoverflow.com/questions/16945518/finding-the-index-of-the-value-which-is-the-min-or-max-in-python

    History
    -------
    Written,  Matthias Cuntz, Jul 2019
    Modified, Matthias Cuntz, May 2020 - numpy docstring, more examples from numpy docs
    """
    if isinstance(a, np.ma.MaskedArray):
        return np.ma.argmax(a, *args, **kwargs)
    elif isinstance(a, np.ndarray):
        return np.argmax(a, *args, **kwargs)
    else:
        return _argmax(a)


def argmin(a, *args, **kwargs):
    """
    Wrapper for numpy.argmin, numpy.ma.argmin, and using min for Python iterables.

    Passes all keywords directly to the individual routines, i.e.
        numpy.argmin(a, axis=None, out=None)

        numpy.ma.argmin(self, axis=None, fill_value=None, out=None)

        No keyword will be passed to min routine for Python iterables.

    Parameters
    ----------
    a : array_like
        input array, masked array, or Python iterable
    *args : optional
        all arguments of numpy.argmin or numpy.ma.argmin
    **kwargs : optional
        all keyword arguments of numpy.argmin or numpy.ma.argmin

    Returns
    -------
    index_array : ndarray, int
        Array of indices of the largest element in input array `a`.
        It has the same shape as `a.shape` with the dimension along `axis` removed.
        a[np.unravel_index(argmin(a), a.shape)] is the minimum value of a.

    Examples
    --------
    >>> import numpy as np

    # One-dimensional array
    >>> a = np.array([0,4,6,2,1,5,3,5])
    >>> ii = argmin(a)
    >>> print(ii)
    0
    >>> print(a[ii])
    0

    # One-dimensional masked array
    >>> a = np.ma.array([0,4,6,2,1,5,3,5], mask=[1,0,1,1,0,0,0,0])
    >>> ii = argmin(a)
    >>> print(ii)
    4
    >>> print(a[ii])
    1
    >>> ii = argmin(a, fill_value=1)
    >>> print(ii)
    0

    # List
    >>> a = [0,4,6,2,1,5,3,5]
    >>> ii = argmin(a)
    >>> print(ii)
    0
    >>> print(a[ii])
    0

    >>> # from numpy.argmin docstring
    >>> a = np.arange(6).reshape(2,3) + 10
    >>> a
    array([[10, 11, 12],
           [13, 14, 15]])
    >>> np.argmin(a)
    0
    >>> np.argmin(a, axis=0)
    array([0, 0, 0])
    >>> np.argmin(a, axis=1)
    array([0, 0])

    # Indices of the minimum elements of a N-dimensional array:
    >>> ind = np.unravel_index(np.argmin(a, axis=None), a.shape)
    >>> ind
    (0, 0)
    >>> a[ind]
    10

    >>> b = np.arange(6) + 10
    >>> b[4] = 10
    >>> b
    array([10, 11, 12, 13, 10, 15])
    >>> np.argmin(b)  # Only the first occurrence is returned.
    0

    Notes
    -----
    argmin for iterables was taken from
    https://stackoverflow.com/questions/16945518/finding-the-index-of-the-value-which-is-the-min-or-max-in-python

    History
    -------
    Written,  Matthias Cuntz, Jul 2019
    Modified, Matthias Cuntz, May 2020 - numpy docstring, more examples from numpy docs
    """
    if isinstance(a, np.ma.MaskedArray):
        return np.ma.argmin(a, *args, **kwargs)
    elif isinstance(a, np.ndarray):
        return np.argmin(a, *args, **kwargs)
    else:
        return _argmin(a)


def argsort(a, *args, **kwargs):
    """
    Wrapper for numpy.argsort, numpy.ma.argsort, and using sorted for Python iterables.

    Passes all keywords directly to the individual routines, i.e.
        numpy.argsort(a, axis=-1, kind='quicksort', order=None)
        numpy.ma.argsort(a, axis=None, kind='quicksort', order=None, fill_value=None)
        sorted(iterable[, cmp[, key[, reverse]]])

    Only key cannot be given for Python iterables because the input array is used
    as key in the sorted function.

    Parameters
    ----------
    a : array_like
        input array, masked array, or Python iterable
    *args : optional
        all arguments of numpy.argsort, numpy.ma.argsort, and sorted (except key argument)
    **kwargs : optional
        all keyword arguments of numpy.argsort, numpy.ma.argsort, and sorted (except key argument)

    Returns
    -------
    index_array : ndarray, int
        Array of indices that sort `a` along the specified `axis`.
        If `a` is one-dimensional, ``a[index_array]`` yields a sorted `a`.

    Examples
    --------
    >>> import numpy as np

    # 1D array
    >>> a = np.array([0,4,6,2,1,5,3,5])
    >>> ii = argsort(a)
    >>> print(a[ii])
    [0 1 2 3 4 5 5 6]

    >>> ii = argsort(a, kind='quicksort')
    >>> print(a[ii])
    [0 1 2 3 4 5 5 6]

    # 1D masked array
    >>> a = np.ma.array([0,4,6,2,1,5,3,5], mask=[0,0,1,1,0,0,0,0])
    >>> ii = argsort(a)
    >>> print(a[ii])
    [0 1 3 4 5 5 -- --]

    >>> ii = argsort(a, fill_value=1)
    >>> print(a[ii])
    [0 -- -- 1 3 4 5 5]

    # list
    >>> a = [0,4,6,2,1,5,3,5]
    >>> ii = argsort(a)
    >>> b = [ a[i] for i in ii ]
    >>> print(b)
    [0, 1, 2, 3, 4, 5, 5, 6]

    >>> a = [0,4,6,2,1,5,3,5]
    >>> ii = argsort(a, reverse=True)
    >>> b = [ a[i] for i in ii ]
    >>> print(b)
    [6, 5, 5, 4, 3, 2, 1, 0]


    # from numpy.argsort docstring
    # One-dimensional array:
    >>> x = np.array([3, 1, 2])
    >>> np.argsort(x)
    array([1, 2, 0])

    >>> # Two-dimensional array:
    >>> x = np.array([[0, 3], [2, 2]])
    >>> x
    array([[0, 3],
           [2, 2]])
    >>> ind = np.argsort(x, axis=0)  # sorts along first axis (down)
    >>> ind
    array([[0, 1],
           [1, 0]])
    >>> np.take_along_axis(x, ind, axis=0)  # same as np.sort(x, axis=0)
    array([[0, 2],
           [2, 3]])
    >>> ind = np.argsort(x, axis=1)  # sorts along last axis (across)
    >>> ind
    array([[0, 1],
           [0, 1]])
    >>> np.take_along_axis(x, ind, axis=1)  # same as np.sort(x, axis=1)
    array([[0, 3],
           [2, 2]])

    # Indices of the sorted elements of a N-dimensional array:
    >>> ind = np.unravel_index(np.argsort(x, axis=None), x.shape)
    >>> ind
    (array([0, 1, 1, 0]), array([0, 0, 1, 1]))
    >>> x[ind]  # same as np.sort(x, axis=None)
    array([0, 2, 2, 3])

    >>> # Sorting with keys:
    >>> x = np.array([(1, 0), (0, 1)], dtype=[('x', '<i4'), ('y', '<i4')])
    >>> x
    array([(1, 0), (0, 1)],
          dtype=[('x', '<i4'), ('y', '<i4')])
    >>> np.argsort(x, order=('x','y'))
    array([1, 0])
    >>> np.argsort(x, order=('y','x'))
    array([0, 1])

    Notes
    -----
    argsort for iterables was taken from
    http://stackoverflow.com/questions/3382352/equivalent-of-numpy-argsort-in-basic-python
    http://stackoverflow.com/questions/3071415/efficient-method-to-calculate-the-rank-vector-of-a-list-in-python

    History
    -------
    Written,  Matthias Cuntz, Dec 2014
    Modified, Matthias Cuntz, May 2020 - numpy docstring, more examples from numpy docs
    """
    if isinstance(a, np.ma.MaskedArray):
        return np.ma.argsort(a, *args, **kwargs)
    elif isinstance(a, np.ndarray):
        return np.argsort(a, *args, **kwargs)
    else:
        return _argsort(a, *args, **kwargs)


# same as numpy.argmax but for python iterables
def _argmax(iterable):
    return max(enumerate(iterable), key=lambda x: x[1])[0]

# same as numpy.argmin but for python iterables
def _argmin(iterable):
    return min(enumerate(iterable), key=lambda x: x[1])[0]


# same as numpy.argsort but for python iterables
def _argsort(seq, *args, **kwargs):
    if 'key' in kwargs:
        raise KeyError('keyword key cannot be given to argsort.')
    return sorted(range(len(seq)), *args, key=seq.__getitem__, **kwargs)


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
