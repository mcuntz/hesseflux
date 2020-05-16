#!/usr/bin/env python
"""
division : Divide two arrays, return `otherwise` if division by 0.

This module was written by Matthias Cuntz while at Department of
Computational Hydrosystems, Helmholtz Centre for Environmental
Research - UFZ, Leipzig, Germany, and continued while at Institut
National de Recherche pour l'Agriculture, l'Alimentation et
l'Environnement (INRAE), Nancy, France.

Copyright (c) 2012-2020 Matthias Cuntz - mc (at) macu (dot) de
Released under the MIT License; see LICENSE file for details.

* Written Jan 2012 by Matthias Cuntz (mc (at) macu (dot) de)
* Added wrapper div, May 2012, Matthias Cuntz
* Ported to Python 3, Feb 2013, Matthias Cuntz
* Do not return masked array if no masked array given, Oct 2014, Matthias Cuntz
* Added two-digit year, Nov 2018, Matthias Cuntz
* Removed bug that non-masked array was returned if masked array given, Sep 2015, Matthias Cuntz
* Using numpy docstring format, May 2020, Matthias Cuntz

.. moduleauthor:: Matthias Cuntz

The following functions are provided

.. autosummary::
   division
   div
"""
from __future__ import division, absolute_import, print_function
import numpy as np


__all__ = ['division', 'div']


def division(a, b, otherwise=np.nan, prec=0.):
    """
    Divide two arrays, return `otherwise` if division by 0.

    Parameters
    ----------
    a : array_like
        enumerator
    b : array_like
        denominator
    otherwise : float
        value to return if `b=0` (default: `np.nan`)
    prec : float
        if `|b|<|prec|` then `otherwise`

    Returns
    -------
    ratio : numpy array or masked array
        a/b        if `|b| >  |prec|`

        otherwise  if `|b| <= |prec|`

        Output is numpy array. It is a masked array if at least one
        of `a` or `b` is a masked array.

    Examples
    --------
    >>> a = [1., 2., 3.]
    >>> b = 2.
    >>> print('{:.1f} {:.1f} {:.1f}'.format(*division(a, b)))
    0.5 1.0 1.5

    >>> a = [1., 1., 1.]
    >>> b = [2., 1., 0.]
    >>> print('{:.1f} {:.1f} {:.1f}'.format(*division(a, b)))
    0.5 1.0 nan
    >>> print('{:.1f} {:.1f} {:.1f}'.format(*division(a, b, 0.)))
    0.5 1.0 0.0
    >>> print('{:.1f} {:.1f} {:.1f}'.format(*division(a, b, otherwise=0.)))
    0.5 1.0 0.0
    >>> print('{:.1f} {:.1f} {:.1f}'.format(*division(a, b, prec=1.)))
    0.5 nan nan

    >>> import numpy as np
    >>> a = np.array([1., 1., 1.])
    >>> b = [2., 1., 0.]
    >>> print('{:.1f} {:.1f} {:.1f}'.format(*division(a, b)))
    0.5 1.0 nan

    >>> b = np.array([2., 1., 0.])
    >>> print('{:.1f} {:.1f} {:.1f}'.format(*division(a, b)))
    0.5 1.0 nan

    >>> mask = [0, 0, 1]
    >>> b = np.ma.array([2., 1., 0.], mask=mask)
    >>> ratio = division(a, b)
    >>> print(ratio)
    [0.5 1.0 --]


    History
    -------
    Written,  Matthias Cuntz, Jan 2012
    Modified, Matthias Cuntz, May 2012 - div
              Matthias Cuntz, Feb 2013 - ported to Python 3
              Matthias Cuntz, Oct 2014 - do not return masked array if no masked array given
              Matthias Cuntz, Sep 2015 - bug: returned non-masked array in case of masked array input
              Matthias Cuntz, May 2020 - numpy docstring format
    """
    oldsettings = np.geterr()
    np.seterr(divide='ignore')

    if isinstance(a, np.ma.masked_array) or isinstance(b, np.ma.masked_array):
        out = np.ma.where(np.ma.abs(np.ma.array(b)) > np.abs(prec), np.ma.array(a)/np.ma.array(b), otherwise)
    else:
        out = np.where(np.abs(np.array(b)) > np.abs(prec), np.array(a)/np.array(b), otherwise)

    np.seterr(**oldsettings)

    return out


def div(*args, **kwargs):
    """
    Wrapper function for :func:`division`.

    Examples
    --------
    >>> a = [1., 2., 3.]
    >>> b = 2.
    >>> print('{:.1f} {:.1f} {:.1f}'.format(*div(a, b)))
    0.5 1.0 1.5

    >>> a = [1., 1., 1.]
    >>> b = [2., 1., 0.]
    >>> print('{:.1f} {:.1f} {:.1f}'.format(*div(a, b)))
    0.5 1.0 nan
    >>> print('{:.1f} {:.1f} {:.1f}'.format(*div(a, b, 0.)))
    0.5 1.0 0.0
    >>> print('{:.1f} {:.1f} {:.1f}'.format(*div(a, b, otherwise=0.)))
    0.5 1.0 0.0
    >>> print('{:.1f} {:.1f} {:.1f}'.format(*div(a, b, prec=1.)))
    0.5 nan nan
    """
    return division(*args, **kwargs)


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

