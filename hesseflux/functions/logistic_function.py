#!/usr/bin/env python
"""
Module with several forms of the logistic function and its first and second derivatives.

The current functions are:
    logistic              Logistic function L/(1+exp(-k(x-x0)))

    logistic_p            logistic(x,*p)

    dlogistic             First derivative of logistic function

    dlogistic_p           dlogistic(x,*p)

    d2logistic            Second derivative of logistic function

    d2logistic_p          d2logistic(x,*p)

    logistic_offset       logistic function with offset L/(1+exp(-k(x-x0))) + a

    logistic_offset_p     logistic_offset(x,*p)

    dlogistic_offset      First derivative of logistic function with offset

    dlogistic_offset_p    dlogistic_offset(x,*p)

    d2logistic_offset     Second derivative of logistic function with offset

    d2logistic_offset_p   d2logistic_offset(x,*p)

    logistic2_offset      Double logistic function with offset L1/(1+exp(-k1(x-x01))) - L2/(1+exp(-k2(x-x02))) + a

    logistic2_offset_p    logistic2_offset(x,*p)

    dlogistic2_offset     First derivative of double logistic function with offset

    dlogistic2_offset_p   dlogistic2_offset(x,*p)

    d2logistic2_offset    Second derivative of double logistic function with offset

    d2logistic2_offset_p  d2logistic2_offset(x,*p)

This module was written by Matthias Cuntz while at Department of
Computational Hydrosystems, Helmholtz Centre for Environmental
Research - UFZ, Leipzig, Germany, and continued while at Institut
National de Recherche pour l'Agriculture, l'Alimentation et
l'Environnement (INRAE), Nancy, France.

Copyright (c) 2015-2020 Matthias Cuntz - mc (at) macu (dot) de
Released under the MIT License; see LICENSE file for details.

* Written Mar 2015 by Matthias Cuntz (mc (at) macu (dot) de)
* Added functions logistic_p and logistic_offset_p, Dec 2017, Matthias Cuntz
* Changed to Sphinx docstring and numpydoc, Dec 2019, Matthias Cuntz
* Distinguish iterable and array_like parameter types, Jan 2020, Matthias Cuntz
* Make systematically function_p versions of all logistic functions and its derivatives, Feb 2020, Matthias Cuntz
* Split logistic and curvature into separate files, May 2020, Matthias Cuntz

.. moduleauthor:: Matthias Cuntz

The following functions are provided:

.. autosummary::
    logistic
    logistic_p
    dlogistic
    dlogistic_p
    d2logistic
    d2logistic_p
    logistic_offset
    logistic_offset_p
    dlogistic_offset
    dlogistic_offset_p
    d2logistic_offset
    d2logistic_offset_p
    logistic2_offset
    logistic2_offset_p
    dlogistic2_offset
    dlogistic2_offset_p
    d2logistic2_offset
    d2logistic2_offset_p
"""
from __future__ import division, absolute_import, print_function
import numpy as np
import scipy.special as sp


__all__ = ['logistic', 'logistic_p',
           'dlogistic', 'dlogistic_p',
           'd2logistic', 'd2logistic_p',
           'logistic_offset', 'logistic_offset_p',
           'dlogistic_offset', 'dlogistic_offset_p',
           'd2logistic_offset', 'd2logistic_offset_p',
           'logistic2_offset', 'logistic2_offset_p',
           'dlogistic2_offset', 'dlogistic2_offset_p',
           'd2logistic2_offset', 'd2logistic2_offset_p']


# -----------------------------------------------------------
# a/(1+exp(-b(x-c))) - logistic function
def logistic(x, L, k, x0):
    """
    Logistic function:

        `L/(1+exp(-k(x-x0)))`

    Parameters
    ----------
    x : array_like
        Independent variable to evalute logistic function
    L : float
        Maximum of logistic function
    k : float
        Steepness of logistic function
    x0 : float
        Inflection point of logistic function

    Returns
    -------
    float or ndarray
        Logistic function at `x` with maximum `L`, steepness `k` and inflection point `x0`
    """
    return L * sp.expit(k * (x - x0))


def logistic_p(x, p):
    """
    Wrapper function for :func:`logistic`: `logistic(x, *p)`.
    """
    return logistic(x, *p)


# -----------------------------------------------------------
# 1st derivative of logistic functions
def dlogistic(x, L, k, x0):
    """
    First derivative of logistic function:

        `L/(1+exp(-k(x-x0)))`

    which is

        `k.L/(2(cosh(k(x-x0))+1))`

    Parameters
    ----------
    x : array_like
        Independent variable to evalute derivative of logistic function
    L : float
        Maximum of logistic function
    k : float
        Steepness of logistic function
    x0 : float
        Inflection point of logistic function

    Returns
    -------
    float or ndarray
        First derivative of logistic function at `x` with maximum `L`, steepness `k` and inflection point `x0`
    """
    return k * L / (2. * (np.cosh(k * (x - x0)) + 1.))


def dlogistic_p(x, p):
    """
    Wrapper function for :func:`dlogistic`: `dlogistic(x, *p)`.
    """
    return dlogistic(x, *p)


# -----------------------------------------------------------
# 2nd derivative of logistic functions
def d2logistic(x, L, k, x0):
    """
    Second derivative of logistic function:

        `L/(1+exp(-k(x-x0)))`

    which is

        `-k^2.L.sinh(k(x-x0))/(2(cosh(k(x-x0))+1)^2)`

    Parameters
    ----------
    x : array_like
        Independent variable to evalute derivative of logistic function
    L : float
        Maximum of logistic function
    k : float
        Steepness of logistic function
    x0 : float
        Inflection point of logistic function

    Returns
    -------
    float or ndarray
        Second derivative of logistic function at `x` with maximum `L`, steepness `k` and inflection point `x0`
    """
    return ( -k**2 * L * np.sinh(k * (x - x0)) /
                 (2. * (np.cosh(k * (x - x0)) + 1.)**2) )


def d2logistic_p(x, p):
    """
    Wrapper function for :func:`d2logistic`: `d2logistic(x, *p)`.
    """
    return d2logistic(x, *p)


# -----------------------------------------------------------
# L/(1+exp(-k(x-x0))) + a - logistic function with offset
def logistic_offset(x, L, k, x0, a):
    """
    Logistic function with offset:

        `L/(1+exp(-k(x-x0))) + a`

    Parameters
    ----------
    x : array_like
        Independent variable to evalute logistic function
    L : float
        Maximum of logistic function
    k : float
        Steepness of logistic function
    x0 : float
        Inflection point of logistic function
    a : float
        Offset of logistic function

    Returns
    -------
    float or ndarray
        Logistic function at `x` with maximum `L`, steepness `k`, inflection point `x0` and offset `a`
    """
    return L * sp.expit(k * (x - x0)) + a


def logistic_offset_p(x, p):
    """
    Wrapper function for :func:`logistic_offset`: `logistic_offset(x, *p)`.
    """
    return logistic_offset(x, *p)


# -----------------------------------------------------------
# 1st derivative of logistic functions with offset
def dlogistic_offset(x, L, k, x0, a):
    """
    First derivative of logistic function with offset:

        `L/(1+exp(-k(x-x0))) + a`

    which is

        `k.L/(2(cosh(k(x-x0))+1))`

    Parameters
    ----------
    x : array_like
        Independent variable to evalute derivative of logistic function
    L : float
        Maximum of logistic function
    k : float
        Steepness of logistic function
    x0 : float
        Inflection point of logistic function
    a : float
        Offset of logistic function

    Returns
    -------
    float or ndarray
        First derivative of logistic function with offset at `x` with maximum `L`, steepness `k`,
        inflection point `x0`, and offset `a`
    """
    return k * L / (2. * (np.cosh(k * (x - x0)) + 1.))


def dlogistic_offset_p(x, p):
    """
    Wrapper function for :func:`dlogistic_offset`: `dlogistic_offset(x, *p)`.
    """
    return dlogistic_offset(x, *p)


# -----------------------------------------------------------
# 2nd derivative of logistic functions with offset
def d2logistic_offset(x, L, k, x0, a):
    """
    Second derivative of logistic function with offset

        `L/(1+exp(-k(x-x0))) + a`

    which is

        `-k^2.L.sinh(k(x-x0))/(2(cosh(k(x-x0))+1)^2)`

    Parameters
    ----------
    x : array_like
        Independent variable to evalute derivative of logistic function
    L : float
        Maximum of logistic function
    k : float
        Steepness of logistic function
    x0 : float
        Inflection point of logistic function
    a : float
        Offset of logistic function

    Returns
    -------
    float or ndarray
        Second derivative of logistic function at `x` with maximum `L`, steepness `k`,
        inflection point `x0`, and offset `a`
    """
    return ( -k**2 * L * np.sinh(k * (x - x0)) /
                 (2. * (np.cosh(k * (x - x0)) + 1.)**2) )


def d2logistic_offset_p(x, p):
    """
    Wrapper function for :func:`d2logistic_offset`: `d2logistic_offset(x, *p)`.
    """
    return d2logistic_offset(x, *p)


# -----------------------------------------------------------
# L/(1+exp(-k(x-x0))) + a - logistic function with offset
def logistic2_offset(x, L1, k1, x01, L2, k2, x02, a):
    """
    Double logistic function with offset:

        `L1/(1+exp(-k1(x-x01))) - L2/(1+exp(-k2(x-x02))) + a`

    Parameters
    ----------
    x : array_like
        Independent variable to evalute logistic function
    L1 : float
        Maximum of first logistic function
    k1 : float
        Steepness of first logistic function
    x01 : float
        Inflection point of first logistic function
    L2 : float
        Maximum of second logistic function
    k2 : float
        Steepness of second logistic function
    x02 : float
        Inflection point of second logistic function
    a : float
        Offset of double logistic function

    Returns
    -------
    float or ndarray
        Double Logistic function at `x`
    """
    return L1 * sp.expit(k1 * (x - x01)) - L2 * sp.expit(k2 * (x - x02)) + a


def logistic2_offset_p(x, p):
    """
    Wrapper function for :func:`logistic2_offset`: `logistic2_offset(x, *p)`.
    """
    return logistic2_offset(x, *p)


# -----------------------------------------------------------
# 1st derivative of logistic functions with offset
def dlogistic2_offset(x, L1, k1, x01, L2, k2, x02, a):
    """
    First derivative of double logistic function with offset:

        `L1/(1+exp(-k1(x-x01))) - L2/(1+exp(-k2(x-x02))) + a`

    which is

        `k1.L1/(2(cosh(k1(x-x01))+1)) - k2.L2/(2(cosh(k2(x-x02))+1))`

    Parameters
    ----------
    x : array_like
        Independent variable to evalute logistic function
    L1 : float
        Maximum of first logistic function
    k1 : float
        Steepness of first logistic function
    x01 : float
        Inflection point of first logistic function
    L2 : float
        Maximum of second logistic function
    k2 : float
        Steepness of second logistic function
    x02 : float
        Inflection point of second logistic function
    a : float
        Offset of double logistic function

    Returns
    -------
    float or ndarray
        First derivative of double logistic function with offset at `x`
    """
    return ( k1 * L1 / (2. * (np.cosh(k1 * (x - x01)) + 1.)) -
             k2 * L2 / (2. * (np.cosh(k2 * (x - x02)) + 1.)) )


def dlogistic2_offset_p(x, p):
    """
    Wrapper function for :func:`dlogistic2_offset`: `dlogistic2_offset(x, *p)`.
    """
    return dlogistic2_offset(x, *p)


# -----------------------------------------------------------
# 2nd derivative of logistic functions with offset
def d2logistic2_offset(x, L1, k1, x01, L2, k2, x02, a):
    """
    Second derivative of double logistic function with offset:

        `L1/(1+exp(-k1(x-x01))) - L2/(1+exp(-k2(x-x02))) + a`

    which is

        `-k1^2.L1.sinh(k1(x-x01))/(2(cosh(k1(x-x01))+1)^2) +k2^2.L2.sinh(k2(x-x02))/(2(cosh(k2(x-x02))+1)^2)`

    Parameters
    ----------
    x : array_like
        Independent variable to evalute logistic function
    L1 : float
        Maximum of first logistic function
    k1 : float
        Steepness of first logistic function
    x01 : float
        Inflection point of first logistic function
    L2 : float
        Maximum of second logistic function
    k2 : float
        Steepness of second logistic function
    x02 : float
        Inflection point of second logistic function
    a : float
        Offset of double logistic function

    Returns
    -------
    float or ndarray
        Second derivative of double logistic function with offset at `x`
    """
    return ( -k1**2 * L1 * np.sinh(k1 * (x - x01)) /
                 (2. * (np.cosh(k1 * (x - x01)) + 1.)**2) +
             k2**2 * L2 * np.sinh(k2 * (x - x02)) /
                 (2. * (np.cosh(k2 * (x - x02)) + 1.)**2) )


def d2logistic2_offset_p(x, p):
    """
    Wrapper function for :func:`d2logistic2_offset`: `d2logistic2_offset(x, *p)`.
    """
    return d2logistic2_offset(x, *p)


# -----------------------------------------------------------

if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # logistic(1., 1., 0., 2.)
    # # 0.5
    # logistic(1., 1., 2., 1.)
    # # 0.5
    # logistic(2., 1., 1., 1.)
    # # 1./(1.+np.exp(-1.))
    # # 0.7310585786300049
    # logistic_p(1., [1., 0., 2.])
    # logistic_p(1., [1., 2., 1.])
    # logistic_p(2., [1., 1., 1.])
    # logistic_offset(1., 1., 0., 2., 1.)
    # logistic_offset(1., 1., 2., 1., 1.)
    # # 1.5
    # logistic_offset(2., 1., 1., 1., 1.)
    # # 1./(1.+np.exp(-1.)) + 1.
    # # 1.7310585786300049
    # logistic_offset_p(1., [1., 0., 2., 1.])
    # logistic_offset_p(1., [1., 2., 1., 1.])
    # logistic_offset_p(2., [1., 1., 1., 1.])
    # logistic2_offset(1.,  1., 2., 1.,  2., 2., 1.,  1.)
    # # 0.5
    # logistic2_offset_p(1.,  [1., 2., 1.,  2., 2., 1.,  1.])
    # dlogistic(1., 1., 2., 1.)
    # # 0.5
    # dlogistic_offset(1., 1., 2., 1., 1.)
    # # 0.5
    # dlogistic2_offset(1.,  1., 2., 1.,  2., 2., 1.,  1.)
    # # -0.5
    # print(np.around(d2logistic(1., 1., 2., 2.),4))
    # # 0.3199
    # print(np.around(d2logistic_offset(1., 1., 2., 2., 1.),4))
    # # 0.3199
    # print(np.around(d2logistic2_offset(1., 1., 2., 2.,  2., 2., 2.,  1.),4))
    # # -0.3199
