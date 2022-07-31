#!/usr/bin/env python
"""
Module provides common test functions for parameter estimation and
optimisation algorithms such as Rosenbrock and Griewank functions.

Current functions are:
    ackley                >=2 params:     Ackley function, global optimum: 0.0 at origin

    goldstein_price       2 params:       Goldstein-Price function, global optimum: 3.0 (0.0,-1.0)

    griewank              2 or 10 params: Griewank function, global optimum: 0 at origin

    rastrigin             2 params:       Rastrigin function, global optimum: -2 (0,0)

    rosenbrock            2 params:       Rosenbrock function, global optimum: 0 (1,1)

    six_hump_camelback    2 params:       Six-hump Camelback function
                                          True Optima: -1.031628453489877
                                          (-0.08983,0.7126) and (0.08983,-0.7126)

This module was written by Matthias Cuntz while at Department of
Computational Hydrosystems, Helmholtz Centre for Environmental
Research - UFZ, Leipzig, Germany, and continued while at Institut
National de Recherche pour l'Agriculture, l'Alimentation et
l'Environnement (INRAE), Nancy, France.

Copyright (c) 2013-2020 Matthias Cuntz - mc (at) macu (dot) de
Released under the MIT License; see LICENSE file for details.

* Written Oct 2013 by Matthias Cuntz (mc (at) macu (dot) de)
* Rearrange function library, Mar 2015, Matthias Cuntz
* Changed to Sphinx docstring and numpydoc, May 2020, Matthias Cuntz

.. moduleauthor:: Matthias Cuntz

The following functions are provided:

.. autosummary::
   ackley
   griewank
   goldstein_price
   rastrigin
   rosenbrock
   six_hump_camelback
"""
from __future__ import division, absolute_import, print_function
import numpy as np


__all__ = ['ackley', 'griewank', 'goldstein_price', 'rastrigin', 'rosenbrock', 'six_hump_camelback']


# -----------------------------------------------------------

def ackley(x):
    """
    Ackley function (>= 2-D).

    Global Optimum: 0.0, at origin.

    Parameters
    ----------
    x : array
        multi-dimensional x-values (len(x) >= 2)

    Returns
    -------
    float
       Value of Ackley function.
    """
    a = 20.0
    b = 0.2
    c = 2.0*np.pi

    n  = np.size(x)
    s1 = np.sum(x**2)
    s2 = np.sum(np.cos(c*x))
    f  = -a * np.exp(-b*np.sqrt(1.0/n*s1)) - np.exp(1.0/n*s2) + a + np.exp(1.0)

    return f

# -----------------------------------------------------------

def griewank(x):
    """
    Griewank function (2-D or 10-D).

    Global Optimum: 0.0, at origin.

    Parameters
    ----------
    x : array
        multi-dimensional x-values.

        `len(x)=2` or `len(x)=10`.
        
        `x[i]` bound to [-600,600] for all i.

    Returns
    -------
    float
       Value of Griewank function.
    """
    nopt = np.size(x)
    #if (nopt == 2) | (nopt == 10):
    xx = x
    if nopt==2:
        d = 200.0
    else:
        d = 4000.0

    u1 = 0.0
    u2 = 1.0
    for j in range(nopt):
        u1 = u1 + xx[j]**2/d
        u2 = u2 * np.cos(xx[j]/np.sqrt(float(j+1)))

    f = u1 - u2 + 1
    return f

# -----------------------------------------------------------

def goldstein_price(x):
    """
    Goldstein-Price function (2-D).

    Global Optimum: 3.0, at (0.0,-1.0).

    Parameters
    ----------
    x : array
        2 x-values. `len(x)=2`.
        
        `x[i]` bound to [-2,2] for i=1 and 2.

    Returns
    -------
    float
       Value of Goldstein-Price function.
    """
    x1 = x[0]
    x2 = x[1]
    u1 = (x1 + x2 + 1.0)**2
    u2 = 19. - 14.*x1 + 3.*x1**2 - 14.*x2 + 6.*x1*x2 +3.*x2**2
    u3 = (2.*x1 - 3.*x2)**2
    u4 = 18. - 32.*x1 + 12.*x1**2 + 48.*x2 -36.*x1*x2 + 27.*x2**2
    u5 = u1 * u2
    u6 = u3 * u4
    f = (1. + u5) * (30. + u6)
    return f

# -----------------------------------------------------------

def rastrigin(x):
    """
    Rastrigin function (2-D).

    Global Optimum: -2.0, at origin.

    Parameters
    ----------
    x : array
        2 x-values. `len(x)=2`.
        
        `x[i]` bound to [-1,1] for i=1 and 2.

    Returns
    -------
    float
       Value of Rastrigin function.
    """
    x1 = x[0]
    x2 = x[1]
    f = x1**2 + x2**2 - np.cos(18.0*x1) - np.cos(18.0*x2)
    return f

# -----------------------------------------------------------

def rosenbrock(x):
    """
    Rosenbrock function (2-D).

    Global Optimum: 0.0, at (1.0,1.0).

    Parameters
    ----------
    x : array
        2 x-values. `len(x)=2`.
        
        `x[1]` bound to [-5,5].

        `x[2]` bound to [-2,8].

    Returns
    -------
    float
       Value of Rosenbrock function.
    """

    x1 = x[0]
    x2 = x[1]
    a = 100.0
    f = a * (x2 - x1**2)**2 + (1 - x1)**2
    return f

# -----------------------------------------------------------

def six_hump_camelback(x):
    """
    Six-hump Camelback function (2-D).

    Global Optima: -1.031628453489877, at (-0.08983,0.7126) and (0.08983,-0.7126).

    Parameters
    ----------
    x : array
        2 x-values. `len(x)=2`.
        
        `x[i]` bound to [-5,5] for i=1 and 2.

    Returns
    -------
    float
       Value of Six-hump Camelback function.
    """
    x1 = x[0]
    x2 = x[1]
    f = (4 - 2.1*x1**2 + x1**4/3)*x1**2 + x1*x2 + (-4 + 4*x2**2)*x2**2
    return f

# -----------------------------------------------------------

if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
