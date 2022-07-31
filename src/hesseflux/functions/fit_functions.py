#!/usr/bin/env python
"""
Module defines common functions that are used in curve_fit or fmin
parameter estimations.

For all fit functions, it defines the functions in two forms (ex. of 3
params):

    `func(x, p1, p2, p3)`

    `func_p(x, p)` with `p[0:3]`

The first form can be used, for example, with
`scipy.optimize.curve_fit` (ex. function f1x=a+b/x):

    `p, cov = scipy.optimize.curve_fit(functions.f1x, x, y, p0=[p0,p1])`

It also defines two cost functions along with the fit functions, one
with the absolute sum, one with the squared sum of the deviations:

    `cost_func`    `sum(abs(obs-func))`

    `cost2_func`   `sum((obs-func)**2)`

These cost functions can be used, for example, with
`scipy.optimize.minimize`:

    `p = scipy.optimize.minimize(jams.functions.cost_f1x, np.array([p1,p2]), args=(x,y), method='Nelder-Mead', options={'disp':False})`

Note the different argument orders:

    `curvefit` needs `f(x,*args)` with the independent variable as the
    first argument and the parameters to fit as separate remaining
    arguments.

    `minimize` is a general minimiser with respect to the first argument,
    i.e. `func(p,*args)`.

The module provides also two common cost functions (absolute and
squared deviations) where any function in the form `func(x, p)` can be
used as second argument:

    `cost_abs(p, func, x, y)`

    `cost_square(p, func, x, y)`

This means, for example `cost_f1x(p, x, y)` is the same as
`cost_abs(p, functions.f1x_p, x, y)`.
For example:

    `p = scipy.optimize.minimize(jams.functions.cost_abs, np.array([p1,p2]), args=(functions.f1x_p,x,y), method='Nelder-Mead', options={'disp':False})`

The current functions are (the functions have the name in the first
column. The seond form has a '\_p' appended to the name. The cost
functions, which have 'cost\_' and 'cost2\_' prepended to the name.):

    arrhenius         1 param:  Arrhenius temperature dependence of biochemical rates: `exp((T-TC25)*E/(T25*R*(T+T0)))`, parameter: E

    f1x               2 params: General 1/x function: `a + b/x`

    fexp              3 params: General exponential function: `a + b * exp(c*x)`

    gauss             2 params: Gauss function: `1/(sig*sqrt(2*pi)) *exp(-(x-mu)**2/(2*sig**2))`, parameter: mu, sig

    lasslop           6 params: Lasslop et al. (2010) a rectangular, hyperbolic light-response GPP with Lloyd & Taylor (1994) respiration and the maximum canopy uptake rate at light saturation decreases exponentially with VPD as in Koerner (1995)

    line0             1 params: Straight line: `a*x`

    line              2 params: Straight line: `a + b*x`

    lloyd_fix         2 params: Lloyd & Taylor (1994) Arrhenius type with T0=-46.02 degC and Tref=10 degC

    lloyd_only_rref   1 param:  Lloyd & Taylor (1994) Arrhenius type with fixed exponential term

    logistic          3 params: Logistic function: `a/(1+exp(-b(x-c)))`

    logistic_offset   4 params: Logistic function with offset: `a/(1+exp(-b(x-c))) + d`

    logistic2_offset  7 params: Double logistic function with offset `L1/(1+exp(-k1(x-x01))) - L2/(1+exp(-k2(x-x02))) + a`

    poly              n params: General polynomial: `c0 + c1*x + c2*x**2 + ... + cn*x**n`

    sabx              2 params: sqrt(f1x), i.e. general sqrt(1/x) function: `sqrt(a + b/x)`

    see               3 params: Sequential Elementary Effects fitting function: `a*(x-b)**c`

This module was written by Matthias Cuntz while at Department of
Computational Hydrosystems, Helmholtz Centre for Environmental
Research - UFZ, Leipzig, Germany, and continued while at Institut
National de Recherche pour l'Agriculture, l'Alimentation et
l'Environnement (INRAE), Nancy, France.

Copyright (c) 2012-2020 Matthias Cuntz - mc (at) macu (dot) de
Released under the MIT License; see LICENSE file for details.

* Written Dec 2012 by Matthias Cuntz (mc (at) macu (dot) de)
* Ported to Python 3, Feb 2013, Matthias Cuntz
* Added general cost functions cost_abs and cost_square, May 2013, Matthias Cuntz
* Added line0, Feb 2014, Matthias Cuntz
* Removed multiline_p, May 2020, Matthias Cuntz
* Changed to Sphinx docstring and numpydoc, May 2020, Matthias Cuntz

.. moduleauthor:: Matthias Cuntz

The following functions are provided:

.. autosummary::
   cost_abs
   cost_square
   arrhenius
   arrhenius_p
   cost_arrhenius
   cost2_arrhenius
   f1x
   f1x_p
   cost_f1x
   cost2_f1x
   fexp
   fexp_p
   cost_fexp
   cost2_fexp
   gauss
   gauss_p
   cost_gauss
   cost2_gauss
   lasslop
   lasslop_p
   cost_lasslop
   cost2_lasslop
   line
   line_p
   cost_line
   cost2_line
   line0
   line0_p
   cost_line0
   cost2_line0
   lloyd_fix
   lloyd_fix_p
   cost_lloyd_fix
   cost2_lloyd_fix
   lloyd_only_rref
   lloyd_only_rref_p
   cost_lloyd_only_rref
   cost2_lloyd_only_rref
   sabx
   sabx_p
   cost_sabx
   cost2_sabx
   poly
   poly_p
   cost_poly
   cost2_poly
   cost_logistic
   cost2_logistic
   cost_logistic_offset
   cost2_logistic_offset
   cost_logistic2_offset
   cost2_logistic2_offset
   see
   see_p
   cost_see
   cost2_see
"""
from __future__ import division, absolute_import, print_function
import numpy as np
import scipy.special as sp
try:        # import package
    from .logistic_function import logistic_p, logistic_offset_p, logistic2_offset_p
    from ..const import T0, T25, R
except:
    try:    # e.g. python nee2gpp.py
        from functions.logistic_function import logistic_p, logistic_offset_p, logistic2_offset_p
        from const import T0, T25, R
    except: # python fit_functions.py
        from logistic_function import logistic_p, logistic_offset_p, logistic2_offset_p
        T0  = 273.15    # Celcius <-> Kelvin [K]
        T25 = 298.15    # Standard ambient temperature [K]
        R   = 8.3144621 # Ideal gas constant [J K^-1 mol^-1]


__all__ = ['cost_abs', 'cost_square',
           'arrhenius', 'arrhenius_p', 'cost_arrhenius', 'cost2_arrhenius',
           'f1x', 'f1x_p', 'cost_f1x', 'cost2_f1x',
           'fexp', 'fexp_p', 'cost_fexp', 'cost2_fexp',
           'gauss', 'gauss_p', 'cost_gauss', 'cost2_gauss',
           'lasslop', 'lasslop_p', 'cost_lasslop', 'cost2_lasslop',
           'line', 'line_p', 'cost_line', 'cost2_line',
           'line0', 'line0_p', 'cost_line0', 'cost2_line0',
           'lloyd_fix', 'lloyd_fix_p', 'cost_lloyd_fix', 'cost2_lloyd_fix',
           'lloyd_only_rref', 'lloyd_only_rref_p', 'cost_lloyd_only_rref', 'cost2_lloyd_only_rref',
           'sabx', 'sabx_p', 'cost_sabx', 'cost2_sabx',
           'poly', 'poly_p', 'cost_poly', 'cost2_poly',
           'cost_logistic', 'cost2_logistic',
           'cost_logistic_offset', 'cost2_logistic_offset',
           'cost_logistic2_offset', 'cost2_logistic2_offset',
           'see', 'see_p', 'cost_see', 'cost2_see']


# -----------------------------------------------------------
# general cost functions
def cost_abs(p, func, x, y):
    """
    General cost function for robust optimising `func(x,p)` vs. `y` with sum of absolute deviations.

    Parameters
    ----------
    p : iterable of floats
        parameters
    func : callable
        `fun(x,p) -> float`
    x : float or array_like of floats
        independent variable
    y : float or array_like of floats
        dependent variable, observations

    Returns
    -------
    float
        sum of absolute deviations
    """
    return np.sum(np.abs(y-func(x,p)))


def cost_square(p, func, x, y):
    """
    General cost function for optimising `func(x,p)` vs. `y` with sum of square deviations.

    Parameters
    ----------
    p : iterable of floats
        parameters
    func : callable
        `fun(x,p) -> float`
    x : float or array_like of floats
        independent variable
    y : float or array_like of floats
        dependent variable, observations

    Returns
    -------
    float
        sum of squared deviations
    """
    return np.sum((y-func(x,p))**2)


# -----------------------------------------------------------
# arrhenius
def arrhenius(T, E):
    """
    Arrhenius temperature dependence of rates.

    Parameters
    ----------
    T : float or array_like of floats
        temperature [degC]
    E : float
        activation energy [J]

    Returns
    -------
    float
        function value(s)
    """
    return np.exp((T-(T25-T0))*E/(T25*R*(T+T0)))


def arrhenius_p(T, p):
    """
    Arrhenius temperature dependence of rates.

    Parameters
    ----------
    T : float or array_like of floats
        temperature [degC]
    p : iterable
        `p[0]` is activation energy [J]

    Returns
    -------
    float
        function value(s)
    """
    return np.exp((T-(T25-T0))*p[0]/(T25*R*(T+T0)))


def cost_arrhenius(p, T, rate):
    """
    Sum of absolute deviations of obs and arrhenius function.

    Parameters
    ----------
    p : iterable of floats
        `p[0]` is activation energy [J]
    x : float or array_like of floats
        independent variable
    y : float or array_like of floats
        dependent variable, observations

    Returns
    -------
    float
        sum of absolute deviations
    """
    return np.sum(np.abs(rate-arrhenius_p(T,p)))


def cost2_arrhenius(p, T, rate):
    """
    Sum of squared deviations of obs and arrhenius.

    Parameters
    ----------
    p : iterable of floats
        `p[0]` is activation energy [J]
    x : float or array_like of floats
        independent variable
    y : float or array_like of floats
        dependent variable, observations

    Returns
    -------
    float
        sum of squared deviations
    """
    return np.sum((rate-arrhenius_p(T,p))**2)


# -----------------------------------------------------------
# a+b/x
def f1x(x,a,b):
    """
    General 1/x function: a + b/x

    Parameters
    ----------
    x : float or array_like of floats
        independent variable
    a : float
        first parameter
    b : float
        second parameter

    Returns
    -------
    float
        function value(s)
    """
    return a+b/x


def f1x_p(x,p):
    """
    General 1/x function: a + b/x

    Parameters
    ----------
    x : float or array_like of floats
        independent variable
    p : iterable of floats
        parameters (`len(p)=2`)

        `p[0]` a

        `p[1]` b

    Returns
    -------
    float
        function value(s)
    """
    return p[0]+p[1]/x


def cost_f1x(p,x,y):
    """
    Sum of absolute deviations of obs and general 1/x function: a + b/x

    Parameters
    ----------
    p : iterable of floats
        parameters (`len(p)=2`)

        `p[0]` a

        `p[1]` b
    x : float or array_like of floats
        independent variable
    y : float or array_like of floats
        dependent variable, observations

    Returns
    -------
    float
        sum of absolute deviations
    """
    return np.sum(np.abs(y-f1x_p(x,p)))


def cost2_f1x(p,x,y):
    """
    Sum of squared deviations of obs and general 1/x function: a + b/x

    Parameters
    ----------
    p : iterable of floats
        parameters (`len(p)=2`)

        `p[0]` a

        `p[1]` b
    x : float or array_like of floats
        independent variable
    y : float or array_like of floats
        dependent variable, observations

    Returns
    -------
    float
        sum of squared deviations
    """
    return np.sum((y-f1x_p(x,p))**2)


# -----------------------------------------------------------
# a+b*exp(c*x)
def fexp(x,a,b,c):
    """
    General exponential function: a + b * exp(c*x)

    Parameters
    ----------
    x : float or array_like of floats
        independent variable
    a : float
        first parameter
    b : float
        second parameter
    c : float
        third parameter

    Returns
    -------
    float
        function value(s)
    """
    return a+b*np.exp(c*x)


def fexp_p(x,p):
    """
    General exponential function: a + b * exp(c*x)

    Parameters
    ----------
    x : float or array_like of floats
        independent variable
    p : iterable of floats
        parameters (`len(p)=3`)

        `p[0]` a

        `p[1]` b

        `p[2]` c

    Returns
    -------
    float
        function value(s)
    """
    return p[0]+p[1]*np.exp(p[2]*x)


def cost_fexp(p,x,y):
    """
    Sum of absolute deviations of obs and general exponential function: a + b * exp(c*x)

    Parameters
    ----------
    p : iterable of floats
        parameters (`len(p)=3`)

        `p[0]` a

        `p[1]` b

        `p[2]` c
    x : float or array_like of floats
        independent variable
    y : float or array_like of floats
        dependent variable, observations

    Returns
    -------
    float
        sum of absolute deviations
    """
    return np.sum(np.abs(y-fexp_p(x,p)))


def cost2_fexp(p,x,y):
    """
    Sum of squared deviations of obs and general exponential function: a + b * exp(c*x)

    Parameters
    ----------
    p : iterable of floats
        parameters (`len(p)=3`)

        `p[0]` a

        `p[1]` b

        `p[2]` c
    x : float or array_like of floats
        independent variable
    y : float or array_like of floats
        dependent variable, observations

    Returns
    -------
    float
        sum of squared deviations
    """
    return np.sum((y-fexp_p(x,p))**2)


# -----------------------------------------------------------
# Gauss: 1/(sig*sqrt(2*pi)) *exp(-(x-mu)**2/(2*sig**2))
def gauss(x,mu,sig):
    """
    Gauss function: 1 / (sqrt(2*pi)*sig) * exp( -(x-mu)**2 / (2*sig**2) )

    Parameters
    ----------
    x : float or array_like of floats
        independent variable
    mu : float
        mean
    sig : float
        width

    Returns
    -------
    float
        function value(s)
    """
    return np.exp(-(x-mu)**2/(2.*sig**2))/(sig*np.sqrt(2.*np.pi))


def gauss_p(x,p):
    """
    Gauss function: 1 / (sqrt(2*pi)*sig) * exp( -(x-mu)**2 / (2*sig**2) )

    Parameters
    ----------
    x : float or array_like of floats
        independent variable
    p : iterable of floats
        parameters (`len(p)=2`)

        `p[0]` mean

        `p[1]` width

    Returns
    -------
    float
        function value(s)
    """
    return np.exp(-(x-p[0])**2/(2.*p[1]**2))/(p[1]*np.sqrt(2.*np.pi))


def cost_gauss(p,x,y):
    """
    Sum of absolute deviations of obs and Gauss function: 1 / (sqrt(2*pi)*sig) * exp( -(x-mu)**2 / (2*sig**2) )

    Parameters
    ----------
    p : iterable of floats
        parameters (`len(p)=2`)

        `p[0]` mean

        `p[1]` width
    x : float or array_like of floats
        independent variable
    y : float or array_like of floats
        dependent variable, observations

    Returns
    -------
    float
        sum of absolute deviations
    """
    return np.sum(np.abs(y-gauss_p(x,p)))


def cost2_gauss(p,x,y):
    """
    Sum of squared deviations of obs and Gauss function: 1 / (sqrt(2*pi)*sig) * exp( -(x-mu)**2 / (2*sig**2) )

    Parameters
    ----------
    p : iterable of floats
        parameters (`len(p)=2`)

        `p[0]` mean

        `p[1]` width
    x : float or array_like of floats
        independent variable
    y : float or array_like of floats
        dependent variable, observations

    Returns
    -------
    float
        sum of squared deviations
    """
    return np.sum((y-gauss_p(x,p))**2)


# -----------------------------------------------------------
# lasslop
def lasslop(Rg, et, VPD, alpha, beta0, k, Rref):
    """
    Lasslop et al. (2010) is the rectangular, hyperbolic light-response
    of NEE as by Falge et al. (2001), where the respiration is calculated
    with Lloyd & Taylor (1994), and the maximum canopy uptake rate at
    light saturation decreases exponentially with VPD as in Koerner (1995).

    Parameters
    ----------
    Rg : float or array_like of floats
        Global radiation [W m-2]
    et : float or array_like of floats
        Exponential in Lloyd & Taylor: np.exp(E0*(1./(Tref-T0)-1./(T-T0))) []
    VPD : float or array_like of floats
        Vapour Pressure Deficit [Pa]
    alpha : float
        Light use efficiency, i.e. initial slope of light response curve [umol(C) J-1]
    beta0 : float
        Maximum CO2 uptake rate at VPD0=10 hPa [umol(C) m-2 s-1]
    k : float
        e-folding of exponential decrease of maximum CO2 uptake with VPD increase [Pa-1]
    Rref : float
        Respiration at Tref (10 degC) [umol(C) m-2 s-1]

    Returns
    -------
    float
        net ecosystem exchange [umol(CO2) m-2 s-1]
    """
    # Lloyd & Taylor (1994)
    gamma = Rref*et
    # Koerner (1995)
    VPD0  = 1000. # 10 hPa
    kk    = np.clip(-k*(VPD-VPD0), -600., 600.)
    beta  = np.where(VPD > VPD0, beta0*np.exp(kk), beta0)
    return -alpha*beta*Rg/(alpha*Rg+beta) + gamma


def lasslop_p(Rg, et, VPD, p):
    """
    Lasslop et al. (2010) is the rectangular, hyperbolic light-response
    of NEE as by Falge et al. (2001), where the respiration is calculated
    with Lloyd & Taylor (1994), and the maximum canopy uptake rate at
    light saturation decreases exponentially with VPD as in Koerner (1995).

    Parameters
    ----------
    Rg : float or array_like of floats
        Global radiation [W m-2]
    et : float or array_like of floats
        Exponential in Lloyd & Taylor: np.exp(E0*(1./(Tref-T0)-1./(T-T0))) []
    VPD : float or array_like of floats
        Vapour Pressure Deficit [Pa]
    p : iterable of floats
        parameters (`len(p)=4`)

        `p[0]` Light use efficiency, i.e. initial slope of light response curve [umol(C) J-1]

        `p[1]` Maximum CO2 uptake rate at VPD0=10 hPa [umol(C) m-2 s-1]

        `p[2]` e-folding of exponential decrease of maximum CO2 uptake with VPD increase [Pa-1]

        `p[3]` Respiration at Tref (10 degC) [umol(C) m-2 s-1]

    Returns
    -------
    float
        net ecosystem exchange [umol(CO2) m-2 s-1]
    """
    # Lloyd & Taylor (1994)
    gamma = p[3]*et
    # Koerner (1995)
    VPD0  = 1000. # 10 hPa
    kk    = np.clip(-p[2]*(VPD-VPD0), -600., 600.)
    beta  = np.where(VPD > VPD0, p[1]*np.exp(kk), p[1])
    return -p[0]*beta*Rg/(p[0]*Rg+beta) + gamma


def cost_lasslop(p, Rg, et, VPD, NEE):
    """
    Sum of absolute deviations of obs and Lasslop.

    Parameters
    ----------
    p : iterable of floats
        parameters (`len(p)=4`)

        `p[0]` Light use efficiency, i.e. initial slope of light response curve [umol(C) J-1]

        `p[1]` Maximum CO2 uptake rate at VPD0=10 hPa [umol(C) m-2 s-1]

        `p[2]` e-folding of exponential decrease of maximum CO2 uptake with VPD increase [Pa-1]

        `p[3]` Respiration at Tref (10 degC) [umol(C) m-2 s-1]
    Rg : float or array_like of floats
        Global radiation [W m-2]
    et : float or array_like of floats
        Exponential in Lloyd & Taylor: np.exp(E0*(1./(Tref-T0)-1./(T-T0))) []
    VPD : float or array_like of floats
        Vapour Pressure Deficit [Pa]
    NEE : float or array_like of floats
        Observed net ecosystem exchange [umol(CO2) m-2 s-1]

    Returns
    -------
    float
        sum of absolute deviations
    """
    return np.sum(np.abs(NEE-lasslop(Rg, et, VPD, p[0], p[1], p[2], p[3])))


def cost2_lasslop(p, Rg, et, VPD, NEE):
    """
    Sum of squared deviations of obs and Lasslop.

    Parameters
    ----------
    p : iterable of floats
        parameters (`len(p)=4`)

        `p[0]` Light use efficiency, i.e. initial slope of light response curve [umol(C) J-1]

        `p[1]` Maximum CO2 uptake rate at VPD0=10 hPa [umol(C) m-2 s-1]

        `p[2]` e-folding of exponential decrease of maximum CO2 uptake with VPD increase [Pa-1]

        `p[3]` Respiration at Tref (10 degC) [umol(C) m-2 s-1]
    Rg : float or array_like of floats
        Global radiation [W m-2]
    et : float or array_like of floats
        Exponential in Lloyd & Taylor: np.exp(E0*(1./(Tref-T0)-1./(T-T0))) []
    VPD : float or array_like of floats
        Vapour Pressure Deficit [Pa]
    NEE : float or array_like of floats
        Observed net ecosystem exchange [umol(CO2) m-2 s-1]

    Returns
    -------
    float
        sum of squared deviations
    """
    return np.sum((NEE-lasslop(Rg, et, VPD, p[0], p[1], p[2], p[3]))**2)


# -----------------------------------------------------------
# a+b*x
def line(x,a,b):
    """
    Straight line: a + b*x

    Parameters
    ----------
    x : float or array_like of floats
        independent variable
    a : float
        first parameter
    b : float
        second parameter

    Returns
    -------
    float
        function value(s)
    """
    return a+b*x


def line_p(x,p):
    """
    Straight line: a + b*x

    Parameters
    ----------
    x : float or array_like of floats
        independent variable
    p : iterable of floats
        parameters (`len(p)=2`)

        `p[0]` a

        `p[1]` b

    Returns
    -------
    float
        function value(s)
    """
    return p[0]+p[1]*x


def cost_line(p,x,y):
    """
    Sum of absolute deviations of obs and straight line: a + b*x

    Parameters
    ----------
    p : iterable of floats
        parameters (`len(p)=2`)

        `p[0]` a

        `p[1]` b
    x : float or array_like of floats
        independent variable
    y : float or array_like of floats
        dependent variable, observations

    Returns
    -------
    float
        sum of absolute deviations
    """
    return np.sum(np.abs(y-line_p(x,p)))


def cost2_line(p,x,y):
    """
    Sum of squared deviations of obs and straight line: a + b*x

    Parameters
    ----------
    p : iterable of floats
        parameters (`len(p)=2`)

        `p[0]` a

        `p[1]` b
    x : float or array_like of floats
        independent variable
    y : float or array_like of floats
        dependent variable, observations

    Returns
    -------
    float
        sum of squared deviations
    """
    return np.sum((y-line_p(x,p))**2)


# -----------------------------------------------------------
# b*x
def line0(x,a):
    """
    Straight line through origin: a*x

    Parameters
    ----------
    x : float or array_like of floats
        independent variable
    a : float
        first parameter

    Returns
    -------
    float
        function value(s)
    """
    return a*x


def line0_p(x,p):
    """
    Straight line through origin: a*x

    Parameters
    ----------
    x : float or array_like of floats
        independent variable
    p : iterable of floats
        `p[0]` is a

    Returns
    -------
    float
        function value(s)
    """
    return p*x

  
def cost_line0(p,x,y):
    """
    Sum of absolute deviations of obs and straight line through origin: a*x

    Parameters
    ----------
    p : iterable of floats
        `p[0]` is a
    x : float or array_like of floats
        independent variable
    y : float or array_like of floats
        dependent variable, observations

    Returns
    -------
    float
        sum of absolute deviations
    """
    return np.sum(np.abs(y-line0_p(x,p)))

  
def cost2_line0(p,x,y):
    """
    Sum of squared deviations of obs and straight line through origin: a*x

    Parameters
    ----------
    p : iterable of floats
        `p[0]` is a
    x : float or array_like of floats
        independent variable
    y : float or array_like of floats
        dependent variable, observations

    Returns
    -------
    float
        sum of squared deviations
    """
    return np.sum((y-line0_p(x,p))**2)


# -----------------------------------------------------------
# lloyd_fix
def lloyd_fix(T, Rref, E0):
    """
    Lloyd & Taylor (1994) Arrhenius type with T0=-46.02 degC and Tref=10 degC

    Parameters
    ----------
    T : float or array_like of floats
        Temperature [K]
    Rref : float
        Respiration at Tref=10 degC [umol(C) m-2 s-1]
    E0 : float
        Activation energy [K]

    Returns
    -------
    float
        Respiration [umol(C) m-2 s-1]
    """
    Tref = 283.15 #  10    [degC]
    T0   = 227.13 # -46.02 [degC]
    return Rref*np.exp(E0*(1./(Tref-T0)-1./(T-T0)))


def lloyd_fix_p(T, p):
    """
    Lloyd & Taylor (1994) Arrhenius type with T0=-46.02 degC and Tref=10 degC

    Parameters
    ----------
    T : float or array_like of floats
        Temperature [K]
    p : iterable of floats
        parameters (`len(p)=2`)

        `p[0]` Respiration at Tref=10 degC [umol(C) m-2 s-1]

        `p[1]` Activation energy [K]

    Returns
    -------
    float
        Respiration [umol(C) m-2 s-1]
    """
    Tref = 283.15 #  10    [degC]
    T0   = 227.13 # -46.02 [degC]
    return p[0]*np.exp(p[1]*(1./(Tref-T0)-1./(T-T0)))


def cost_lloyd_fix(p, T, resp):
    """
    Sum of absolute deviations of obs and Lloyd & Taylor (1994) Arrhenius type.

    Parameters
    ----------
    p : iterable of floats
        parameters (`len(p)=2`)

        `p[0]` Respiration at Tref=10 degC [umol(C) m-2 s-1]

        `p[1]` Activation energy [K]
    T : float or array_like of floats
        Temperature [K]
    resp : float or array_like of floats
        Observed respiration [umol(C) m-2 s-1]

    Returns
    -------
    float
        sum of absolute deviations
    """
    return np.sum(np.abs(resp-lloyd_fix_p(T,p)))


def cost2_lloyd_fix(p, T, resp):
    """
    Sum of squared deviations of obs and Lloyd & Taylor (1994) Arrhenius type.

    Parameters
    ----------
    p : iterable of floats
        parameters (`len(p)=2`)

        `p[0]` Respiration at Tref=10 degC [umol(C) m-2 s-1]

        `p[1]` Activation energy [K]
    T : float or array_like of floats
        Temperature [K]
    resp : float or array_like of floats
        Observed respiration [umol(C) m-2 s-1]

    Returns
    -------
    float
        sum of squared deviations
    """
    return np.sum((resp-lloyd_fix_p(T,p))**2)


# -----------------------------------------------------------
# lloyd_only_rref
def lloyd_only_rref(et, Rref):
    """
    If E0 is know in Lloyd & Taylor (1994) then one can calc
    the exponential term outside the routine and the fitting
    becomes linear. One could also use functions.line0.

    Parameters
    ----------
    et : float or array_like of floats
        exp-term in Lloyd & Taylor
    Rref : float
        Respiration at Tref=10 degC [umol(C) m-2 s-1]

    Returns
    -------
    float
        Respiration [umol(C) m-2 s-1]
    """
    return Rref*et


def lloyd_only_rref_p(et, p):
    """
    If E0 is know in Lloyd & Taylor (1994) then one can calc
    the exponential term outside the routine and the fitting
    becomes linear. One could also use functions.line0.

    Parameters
    ----------
    et : float or array_like of floats
        exp-term in Lloyd & Taylor
    p : iterable of floats
        `p[0]` is respiration at Tref=10 degC [umol(C) m-2 s-1]

    Returns
    -------
    float
        Respiration [umol(C) m-2 s-1]
    """
    return p[0]*et


def cost_lloyd_only_rref(p, et, resp):
    """
    Sum of absolute deviations of obs and Lloyd & Taylor with known exponential term.

    Parameters
    ----------
    p : iterable of floats
        `p[0]` is respiration at Tref=10 degC [umol(C) m-2 s-1]
    et : float or array_like of floats
        exp-term in Lloyd & Taylor
    resp : float or array_like of floats
        Observed respiration [umol(C) m-2 s-1]

    Returns
    -------
    float
        sum of absolute deviations
    """
    return np.sum(np.abs(resp-lloyd_only_rref_p(et,p)))


def cost2_lloyd_only_rref(p, et, resp):
    """
    Sum of squared deviations of obs and Lloyd & Taylor with known exponential term.

    Parameters
    ----------
    p : iterable of floats
        `p[0]` is respiration at Tref=10 degC [umol(C) m-2 s-1]
    et : float or array_like of floats
        exp-term in Lloyd & Taylor
    resp : float or array_like of floats
        Observed respiration [umol(C) m-2 s-1]

    Returns
    -------
    float
        sum of squared deviations
    """
    return np.sum((resp-lloyd_only_rref_p(et,p))**2)


# -----------------------------------------------------------
# sqrt(a + b/x) - theoretical form of Jackknife-after-bootstrap

def sabx(x, a, b):
    """
    Square root of general 1/x function: sqrt(a + b/x)

    Parameters
    ----------
    x : float or array_like of floats
        independent variable
    a : float
        first parameter
    b : float
        second parameter

    Returns
    -------
    float
        function value(s)
    """
    return np.sqrt(a+b/x)


def sabx_p(x, p):
    """
    Square root of general 1/x function: sqrt(a + b/x)

    Parameters
    ----------
    x : float or array_like of floats
        independent variable
    p : iterable of floats
        parameters (`len(p)=2`)

        `p[0]` a

        `p[1]` b

    Returns
    -------
    float
        function value(s)
    """
    return np.sqrt(p[0]+p[1]/x)

  
def cost_sabx(p,x,y):
    """
    Sum of absolute deviations of obs and square root of general 1/x function: sqrt(a + b/x)

    Parameters
    ----------
    p : iterable of floats
        parameters (`len(p)=2`)

        `p[0]` a

        `p[1]` b
    x : float or array_like of floats
        independent variable
    y : float or array_like of floats
        dependent variable, observations

    Returns
    -------
    float
        sum of absolute deviations
    """
    return np.sum(np.abs(y-sabx_p(x,p)))

def cost2_sabx(p,x,y):
    """
    Sum of squared deviations of obs and square root of general 1/x function: sqrt(a + b/x)

    Parameters
    ----------
    p : iterable of floats
        parameters (`len(p)=2`)

        `p[0]` a

        `p[1]` b
    x : float or array_like of floats
        independent variable
    y : float or array_like of floats
        dependent variable, observations

    Returns
    -------
    float
        sum of squared deviations
    """
    return np.sum((y-sabx_p(x,p))**2)


# -----------------------------------------------------------
# c0 + c1*x + c2*x**2 + ... + cn*x**n
def poly(x,*args):
    """
    General polynomial: c0 + c1*x + c2*x**2 + ... + cn*x**n

    Parameters
    ----------
    x : float or array_like of floats
        independent variable
    *args : float
        parameters `len(args)=n+1`

    Returns
    -------
    float
        function value(s)
    """
    return np.polynomial.polynomial.polyval(x, list(args))

  
def poly_p(x,p):
    """
    General polynomial: c0 + c1*x + c2*x**2 + ... + cn*x**n

    Parameters
    ----------
    x : float or array_like of floats
        independent variable
    p : iterable of floats
        parameters (`len(p)=n+1`)

    Returns
    -------
    float
        function value(s)
    """
    return np.polynomial.polynomial.polyval(x, p)

  
def cost_poly(p,x,y):
    """
    Sum of absolute deviations of obs and general polynomial: c0 + c1*x + c2*x**2 + ... + cn*x**n

    Parameters
    ----------
    p : iterable of floats
        parameters (`len(p)=n+1`)
    x : float or array_like of floats
        independent variable
    y : float or array_like of floats
        dependent variable, observations

    Returns
    -------
    float
        sum of absolute deviations
    """
    return np.sum(np.abs(y-poly_p(x,p)))


def cost2_poly(p,x,y):
    """
    Sum of squared deviations of obs and general polynomial: c0 + c1*x + c2*x**2 + ... + cn*x**n

    Parameters
    ----------
    p : iterable of floats
        parameters (`len(p)=n+1`)
    x : float or array_like of floats
        independent variable
    y : float or array_like of floats
        dependent variable, observations

    Returns
    -------
    float
        sum of squared deviations
    """
    return np.sum((y-poly_p(x,p))**2)


# -----------------------------------------------------------
# a/(1+exp(-b(x-c))) - logistic function
def cost_logistic(p, x, y):
    """
    Sum of absolute deviations of obs and logistic function L/(1+exp(-k(x-x0)))

    Parameters
    ----------
    p : iterable of floats
        parameters (`len(p)=3`)
        
        `p[0]` L  - Maximum of logistic function

        `p[1]` k  - Steepness of logistic function

        `p[2]` x0 - Inflection point of logistic function
    x : float or array_like of floats
        independent variable
    y : float or array_like of floats
        dependent variable, observations

    Returns
    -------
    float
        sum of absolute deviations
    """
    return np.sum(np.abs(y-logistic_p(x,p)))


def cost2_logistic(p,x,y):
    """
    Sum of squared deviations of obs and logistic function L/(1+exp(-k(x-x0)))

    Parameters
    ----------
    p : iterable of floats
        parameters (`len(p)=3`)
        
        `p[0]` L  - Maximum of logistic function

        `p[1]` k  - Steepness of logistic function

        `p[2]` x0 - Inflection point of logistic function
    x : float or array_like of floats
        independent variable
    y : float or array_like of floats
        dependent variable, observations

    Returns
    -------
    float
        sum of squared deviations
    """
    return np.sum((y-logistic_p(x,p))**2)


# -----------------------------------------------------------
# a/(1+exp(-b(x-c))) + d - logistic function with offset
def cost_logistic_offset(p, x, y):
    """
    Sum of absolute deviations of obs and logistic function 1/x function: L/(1+exp(-k(x-x0))) + a

    Parameters
    ----------
    p : iterable of floats
        parameters (`len(p)=4`)
        
        `p[0]` L  - Maximum of logistic function

        `p[1]` k  - Steepness of logistic function

        `p[2]` x0 - Inflection point of logistic function

        `p[3]` a  - Offset of logistic function
    x : float or array_like of floats
        independent variable
    y : float or array_like of floats
        dependent variable, observations

    Returns
    -------
    float
        sum of absolute deviations
    """
    return np.sum(np.abs(y-logistic_offset_p(x,p)))


def cost2_logistic_offset(p,x,y):
    """
    Sum of squared deviations of obs and logistic function 1/x function: L/(1+exp(-k(x-x0))) + a

    Parameters
    ----------
    p : iterable of floats
        parameters (`len(p)=4`)
        
        `p[0]` L  - Maximum of logistic function

        `p[1]` k  - Steepness of logistic function

        `p[2]` x0 - Inflection point of logistic function

        `p[3]` a  - Offset of logistic function
    x : float or array_like of floats
        independent variable
    y : float or array_like of floats
        dependent variable, observations

    Returns
    -------
    float
        sum of squared deviations
    """
    return np.sum((y-logistic_offset_p(x,p))**2)


# -----------------------------------------------------------
# L1/(1+exp(-k1(x-x01))) - L2/(1+exp(-k2(x-x02))) + a2 - double logistic function with offset
def cost_logistic2_offset(p, x, y):
    """
    Sum of absolute deviations of obs and double logistic function with offset:
    L1/(1+exp(-k1(x-x01))) - L2/(1+exp(-k2(x-x02))) + a

    Parameters
    ----------
    p : iterable of floats
        parameters (`len(p)=7`)
        
        `p[0]` L1  - Maximum of first logistic function

        `p[1]` k1  - Steepness of first logistic function

        `p[2]` x01 - Inflection point of first logistic function
        
        `p[3]` L2  - Maximum of second logistic function

        `p[4]` k2  - Steepness of second logistic function

        `p[5]` x02 - Inflection point of second logistic function

        `p[6]` a   - Offset of double logistic function
    x : float or array_like of floats
        independent variable
    y : float or array_like of floats
        dependent variable, observations

    Returns
    -------
    float
        sum of absolute deviations
    """
    return np.sum(np.abs(y-logistic2_offset_p(x,p)))


def cost2_logistic2_offset(p,x,y):
    """
    Sum of squared deviations of obs and double logistic function with offset:
    L1/(1+exp(-k1(x-x01))) - L2/(1+exp(-k2(x-x02))) + a

    Parameters
    ----------
    p : iterable of floats
        parameters (`len(p)=7`)
        
        `p[0]` L1  - Maximum of first logistic function

        `p[1]` k1  - Steepness of first logistic function

        `p[2]` x01 - Inflection point of first logistic function
        
        `p[3]` L2  - Maximum of second logistic function

        `p[4]` k2  - Steepness of second logistic function

        `p[5]` x02 - Inflection point of second logistic function

        `p[6]` a   - Offset of double logistic function
    x : float or array_like of floats
        independent variable
    y : float or array_like of floats
        dependent variable, observations

    Returns
    -------
    float
        sum of squared deviations
    """
    return np.sum((y-logistic2_offset_p(x,p))**2)


# -----------------------------------------------------------
# a*(x-b)**c - Sequential Elementary Effects fitting function
def see(x, a, b, c):
    """
    Fit function of Sequential Elementary Effects: a * (x-b)**c

    Parameters
    ----------
    x : float or array_like of floats
        independent variable
    a : float
        first parameter
    b : float
        second parameter
    c : float
        third parameter

    Returns
    -------
    float
        function value(s)
    """
    return np.where((x-b)<0., 0., a*(x-b)**c)


def see_p(x, p):
    """
    Fit function of Sequential Elementary Effects: a * (x-b)**c

    Parameters
    ----------
    x : float or array_like of floats
        independent variable
    p : iterable of floats
        parameters (`len(p)=3`)

        `p[0]` a

        `p[1]` b

        `p[2]` c

    Returns
    -------
    float
        function value(s)
    """
    return np.where((x-p[1]) < 0., 0., p[0] * (x-p[1])**p[2])

  
def cost_see(p, x, y):
    """
    Sum of absolute deviations of obs and fit function of Sequential Elementary Effects: a * (x-b)**c

    Parameters
    ----------
    p : iterable of floats
        parameters (`len(p)=3`)

        `p[0]` a

        `p[1]` b

        `p[2]` c
    x : float or array_like of floats
        independent variable
    y : float or array_like of floats
        dependent variable, observations

    Returns
    -------
    float
        sum of absolute deviations
    """
    return np.sum(np.abs(y-see_p(x,p)))


def cost2_see(p,x,y):
    """
    Sum of squared deviations of obs and fit function of Sequential Elementary Effects: a * (x-b)**c

    Parameters
    ----------
    p : iterable of floats
        parameters (`len(p)=3`)

        `p[0]` a

        `p[1]` b

        `p[2]` c
    x : float or array_like of floats
        independent variable
    y : float or array_like of floats
        dependent variable, observations

    Returns
    -------
    float
        sum of squared deviations
    """
    return np.sum((y-see_p(x,p))**2)


# -----------------------------------------------------------

if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)

    # Rref = 1.0
    # E0   = 126.
    # T    = 293.15
    # resp = 2.0
    # print(lloyd_fix(T, Rref, E0))
    # #1.40590910521
    # print(lloyd_fix_p(T, [Rref, E0]))
    # #1.40590910521
    # print(cost_lloyd_fix([Rref, E0], T, resp))
    # #0.59409089479
    # print(cost2_lloyd_fix([Rref, E0], T, resp))
    # #0.352943991272

    # print(poly(T,2,1))
    # #295.15
    # print(poly_p(T,[2,1]))
    # #295.15
