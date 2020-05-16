#!/usr/bin/env python
"""
logtools is a Python port of the Control File Functions of
Logtools, the Logger Tools Software of Olaf Kolle, MPI-BGC Jena, (c) 2012.

Some functions are renamed compared to the original logger tools:

    `chs` -> `varchs`

    `add` -> `varadd`

    `sub` -> `varsub`

    `mul` -> `varmul`

    `div` -> `vardiv`

    `sqr` -> `varsqr`/`varsqrt`

    `exp` -> `varexp`

    `log` -> `varlog`

    `pot` -> `varpot`

Not all functions are implemented (yet). Missing functions are:

    `varset`

    `met_torad`

    `met_psy_rh`

    `met_dpt_rh`

    `write`

Some functions are slightly enhanced, which is reflected in the
documentation of the indidual functions.

All functions have an additional keyword `undef`, which defaults to -9999.:
elements are excluded from the calculations if any of the inputs equals `undef`.

Only bit_test and the if-statements `ifeq`, `ifne`, `ifle`, `ifge`, `iflt`, `igt`
do not have the `undef` keyword.

The Looger Tools control functions are:

    1. Assignment # not implemented


    2. Change sign
        x = varchs(a) means x = -a, where a is a variable or a number.

        def varchs(var1, undef=-9999.):


    3. Addition
        x = varadd(a, b) means x = a + b, where a and b are ndarray or float.

        def varadd(var1, var2, undef=-9999.):


    4. Subtraction
        x = varsub(a, b) means x = a - b, where a and b are ndarray or float.

        def varsub(var1, var2, undef=-9999.):


    5. Multiplication
        x = varmul(a, b) means x = a * b, where a and b are ndarray or float.

        def varmul(var1, var2, undef=-9999.):


    6. Division
        x = vardiv(a, b) means x = a/b, where a and b are ndarray or float.

        def vardiv(var1, var2, undef=-9999.):


    7. Square root
        x = varsqr(a) means x = sqrt(a), where a is a variable or a number.
        x = varsqrt(a) means x = sqrt(a), where a is a variable or a number.

        def varsqr(var1, undef=-9999.):
        def varsqrt(var1, undef=-9999.):


    8. Exponentiation of e
        x = varexp(a) means x = exp(a), where a is a variable or a number.

        def varexp(var1, undef=-9999.):


    9. Natural logarithm
        x = varlog(a) means x = ln(a), where a is a variable or a number.

        def varlog(var1, undef=-9999.):


    10. Exponentiation
        x = varpot(a, b) means x = a**b, where a and b are ndarray or float.

        def varpot(var1, var2, undef=-9999.):


    11. Apply linear function
        x = lin(y, a0, a1) means x = a0 + a1 * y,
        where a0 and a1 are ndarray or float.

        def lin(var1, a, b, undef=-9999.):


    12. Apply 2nd order function
        x=quad(y,a0,a1,a2) means x = a0 +a1*y + a2*y**2,
        where a0, a1 and a2 are ndarray or float.

        def quad(var1, a, b, c, undef=-9999.):


    13. Apply 3rd order function
        x=cubic(y,a0,a1,a2,a3) means x = a0 +a1*y+a2*y**2+a3*y**3,
        where a0, a1, a2 and a3 are ndarray or float.

        def cubic(var1, a, b, c, d, undef=-9999.):


    14. Calculate fraction of day from hours, minutes and seconds
        x = hms(h, m, s) means x = (h + m/60 + s/3600)/24,
        where h, m and s (hours, minutes and seconds) are ndarray or float.

        def hms(h, m, s, undef=-9999.):


    15. Bitwise test
        x = bit_test(y, b, start=0) means x = 1 if bit b ist set in y otherwise x = 0.
        Returns a list of b is an array.
        Counting of b starts at start.
        For the behaviour of the original logger tools, set start=1.
        Negative b's is not implemented.

        def bit_test(var1, var2, start=0):


    16. Replacement of underflows by new value
        x = setlow(y,lo,ln=None) means IF (y > lo) THEN x = ln ELSE x = y,
        where lo and ln are ndarray or float.
        ln is optional. If not given lo will be used.
        This function may be used to adjust small negative values of short wave radiation
        during nighttime to zero values.

        def setlow(dat, low, islow=None, undef=-9999.):


    17. Replacement of overflows by new value
        x = sethigh(y,lo,ln=None) means IF (y < lo) THEN x = ln ELSE x = y,
        where lo and ln are ndarray or float.
        ln is optional. If not given lo will be used.
        This function may be used to adjust relative humidity values of a little bit more than 100 % to 100 %.

        def sethigh(dat, high, ishigh=None, undef=-9999.):


    18. Replacement of underflows or overflows by the undef
        x = limits(y, ll, lh) means
        IF (y > ll) OR (y < lh) THEN x = undef ELSE x = y,
        where ll and lh are ndarray or float.
        This function may be used to check values lying in between certain limits.
        If one of the limits is exceeded the value is set to undef.

        def limits(dat, mini, maxi, undef=-9999.):


    19. Calculation of mean value
        x = mean(y1, y2, ..., yn) means x = (y1 + y2 + ... + yn)/n,
        where y1, y2, ..., yn are ndarray or float.

        def mean(var1, axis=None, undef=-9999.):


    20. Calculation of minimum value
        x = mini(y1,y2,...,yn) means x = min(y1,y2,...,yn),
        where y1, y2, ..., yn are ndarray or float.

        def mini(var1, axis=None, undef=-9999.):


    21. Calculation of maximum value
        x = maxi(y1,y2,...,yn) means x = max(y1,y2,...,yn),
        where y1, y2, ..., yn are ndarray or float.

        def maxi(var1, axis=None, undef=-9999.):


    22. Calculation of total radiation from net radiometer # no implemented


    23. Calculation of long wave radiation from net radiometer
        x = met_lwrad(y, Tp) where
        y is the output voltage of the net radiometer in mV,
        Tp is the temperature of the net radiometer body in degC.
        The total radiation in W m-2 is calculated according to the following formula:
        x=y*fl +sigma*(Tp +273.16)**4
        where sigma = 5.67051 * 10**8 W m-2 K-4 is the Stephan-Boltzmann-Constant and
        fl is the factor for long wave radiation (reciprocal value of sensitivity) in W m-2 per mV.
        The function assumes that fl was already applied before.
        All parameters may be ndarray or float.

        def met_lwrad(dat, tpyr, undef=-9999.): # assumes that dat was already multiplied with calibration factor


    24. Calculation of radiation temperature from long wave radiation
        x = met_trad(Rl, epsilon) where
        Rl is the long wave radiation in W m-2,
        epsilon is the long wave emissivity of the surface (between 0 and 1).
        The radiation temperature in degC is calculated according to the following formula:
        x= sqrt4(Rl/(sigma*epsilon)) - 273.16
        where sigma = 5.67051 * 108 W m-2 K-4 is the Stephan-Boltzmann-Constant.
        Both parameters may be ndarray or float.

        def met_trad(dat, eps, undef=-9999.):


    25. Calculation of albedo from short wave downward and upward radiation
        x = met_alb(Rsd, Rsu) where
        Rsd is the short wave downward radiation in Wm-2, Rsu is the short wave upward radiation in Wm-2,
        The albedo in % is calculated according to the following formula:
        x = 100 * ( Rsu / Rsd )
        If Rsd > 50 W m-2 or Rsu > 10 W m-2 the result is undef.
        Both parameters may be ndarray or float.

        def met_alb(swd, swu, swdmin=50., swumin=10., undef=-9999.):


    26. Calculation of albedo from short wave downward and upward radiation with limits
        x=met_albl(Rsd,Rsu,Rsd_limit,Rsu_limit)where
        Rsd is the short wave downward radiation in Wm-2,
        Rsu is the short wave upward radiation in Wm-2,
        Rsd_limit is the short wave downward radiation limit in Wm-2,
        Rsu_limit is the short wave upward radiation limit in Wm-2,
        The albedo in % is calculated according to the following formula:
        x = 100 * ( Rsu / Rsd )
        If Rsd > Rsd_limit or Rsu > Rsu_limit the result is undef.
        All four parameters may be ndarray or float.

        def met_albl(swd, swu, swdmin, swumin, undef=-9999.):


    27. Calculation of saturation water vapour pressure
        x = met_vpmax(T) where
        T is the air temperature in degC.
        The saturation water vapour pressure in mbar (hPa) is calculated according to the following formula:
        x = 6.1078 * exp(17.08085 * T / (234.175 + T))
        The parameter may be a variable or a number.

        def met_vpmax(temp, undef=-9999.):


    28. Calculation of actual water vapour pressure
        x = met_vpact(T,rh) where T is the air temperature in degC, rh is the relative humidity in %.
        The actual water vapour pressure in mbar (hPa) is calculated according to the following for- mulas:
        Es = 6.1078*exp(17.08085*T/ (234.175 + T))
        x = Es * rh/100
        Both parameters may be ndarray or float.

        def met_vpact(temp, rh, undef=-9999.):


    29. Calculation of water vapour pressure deficit
        x = met_vpdef(T, rh) where T is the air temperature in degC, rh is the relative humidity in %.
        The water vapour pressure deficit in mbar (hPa) is calculated according to the following for- mulas:
        Es = 6.1078*exp(17.08085*T/ (234.175 + T))
        E = Es * rh/100
        x = Es - E
        Both parameters may be ndarray or float.

        def met_vpdef(temp, rh, undef=-9999.):


    30. Calculation of specific humidity
        x = met_sh(T, rh, p) where
        T is the air temperature in degC,
        rh is the relative humidity in %,
        p is the air pressure in mbar (hPa).
        The specific humidity in g kg-1 is calculated according to the following formulas:
        Es = 6.1078*exp(17.08085*T/ (234.175 + T))
        E = Es * rh/100
        x = 622 * E/(p-0.378*E)
        All parameters may be ndarray or float.

        def met_sh(temp, rh, p, undef=-9999.):


    31. Calculation of potential temperature
        x = met_tpot(T, p) where
        T is the air temperature in degC,
        p is the air pressure in mbar (hPa).
        The potential temperature in K is calculated according to the following formula:
        x = (T + 273.16) * (1000/p)**0.286
        Both parameters may be ndarray or float.

        def met_tpot(temp, p, undef=-9999.):


    32. Calculation of air density
        x = met_rho(T, rh, p) where
        T is the air temperature in degC,
        rh is the relative humidity in %,
        p is the air pressure in mbar (hPa).
        The air density in kg m-3 is calculated according to the following formulas:
        Es = 6.1078*exp(17.08085*T/ (234.175 + T))
        E = Es * rh/100
        sh = 622 * E/(p-0.378*E)
        Tv = ((T + 273.16) * (1 + 0.000608 * sh)) - 273.16
        x = p * 100 / (287.05 * (Tv + 273.16))
        All parameters may be ndarray or float.

        def met_rho(temp, rh, p, undef=-9999.):


    33. Calculation of dew point temperature
        x = met_dpt(T, rh) where
        T is the air temperature in degC, rh is the relative humidity in %.
        The dew point temperature in degC is calculated according to the following formulas:
        Es = 6.1078*exp(17.08085*T/(234.175 + T))
        E = Es * rh/100
        x = 234.175 * ln(E/6.1078)/(17.08085 - ln(E/6.1078))
        Both parameters may be ndarray or float.

        def met_dpt(temp, rh, undef=-9999.):


    34. Calculation of water vapour concentration
        x = met_h2oc(T, rh, p) where T is the air temperature in degC,
        rh is the relative humidity in %,
        p is the air pressure in mbar (hPa).
        The water vapour concentration in mmol mol-1 is calculated according to the following formu- las:
        Es = 6.1078*exp(17.08085*T/ (234.175 + T))
        E = Es * rh/100
        x = 0.1 * E /(0.001*p*100*0.001)
        All parameters may be ndarray or float.

        def met_h2oc(temp, rh, p, undef=-9999.):


    35. Calculation of relative humidity from dry and wet bulb temperature # not implemented


    36. Calculation of relative humidity from dew point temperature # not implemented


    37. Calculation of relative humidity from water vapour concentration
        x = met_h2oc_rh(T, [H2O], p) where
        T is the air temperature in degC,
        [H2O] is the water vapour concentration in mmolmol-1, p is the air pressure in mbar (hPa).
        The relative humidity in % is calculated according to the following formulas:
        Es = 6.1078*exp(17.08085*T/(234.175 + T))
        E = 10 * [H2O] * 0.001 * p * 100 * 0.001
        x = 100 * E / Es
        All parameters may be ndarray or float.

        def met_h2oc_rh(temp, h, p, undef=-9999.):


    38. Rotation of wind direction
        x = met_wdrot(wd, a) where
        wd is the wind direction in degree,
        a is the rotation angle in degree (positive is clockwise).
        The rotated wind direction is calculated according to the following formulas:
        x = wd + a
        IF x > 0 THEN x = x + 360
        IF x >= 360 THEN x = x - 360
        Both parameters may be ndarray or float.

        def met_wdrot(wd, a, undef=-9999.):


    39. Rotation of u-component of wind vector
        x = met_urot(u, v, a) where
        u is the u-component of the wind vector,
        v is the v-component of the wind vector,
        a is the rotation angle in degree (positive is clockwise).
        The rotated u-component is calculated according to the following formula:
        x = u * cos (a) + v * sin (a)
        All three parameters may be ndarray or float.

        def met_urot(u, v, a, undef=-9999.):


    40. Rotation of v-component of wind vector
        x = met_vrot(u, v, a) where
        u is the u-component of the wind vector,
        v is the v-component of the wind vector,
        a is the rotation angle in degree (positive is clockwise).
        The rotated v-component is calculated according to the following formula:
        x = -u * sin (a) + v * cos (a)
        All three parameters may be ndarray or float.

        def met_vrot(u, v, a, undef=-9999.):


    41. Calculation of wind velocity from u- and v-component of wind vector
        x = met_uv_wv(u, v) where
        u is the u-component of the wind vector, v is the v-component of the wind vector.
        The horizontal wind velocity is calculated according to the following formula:
        x = sqrt(u**2 + v**2)
        Both parameters may be ndarray or float.

        def met_uv_wv(u, v, undef=-9999.):


    42. Calculation of wind direction from u- and v-component of wind vector
        x = met_uv_wd(u, v) where
        u is the u-component of the wind vector, v is the v-component of the wind vector.
        The horizontal wind velocity is calculated according to the following formulas:
        IF u = 0 AND v = 0 THEN x = 0
        IF u = 0 AND v > 0 THEN x = 360
        IF u = 0 AND v < 0 THEN x = 180
        IF u < 0 THEN x = 270 - arctan(v/u)
        IF u > 0 THEN x = 90 - arctan(v/u)
        Both parameters may be ndarray or float.

        def met_uv_wd(u, v, undef=-9999.):


    43. Calculation of u-component of wind vector from wind velocity and wind direction
        x = met_wvwd_u(wv, wd) where wv is the horizontal wind velocity, wd is the horizontal wind direction.
        The u-component of the wind vector is calculated according to the following formula:
        x = -wv * sin (wd)
        Both parameters may be ndarray or float.

        def met_wvwd_u(wv, wd, undef=-9999.):


    44. Calculation of v-component of wind vector from wind velocity and wind direction
        x = met_wvwd_v(wv, wd) where wv is the horizontal wind velocity, wd is the horizontal wind direction.
        The v-component of the wind vector is calculated according to the following formula:
        x = -wv * cos (wd)
        Both parameters may be ndarray or float.

        def met_wvwd_v(wv, wd, undef=-9999.):


    45. If-statements
        x = ifeq(y,a0,a1,a2) means IF y == a0 THEN x = a1 ELSE x = a2
        x = ifne(y,a0,a1,a2) means IF y != a0 THEN x = a1 ELSE x = a2
        x = ifle(y,a0,a1,a2) means IF y <= a0 THEN x = a1 ELSE x = a2
        x = ifge(y,a0,a1,a2) means IF y >= a0 THEN x = a1 ELSE x = a2
        x = iflt(y,a0,a1,a2) means IF y > a0 THEN x = a1 ELSE x = a2
        x = ifgt(y,a0,a1,a2) means IF y < a0 THEN x = a1 ELSE x = a2
        All parameters may be ndarray or float.

        def ifeq(var1, iif, ithen, ielse):
        def ifne(var1, iif, ithen, ielse):
        def ifle(var1, iif, ithen, ielse):
        def ifge(var1, iif, ithen, ielse):
        def iflt(var1, iif, ithen, ielse):
        def ifgt(var1, iif, ithen, ielse):


    46. Write variables to a file # not implemented


This module was written by Matthias Cuntz while at Department of
Computational Hydrosystems, Helmholtz Centre for Environmental
Research - UFZ, Leipzig, Germany, and continued while at Institut
National de Recherche pour l'Agriculture, l'Alimentation et
l'Environnement (INRAE), Nancy, France.

Copyright (c) 2014-2020 Matthias Cuntz - mc (at) macu (dot) de
Released under the MIT License; see LICENSE file for details.

* Written Jun-Dec 2014 by Matthias Cuntz (mc (at) macu (dot) de)
* Corrected type in met_tpot, Jun 2014, Corinna Rebmann
* Changed to Sphinx docstring and numpydoc, May 2020, Matthias Cuntz

.. moduleauthor:: Matthias Cuntz

.. autosummary::
   varchs
   varadd
   varsub
   varmul
   vardiv
   varsqr
   varsqrt
   varexp
   varlog
   varpot
   lin
   quad
   cubic
   hms
   bit_test
   setlow
   sethigh
   limits
   mean
   mini
   maxi
   met_lwrad
   met_trad
   met_alb
   met_albl
   met_vpmax
   met_vpact
   met_vpdef
   met_sh
   met_tpot
   met_rho
   met_dpt
   met_h2oc
   met_h2oc_rh
   met_wdrot
   met_urot
   met_vrot
   met_uv_wv
   met_uv_wd
   met_wvwd_u
   met_wvwd_v
   ifeq
   ifne
   ifle
   ifge
   iflt
   ifgt
"""
from __future__ import division, absolute_import, print_function
import numpy as np
try:        # import package
    from ..division import division
    from ..esat     import esat
    from ..const    import sigma, T0
except:
    try:    # e.g. python module.py at main package level
        from division import division
        from esat     import esat
        from const    import sigma, T0
    except: # python logtools.py
        division = _div
        esat     = _esat
        sigma = 5.67e-08 # Stefan-Boltzmann constant [W m^-2 K^-4]
        T0    = 273.15   # Celcius <-> Kelvin [K]


__all__ = ['varchs', 'varadd', 'varsub', 'varmul', 'vardiv', 'varsqr',
           'varsqrt', 'varexp', 'varlog', 'varpot', 'lin', 'quad',
           'cubic', 'hms', 'bit_test', 'setlow', 'sethigh', 'limits',
           'mean', 'mini', 'maxi', 'met_lwrad', 'met_trad', 'met_alb',
           'met_albl', 'met_vpmax', 'met_vpact', 'met_vpdef', 'met_sh',
           'met_tpot', 'met_rho', 'met_dpt', 'met_h2oc', 'met_h2oc_rh',
           'met_wdrot', 'met_urot', 'met_vrot', 'met_uv_wv', 'met_uv_wd',
           'met_wvwd_u', 'met_wvwd_v', 'ifeq', 'ifne', 'ifle', 'ifge',
           'iflt', 'ifgt']


# Not implemented: varset


def varchs(a, undef=-9999.):
    """
    Change sign:
        x = varchs(a) means x = -a, where a is ndarray or float.

    Parameters
    ----------
    a : ndarray
        input variable
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        Changed sign

    History
    -------
    Written,  Matthias Cuntz, Dec 2014
    """
    return np.where(a==undef, undef, -a)


def varadd(a, b, undef=-9999.):
    """
    Addition:
        x = varadd(a, b) means x = a + b, where a and b are ndarray or float.

    Parameters
    ----------
    a : ndarray
        input variable 1
    b : ndarray
        input variable 2
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        Addition

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    return np.where((a==undef) | (b==undef), undef, a + b)


def varsub(a, b, undef=-9999.):
    """
    Subtraction:
        x = varsub(a, b) means x = a - b, where a and b are ndarray or float.

    Parameters
    ----------
    a : ndarray
        input variable 1
    b : ndarray
        input variable 2
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        Subtraction

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    return np.where((a==undef) | (b==undef), undef, a - b)


def varmul(a, b, undef=-9999.):
    """
    Multiplication:
        x = varmul(a, b) means x = a * b, where a and b are ndarray or float.

    Parameters
    ----------
    a : ndarray
        input variable 1
    b : ndarray
        input variable 2
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        Multiplication

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    return np.where((a==undef) | (b==undef), undef, a * b)


def vardiv(a, b, undef=-9999.):
    """
    Division:
        x = vardiv(a, b) means x = a/b, where a and b are ndarray or float.

    Parameters
    ----------
    a : ndarray
        dividend
    b : ndarray
        divisor
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        Division

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    return np.where((a==undef) | (b==undef), undef, division(a, b, undef))


def varsqr(a, undef=-9999.):
    """
    Square root:
        x = varsqr(a) means x = sqrt(a), where a is ndarray or float.

    Parameters
    ----------
    a : ndarray
        input variable
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        Square root

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    return np.where((a==undef), undef, np.sqrt(a))


def varsqrt(a, undef=-9999.):
    """
    Square root:
        x = varsqrt(a) means x = sqrt(a), where a is ndarray or float.

    Parameters
    ----------
    a : ndarray
        input variable
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        Square root

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    return np.where((a==undef), undef, np.sqrt(a))


def varexp(a, undef=-9999.):
    """
    Exponentiation of e:
        x = varexp(a) means x = exp(a), where a is ndarray or float.

    Parameters
    ----------
    a : ndarray
        exponent
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        Exponentiation

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    return np.where((a==undef), undef, np.exp(a))


def varlog(a, undef=-9999.):
    """
    Natural logarithm:
        x = varlog(a) means x = ln(a), where a is ndarray or float.

    Parameters
    ----------
    a : ndarray
        input variable
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        Natural logarithm

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    return np.where((a==undef), undef, np.log(a))


def varpot(a, b, undef=-9999.):
    """
    Exponentiation:
        x = varpot(a, b) means x = a**b, where a and b are ndarray or float.

    Parameters
    ----------
    a : ndarray
        base
    b : ndarray
        exponent
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        Exponentiation

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    return np.where((a==undef) | (b==undef), undef, a**b)


def lin(y, a0, a1, undef=-9999.):
    """
    Apply linear function:
        x = lin(y, a0, a1) means x = a0 + a1 * y

    Parameters
    ----------
    y : ndarray
        input variable
    a0 : ndarray or float
        parameter 1
    a1 : ndarray or float
        parameter 2
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        linear function

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    return np.where((y==undef), undef, a0 + a1*y)


def quad(y, a0, a1, a2, undef=-9999.):
    """
    Apply 2nd order function:
        x=quad(y,a0,a1,a2) means x = a0 + a1*y + a2*y**2

    Parameters
    ----------
    y : ndarray
        input variable
    a0 : ndarray or float
        parameter 1
    a1 : ndarray or float
        parameter 1
    a2 : ndarray or float
        parameter 1
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        2nd order function

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    return np.where((y==undef), undef, a0 + a1*y + a2*y*y)


def cubic(y, a0, a1, a2, a3, undef=-9999.):
    """
    Apply 3rd order function:
        x=cubic(y,a0,a1,a2,a3) means x = a0 + a1*y + a2*y**2 + a3*y**3

    Parameters
    ----------
    y : ndarray
        input variable
    a0 : ndarray or float
        parameter 1
    a1 : ndarray or float
        parameter 2
    a2 : ndarray or float
        parameter 3
    a3 : ndarray or float
        parameter 4
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        3rd order function

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    return np.where((y==undef), undef, a0 + a1*y + a2*y*y + a3*y*y*y)


def hms(h, m, s, undef=-9999.):
    """
    Calculate fraction of day from hours, minutes and seconds:
        x = hms(h, m, s) means x = (h + m/60 + s/3600)/24

    Parameters
    ----------
    h : ndarray
        hour
    m : ndarray
        minute
    s : ndarray
        second
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        fraction of day

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    return np.where((h==undef) | (m==undef) | (s==undef), undef, (h+m/60.+s/3600.)/24.)


def bit_test(y, b, start=0):
    """
    Bitwise test:
        x = bit_test(y, b, start=0) means x = 1 if bit b is set in y otherwise x = 0.

    Returns a list if b is an array.

    Counting of b starts at start.

    For the behaviour of the original logger tools, set start=1.
    Negative b's is not implemented.

    Parameters
    ----------
    y : ndarray
        input variable 1
    b : int or ndarray
        input variable 2
    start : int, optional
        Counting of `b` starts at start (default: 0)

    Returns
    -------
    int or list
        Bitwise test

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    if np.size(b) > 1:
        return [ (y >> i+start)%2 for i in b ]
    else:
        return (y >> b+start)%2


def setlow(y, low, islow=None, undef=-9999.):
    """
    Replacement of underflows by new value:
        x = setlow(y,low,islow) means IF (y < low) THEN x = islow ELSE x = y

    islow is optional. If not given low will be used.

    This function may be used to adjust small negative values of short wave radiation
    during nighttime to zero values.

    Parameters
    ----------
    y : ndarray
        input variable
    low : ndarray
        lower threshold
    islow : None or ndarray, optional
        if not None, use islow in case of y < low
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        underflows replaced by new value

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    if islow is None:
        out = np.maximum(y, low)
    else:
        out = np.where(y < low, islow, y)
    return np.where(y == undef, undef, out)


def sethigh(y, high, ishigh=None, undef=-9999.):
    """
    Replacement of overflows by new value:
        x = sethigh(y,high,ishigh) means IF (y > high) THEN x = ishigh ELSE x = y

    ishigh is optional. If not given high will be used.

    This function may be used to adjust relative humidity values of a little bit more than 100 % to 100 %.

    Parameters
    ----------
    y : ndarray
        input variable
    high : ndarray
        upper threshold
    ishigh : None or ndarray, optional
        if not None, use ishigh in case of y > high
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        overflows replaced by new value

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    if ishigh is None:
        out = np.minimum(y, high)
    else:
        out = np.where(y > high, ishigh, y)
    return np.where(y == undef, undef, out)


def limits(y, mini, maxi, undef=-9999.):
    """
    Replacement of underflows or overflows by undef:
        x = limits(y, mini, maxi) means IF (y < mini) OR (y > maxi) THEN x = undef ELSE x = y

    This function may be used to check values lying in between certain limits.

    If one of the limits is exceeded the value is set to undef.

    Parameters
    ----------
    y : ndarray
        input variable
    mini : ndarray
        lower threshold
    maxi : ndarray
        upper threshold
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        underflows or overflows replaced by `undef`

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    return np.where((y >= mini) & (y <= maxi), y, undef)


def mean(y, axis=None, undef=-9999.):
    """
    Calculation of mean value:
        x = mean(y) means x = (y[0] + y[1] + ... + y[n-1])/n

    Parameters
    ----------
    y : ndarray
        input variable
    axis : None or int or tuple of ints, optional
        Axis or axes along which the means are computed.
        The default is to compute the mean of the flattened array.

        If this is a tuple of ints, a mean is performed over multiple axes,
        instead of a single axis or all the axes as before.
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        mean value

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    return np.ma.mean(np.ma.array(y, mask=(y==undef)), axis=axis).filled(undef)


def mini(y, axis=None, undef=-9999.):
    """
    Calculation of minimum value:
        x = mini(y) means x = min(y[0],y[1],...,y[n-1])

    Parameters
    ----------
    y : ndarray
        input variable
    axis : None or int or tuple of ints, optional
        Axis or axes along which the minimum are computed.
        The default is to compute the minimum of the flattened array.

        If this is a tuple of ints, a minimum is performed over multiple axes,
        instead of a single axis or all the axes as before.
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        minimum value

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    return np.ma.amin(np.ma.array(y, mask=(y==undef)), axis=axis).filled(undef)


def maxi(y, axis=None, undef=-9999.):
    """
    Calculation of maximum value:
        x = maxi(y) means x = max(y[0],y[1],...,y[n-1])

    Parameters
    ----------
    y : ndarray
        input variable
    axis : None or int or tuple of ints, optional
        Axis or axes along which the maximum are computed.
        The default is to compute the maximum of the flattened array.

        If this is a tuple of ints, a maximum is performed over multiple axes,
        instead of a single axis or all the axes as before.
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        maximum value

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    return np.ma.amax(np.ma.array(y, mask=(y==undef)), axis=axis).filled(undef)


# Not implemented: met_torad


def met_lwrad(y, Tp, undef=-9999.): # assumes that y was already multiplied with calibration factor
    """
    Calculation of long wave radiation from net radiometer:
        x = met_lwrad(y, Tp)

    The total radiation in W m-2 is calculated according to the following formula:
        x = y*fl + sigma*(Tp+T0)**4

    where sigma = 5.67051 * 10**8 W m-2 K-4 is the Stephan-Boltzmann-Constant and
    fl is the factor for long wave radiation (reciprocal value of sensitivity) in W m-2 per mV.

    The function assumes that fl was already applied before.

    Parameters
    ----------
    y : ndarray
        output voltage of the net radiometer [mV]
    Tp : ndarray
        pyranometer temperature, i.e. the temperature of the net radiometer body [degC]
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        total radiation in W m-2

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    return np.where((y==undef) | (Tp==undef), undef, y + sigma * (Tp+T0)**4)


def met_trad(Rl, epsilon, undef=-9999.):
    """
    Calculation of radiation temperature from long wave radiation:
        x = met_trad(Rl, epsilon)

    The radiation temperature in degC is calculated according to the following formula:
        x= sqrt4(Rl/(sigma*epsilon)) - T0

    where sigma = 5.67051 * 108 W m-2 K-4 is the Stephan-Boltzmann-Constant.

    Parameters
    ----------
    Rl : ndarray
        longwave radiation [W m-2]
    epsilon : ndarray
        long wave emissivity of the surface [0-1]
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        radiation temperature in degC

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    const = 1./(epsilon*sigma)
    trad  = np.ma.sqrt(np.ma.sqrt(const*np.ma.array(Rl, mask=(Rl==undef)))) - T0
    return trad.filled(undef)


def met_alb(swd, swu, swdmin=50., swumin=10., undef=-9999.):
    """
    Calculation of albedo from short wave downward and upward radiation:
        x = met_alb(swd, swu)

    The albedo in % is calculated according to the following formula:
        x = 100 * ( swu / swd )

    If swd < swdmin (50 W m-2) or swu < swumin (10 W m-2) the result is undef.

    Parameters
    ----------
    swd : ndarray
        shortwave downward radiation [W m-2]
    swu : ndarray
        shortwave upward radiation [W m-2]
    swdmin : float, optional
        If `swd` < `swdmin` the result is undef (default: 50).
    swumin : float, optional
        If `swu` < `swumin` the result is undef (default: 10).
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        albedo in %

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    return np.where((swd==undef) | (swu==undef) | (swd<swdmin) | (swu<swumin),
                    undef, division(swu*100., swd, undef))


def met_albl(swd, swu, swdmin, swumin, undef=-9999.):
    """
    Calculation of albedo from short wave downward and upward radiation with limits:
        x=met_albl(swd,swu,swdmin,swumin)

    The albedo in % is calculated according to the following formula:
        x = 100 * ( swu / swd )

    If swd < swdmin or swu < swumin the result is `undef`.

    Parameters
    ----------
    swd : ndarray
        shortwave downward radiation [W m-2]
    swu : ndarray
        shortwave upward radiation [W m-2]
    swdmin : float
        If `swd` < `swdmin` the result is undef.
    swumin : float
        If `swu` < `swumin` the result is undef.
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        albedo in %

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    return np.where((swd==undef) | (swu==undef) | (swd<swdmin) | (swu<swumin),
                    undef, division(swu*100., swd, undef))


def met_vpmax(temp, undef=-9999.):
    """
    Calculation of saturation water vapour pressure:
        x = met_vpmax(T) where

    The saturation water vapour pressure in mbar (hPa) is calculated according to the following formula:
        x = 6.1078 * exp(17.08085 * T / (234.175 + T))

    Parameters
    ----------
    temp : ndarray
        air temperature [degC]
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        saturation water vapour pressure in mbar (hPa)

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    es = esat(np.ma.array(temp+T0, mask=(temp==undef)))*0.01
    return es.filled(undef)


def met_vpact(temp, rh, undef=-9999.):
    """
    Calculation of actual water vapour pressure:
        x = met_vpact(T,rh)

    The actual water vapour pressure in mbar (hPa) is calculated according to the following formulas:
        Es = 6.1078*exp(17.08085*T/ (234.175 + T))

        x = Es * rh/100

    Parameters
    ----------
    temp : ndarray
        air temperature [degC]
    rh : ndarray
        relative humidity [%]
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        actual water vapour pressure in mbar (hPa)

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    es = esat(np.ma.array(temp+T0, mask=((temp==undef)|(rh==undef))))*0.01
    ea = es*rh*0.01
    return ea.filled(undef)


def met_vpdef(temp, rh, undef=-9999.):
    """
    Calculation of water vapour pressure deficit:
        x = met_vpdef(T, rh)

    The water vapour pressure deficit in mbar (hPa) is calculated according to the following formulas:
        Es = 6.1078*exp(17.08085*T/ (234.175 + T))

        E = Es * rh/100

        x = Es - E

    Parameters
    ----------
    temp : ndarray
        air temperature [degC]
    rh : ndarray
        relative humidity [%]
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        water vapour pressure deficit in mbar (hPa)

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    es  = esat(np.ma.array(temp+T0, mask=((temp==undef)|(rh==undef))))*0.01
    ea  = es*rh*0.01
    vpd = es - ea
    return vpd.filled(undef)


def met_sh(temp, rh, p, undef=-9999.):
    """
    Calculation of specific humidity:
        x = met_sh(T, rh, p)

    The specific humidity in g kg-1 is calculated according to the following formulas:
        Es = 6.1078*exp(17.08085*T/ (234.175 + T))

        E = Es * rh/100

        x = 622 * E/(p-0.378*E)

    Parameters
    ----------
    temp : ndarray
        air temperature [degC]
    rh : ndarray
        relative humidity [%]
    p : ndarray
        air pressure [hPa = mbar]
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        specific humidity in g kg-1

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    es = esat(np.ma.array(temp+T0, mask=((temp==undef)|(rh==undef)|(p==undef))))*0.01
    ea = es*rh*0.01
    sh =  division(622.*ea, (p-0.378*ea), undef)
    return sh.filled(undef)


def met_tpot(temp, p, undef=-9999.):
    """
    Calculation of potential temperature:
        x = met_tpot(T, p)

    The potential temperature in K is calculated according to the following formula:
        x = (T + T0) * (1000/p)**0.286

    Parameters
    ----------
    temp : ndarray
        air temperature [degC]
    p : ndarray
        air pressure [hPa = mbar]
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        potential temperature in K

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    return np.where((temp==undef) | (p==undef), undef, (temp+T0)*division(1000.,p)**0.286)


def met_rho(temp, rh, p, undef=-9999.):
    """
    Calculation of air density:
        x = met_rho(T, rh, p)

    The air density in kg m-3 is calculated according to the following formulas:
        Es = 6.1078*exp(17.08085*T/ (234.175 + T))

        E = Es * rh/100

        sh = 622 * E/(p-0.378*E)

        Tv = ((T + T0) * (1 + 0.000608 * sh)) - T0

        x = p * 100 / (287.05 * (Tv + T0))

    Parameters
    ----------
    temp : ndarray
        air temperature [degC]
    rh : ndarray
        relative humidity [%]
    p : ndarray
        air pressure [hPa = mbar]
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        air density in kg m-3

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    es  = esat(np.ma.array(temp+T0, mask=((temp==undef)|(rh==undef)|(p==undef))))*0.01
    ea  = es*rh*0.01
    sh  = division(622.*ea, (p-0.378*ea), undef)
    Tv  = ((temp+T0)*(1+0.000608*sh)) - T0
    rho = division(p*100., (287.05*(Tv+T0)), undef)
    return rho.filled(undef)


def met_dpt(temp, rh, undef=-9999.):
    """
    Calculation of dew point temperature:
        x = met_dpt(T, rh)

    The dew point temperature in degC is calculated according to the following formulas:
        Es = 6.1078*exp(17.08085*T/(234.175 + T))

        E = Es * rh/100

        x = 234.175 * ln(E/6.1078)/(17.08085 - ln(E/6.1078))

    Parameters
    ----------
    temp : ndarray
        air temperature [degC]
    rh : ndarray
        relative humidity [%]
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        dew point temperature in degC

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    es = esat(np.ma.array(temp+T0, mask=((temp==undef)|(rh==undef))))*0.01
    ea = es*rh*0.01
    dpt = 234.175 * np.ma.log(ea/6.1078) / (17.08085 - np.ma.log(ea/6.1078))
    return dpt.filled(undef)


def met_h2oc(temp, rh, p, undef=-9999.):
    """
    Calculation of water vapour concentration:
        x = met_h2oc(T, rh, p)

    The water vapour concentration in mmol mol-1 is calculated according to the following formulas:
        Es = 6.1078*exp(17.08085*T/ (234.175 + T))

        E = Es * rh/100

        x = 0.1 * E /(0.001*p*100*0.001)

    Parameters
    ----------
    temp : ndarray
        air temperature [degC]
    rh : ndarray
        relative humidity [%]
    p : ndarray
        air pressure [hPa = mbar]
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        water vapour concentration in mmol mol-1

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    es = esat(np.ma.array(temp+T0, mask=((temp==undef)|(rh==undef)|(p==undef))))*0.01
    ea = es*rh*0.01
    c  = division(1000.*ea, p, undef)
    return c.filled(undef)


# Not implemented: met_psy_rh


# Not implemented: met_dpt_rh


def met_h2oc_rh(temp, h, p, undef=-9999.):
    """
    Calculation of relative humidity from water vapour concentration:
        x = met_h2oc_rh(T, [H2O], p)

    The relative humidity in % is calculated according to the following formulas:
        Es = 6.1078*exp(17.08085*T/(234.175 + T))

        E = 10 * [H2O] * 0.001 * p * 100 * 0.001

        x = 100 * E / Es

    Parameters
    ----------
    temp : ndarray
        air temperature [degC]
    h : ndarray
        water vapour concentration [mmol mol-1]
    p : ndarray
        air pressure [hPa = mbar]
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        relative humidity in %

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    es = esat(np.ma.array(temp+T0, mask=((temp==undef)|(h==undef)|(p==undef))))*0.01
    ea = 0.001 * h * p
    c  = 100.*ea/es
    return c.filled(undef)


def met_wdrot(wd, a, undef=-9999.):
    """
    Rotation of wind direction:
        x = met_wdrot(wd, a)

    The rotated wind direction is calculated according to the following formulas:
        x = wd + a

        IF x < 0 THEN x = x + 360

        IF x <= 360 THEN x = x - 360

    Parameters
    ----------
    wd : ndarray
        wind direction [degree]
    a : ndarray
        rotation angle (positive is clockwise) [degree]
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        rotated wind direction

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    rot = np.ma.array(wd+a, mask=((wd==undef)|(a==undef)))
    rot = np.ma.where(rot < 0., rot+360., rot)
    rot = np.ma.where(rot >= 360., rot-360., rot)
    return rot.filled(undef)


def met_urot(u, v, a, undef=-9999.):
    """
    Rotation of u-component of wind vector:
        x = met_urot(u, v, a)

    The rotated u-component is calculated according to the following formula:
        x = u * cos (a) + v * sin (a)

    Parameters
    ----------
    u : ndarray
        u-component of the wind vector
    v : ndarray
        v-component of the wind vector
    a : ndarray
        rotation angle (positive is clockwise) [degree]
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        rotated u-component

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    return np.where((u==undef) | (v==undef) | (a==undef), undef, u*np.cos(np.deg2rad(a)) + v*np.sin(np.deg2rad(a)))


def met_vrot(u, v, a, undef=-9999.):
    """
    Rotation of v-component of wind vector:
        x = met_vrot(u, v, a)

    The rotated v-component is calculated according to the following formula:
        x = -u * sin (a) + v * cos (a)

    Parameters
    ----------
    u : ndarray
        u-component of the wind vector
    v : ndarray
        v-component of the wind vector
    a : ndarray
        rotation angle (positive is clockwise) [degree]
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        rotated v-component

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    return np.where((u==undef) | (v==undef) | (a==undef), undef, -u*np.sin(np.deg2rad(a)) + v*np.cos(np.deg2rad(a)))


def met_uv_wv(u, v, undef=-9999.):
    """
    Calculation of wind velocity from u- and v-component of wind vector:
        x = met_uv_wv(u, v)

    The horizontal wind velocity is calculated according to the following formula:
        x = sqrt(u**2 + v**2)

    Parameters
    ----------
    u : ndarray
        u-component of the wind vector
    v : ndarray
        v-component of the wind vector
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        horizontal wind velocity

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    return np.where((u==undef) | (v==undef), undef, np.sqrt(u*u + v*v))


def met_uv_wd(u, v, undef=-9999.):
    """
    Calculation of wind direction from u- and v-component of wind vector:
        x = met_uv_wd(u, v)

    The horizontal wind velocity is calculated according to the following formulas:
        IF u = 0 AND v = 0 THEN x = 0

        IF u = 0 AND v < 0 THEN x = 360

        IF u = 0 AND v > 0 THEN x = 180

        IF u > 0 THEN x = 270 - arctan(v/u)

        IF u < 0 THEN x = 90 - arctan(v/u)

    Parameters
    ----------
    u : ndarray
        u-component of the wind vector
    v : ndarray
        v-component of the wind vector
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        horizontal wind velocity

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    wd = np.ma.zeros(u.shape)
    wd.mask = (u==undef) | (v==undef)
    wd = np.ma.where((u==0.) & (v==0.), 0., wd)
    wd = np.ma.where((u==0.) & (v<0.), 360., wd)
    wd = np.ma.where((u==0.) & (v>0.), 180., wd)
    wd = np.ma.where((u>0.), 270.-np.rad2deg(np.ma.arctan(v/u)), wd)
    wd = np.ma.where((u<0.), 90.-np.rad2deg(np.ma.arctan(v/u)), wd)
    return wd.filled(undef)


def met_wvwd_u(wv, wd, undef=-9999.):
    """
    Calculation of u-component of wind vector from wind velocity and wind direction:
        x = met_wvwd_u(wv, wd)

    The u-component of the wind vector is calculated according to the following formula:
        x = -wv * sin (wd)

    Parameters
    ----------
    wv : ndarray
        horizontal wind velocity [m s-1]
    wd : ndarray
        horizontal wind direction [degree]
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        u-component of the wind vector

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    return np.where((wv==undef) | (wd==undef), undef, -wv*np.sin(np.deg2rad(wd)))


def met_wvwd_v(wv, wd, undef=-9999.):
    """
    Calculation of v-component of wind vector from wind velocity and wind direction:
        x = met_wvwd_v(wv, wd)

    The v-component of the wind vector is calculated according to the following formula:
        x = -wv * cos (wd)

    Parameters
    ----------
    wv : ndarray
        horizontal wind velocity [m s-1]
    wd : ndarray
        horizontal wind direction [degree]
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        v-component of the wind vector

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    return np.where((wv==undef) | (wd==undef), undef, -wv*np.cos(np.deg2rad(wd)))


def ifeq(y, a0, a1, a2):
    """
    If-statements:
        x = ifeq(y,a0,a1,a2) means IF y == a0 THEN x = a1 ELSE x = a2

    Parameters
    ----------
    y : ndarray
        input variable
    a0 : ndarray
        compare to input `y == a0`
    a1 : ndarray
        result if `y == a0`
    a2 : ndarray
        result if `y != a0`
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        IF y == a0 THEN x = a1 ELSE x = a2

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    out = np.where(y == a0, a1, a2)
    return np.where((y==undef) | (a0==undef) | (a1==undef) | (a2==undef), undef, out)


def ifne(y, a0, a1, a2):
    """
    If-statements:
        x = ifne(y,a0,a1,a2) means IF y != a0 THEN x = a1 ELSE x = a2

    Parameters
    ----------
    y : ndarray
        input variable
    a0 : ndarray
        compare to input `y != a0`
    a1 : ndarray
        result if `y != a0`
    a2 : ndarray
        result if `y == a0`
    y : ndarray
    a0 : ndarray
    a1 : ndarray
    a2 : ndarray
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        IF y != a0 THEN x = a1 ELSE x = a2

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    out = np.where(y != a0, a1, a2)
    return np.where((y==undef) | (a0==undef) | (a1==undef) | (a2==undef), undef, out)


def ifle(y, a0, a1, a2):
    """
    If-statements:
        x = ifle(y,a0,a1,a2) means IF y >= a0 THEN x = a1 ELSE x = a2

    Parameters
    ----------
    y : ndarray
        input variable
    a0 : ndarray
        compare to input `y <= a0`
    a1 : ndarray
        result if `y <= a0`
    a2 : ndarray
        result if `y > a0`
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        IF y >= a0 THEN x = a1 ELSE x = a2

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    out = np.where(y <= a0, a1, a2)
    return np.where((y==undef) | (a0==undef) | (a1==undef) | (a2==undef), undef, out)


def ifge(y, a0, a1, a2):
    """
    If-statements:
        x = ifge(y,a0,a1,a2) means IF y <= a0 THEN x = a1 ELSE x = a2

    Parameters
    ----------
    y : ndarray
        input variable
    a0 : ndarray
        compare to input `y >= a0`
    a1 : ndarray
        result if `y >= a0`
    a2 : ndarray
        result if `y < a0`
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        IF y <= a0 THEN x = a1 ELSE x = a2

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    out = np.where(y <= a0, a1, a2)
    return np.where((y==undef) | (a0==undef) | (a1==undef) | (a2==undef), undef, out)


def iflt(y, a0, a1, a2):
    """
    If-statements:
        x = iflt(y,a0,a1,a2) means IF y < a0 THEN x = a1 ELSE x = a2

    Parameters
    ----------
    y : ndarray
        input variable
    a0 : ndarray
        compare to input `y < a0`
    a1 : ndarray
        result if `y < a0`
    a2 : ndarray
        result if `y >= a0`
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        IF y < a0 THEN x = a1 ELSE x = a2

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    out = np.where(y < a0, a1, a2)
    return np.where((y==undef) | (a0==undef) | (a1==undef) | (a2==undef), undef, out)


def ifgt(y, a0, a1, a2):
    """
    If-statements:
        x = ifgt(y,a0,a1,a2) means IF y > a0 THEN x = a1 ELSE x = a2

    Parameters
    ----------
    y : ndarray
        input variable
    a0 : ndarray
        compare to input `y > a0`
    a1 : ndarray
        result if `y > a0`
    a2 : ndarray
        result if `y <= a0`
    undef : float, optional
        elements are excluded from the calculations if any of the inputs equals `undef` (default: -9999.)

    Returns
    -------
    ndarray
        IF y > a0 THEN x = a1 ELSE x = a2

    History
    -------
    Written,  Matthias Cuntz, Jun 2014
    """
    out = np.where(y > a0, a1, a2)
    return np.where((y==undef) | (a0==undef) | (a1==undef) | (a2==undef), undef, out)


# Not implemented: write


#
# Local replacement functions if helper functions do not exist in library
#

def _div(a, b, otherwise=np.nan, prec=0.):
    """
    Divide two arrays, return `otherwise` if division by 0.

    Copy of ..division.py

    Parameters
    ----------
    a : array_like
        enumerator
    b : array_like
        denominator
    otherwise : float
        value to return if `b=0` (default: `np.nan`)
    prec : float
        if |b|<|prec| then `otherwise`

    Returns
    -------
    ndarray
        ratio : numpy array or masked array
        a/b        if |b| >  |prec|

        otherwise  if |b| <= |prec|

        Output is numpy array. It is a masked array if at least one
        of `a` or `b` is a masked array.
    """
    oldsettings = np.geterr()
    np.seterr(divide='ignore')

    if isinstance(a, np.ma.masked_array) or isinstance(b, np.ma.masked_array):
        out = np.ma.where(np.ma.abs(np.ma.array(b)) > np.abs(prec), np.ma.array(a)/np.ma.array(b), otherwise)
    else:
        out = np.where(np.abs(np.array(b)) > np.abs(prec), np.array(a)/np.array(b), otherwise)

    np.seterr(**oldsettings)

    return out


def _esat(T):
    """
    Calculates the saturation vapour pressure of water with the Goff-Gratch formula.

    From ..esat.py

    Parameters
    ----------
    T : ndarray
        Temperature [K]

    Returns
    -------
    ndarray
        Saturation water pressure at temperature T in Pascal [Pa].
    """
    Ts  = 373.16    # steam point temperature in K
    ews = 1013.246  # saturation pressure at steam point temperature, normal atmosphere
    esat_liq = 10.**(-7.90298*(Ts/T-1.) + 5.02808 * np.ma.log10(Ts/T)
                     - 1.3816e-7 * (10.**(11.344*(1.-T/Ts))-1.)
                     + 8.1328e-3 * (10.**(-3.49149*(Ts/T-1.))-1.) + np.ma.log10(ews)) # [hPa]
    return esat_liq * 100. # [Pa]


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
