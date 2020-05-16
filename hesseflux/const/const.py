#!/usr/bin/env python
"""
const : Provides physical, mathematical, computational, isotope, and material constants.

Defines the following constants:
    Mathematical
        Pi, Pi2, Pi3, TwoPi, Sqrt2, pi, pi2, pi3, Twopi

    Physical
        Gravity, T0, P0, T25, sigma, R, R_air, R_H2O, Na, REarth

    Isotope
        R13VPDB, R18VSMOW, R2VSMOW

    Computational
        tiny, huge, eps

    Material
        mmol_co2, mmol_h2o, mmol_air,
        density_quartz, cheat_quartz, cheat_water, cheat_air, latentheat_vaporization

This module was written by Matthias Cuntz while at Department of
Computational Hydrosystems, Helmholtz Centre for Environmental
Research - UFZ, Leipzig, Germany, and continued while at Institut
National de Recherche pour l'Agriculture, l'Alimentation et
l'Environnement (INRAE), Nancy, France.

Examples
--------
>>> print({:.3f}.format(Pi))
3.142

>>> print({:.3f}.format(Sqrt2))
1.414

>>> print({:.3f}.format(Gravity))
9.810

>>> print({:.3f}.format(T0))
273.150

>>> print({:.3f}.format(sigma))
5.670e-08

>>> print({:.3f}.format(R13VPDB))
0.011

>>> print({:.3f}.format(tiny))
1.000e-06

>>> print({:.3f}.format(REarth))
6371000.000

>>> print({:.3f}.format(mmol_h2o))
18.015

>>> print({:.3f}.format(mmol_air))
28.964

>>> print({:.3f}.format(density_quartz))
2.650

>>> print({:.3f}.format(cheat_quartz))
800.000

>>> print({:.3f}.format(cheat_water))
4180.000

>>> print({:.3f}.format(cheat_air))
1010.000

>>> print({:.3f}.format(latentheat_vaporization))
2.450e+06

Copyright (c) 2012-2020 Matthias Cuntz - mc (at) macu (dot) de
Released under the MIT License; see LICENSE file for details.

* Written Jan 2012 by Matthias Cuntz (mc (at) macu (dot) de)
* Ported to Python 3, Feb 2013, Matthias Cuntz
* Added dielectric constant for water, Mar 2014, Arndt Piayda
* Added heat capacities for air, water and quartz as well as density of quartz, Sep 2014, Arndt Piayda
* Added Pi3=pi/3, R13VPDB, R18VSMOW, R2VSMOW, Mar 2015, Matthias Cuntz
* Renamed heat capacities, molar masses, density of quartz, Mar 2015, Matthias Cuntz
* Moved calculation of dielectric constant of water to own routine, Mar 2015, Matthias Cuntz
* Added computational constants such as tiny=np.finfo(np.float).tiny, Nov 2016, Matthias Cuntz
* Added gas constants for dry air and water, May 2017, RL
* Using numpy docstring format, May 2020, Matthias Cuntz
* Added lowercase version of pi constants, May 2020, Matthias Cuntz

.. moduleauthor:: Matthias Cuntz

.. autosummary::
   Pi
   Pi2
   Pi3
   TwoPi
   pi
   pi2
   pi3
   Twopi
   Sqrt2
   Gravity
   T0
   P0
   T25
   sigma
   R
   R_air
   R_H2O
   Na
   REarth
   mmol_co2
   mmol_h2o
   mmol_air
   density_quartz
   cheat_quartz
   cheat_water
   cheat_air
   latentheat_vaporization
   R13VPDB
   R18VSMOW
   R2VSMOW
   tiny
   huge
   eps
"""
from __future__ import division, absolute_import, print_function
import numpy as np


__all__ = ['Pi', 'Pi2', 'Pi3', 'TwoPi', 'pi', 'pi2', 'pi3', 'Twopi', 'Sqrt2',
           'Gravity', 'T0', 'P0', 'T25', 'sigma', 'R', 'R_air', 'R_H2O', 'Na', 'REarth',
           'mmol_co2', 'mmol_h2o', 'mmol_air',
           'density_quartz', 'cheat_quartz', 'cheat_water', 'cheat_air', 'latentheat_vaporization',
           'R13VPDB', 'R18VSMOW', 'R2VSMOW',
           'tiny', 'huge', 'eps']


# Mathematical
Pi    = 3.141592653589793238462643383279502884197    # Pi
pi    = 3.141592653589793238462643383279502884197
Pi2   = 1.57079632679489661923132169163975144209858  # Pi/2
pi2   = 1.57079632679489661923132169163975144209858
Pi3   = 1.0471975511965977461542144610931676280656   # Pi/3
pi3   = 1.0471975511965977461542144610931676280656
TwoPi = 6.283185307179586476925286766559005768394    # 2*Pi
Twopi = 6.283185307179586476925286766559005768394
Sqrt2 = 1.41421356237309504880168872420969807856967  # Sqrt(2)

# Physical
Gravity  = 9.81          # Standard average Earth's gravity [m^2 s^-1]
T0       = 273.15        # Celcius <-> Kelvin [K]
P0       = 101325.       # Standard pressure [Pa]
T25      = 298.15        # Standard ambient temperature [K]
sigma    = 5.67e-08      # Stefan-Boltzmann constant [W m^-2 K^-4]
R        = 8.3144621     # Ideal gas constant [J K^-1 mol^-1]
R_air    = 287.06        # gas constant of dry air [J K^-1 kg^-1]
R_H2O    = 461.4         # gas constant of water vapour [J K^-1 kg^-1]
Na       = 6.02214129e23 # Avogrado number [mol^-1]
REarth   = 6371009.      # Radius of Earth [m]

# Material
mmol_co2 = 44.01         # Molar mass CO2 [g mol^-1]
mmol_h2o = 18.01528      # Molar mass water [g mol^-1]
mmol_air = 28.9644       # Molar mass of dry air [g mol^-1]
# from Cambell G (1985) Soil Physics with BASIC, Elsevier Science
density_quartz = 2.65    # density of quartz [g cm^-3]
cheat_quartz   = 800.    # heat capacity of quartz [J kg^-1 K^-1]
cheat_water    = 4180.   # heat capacity of water [J kg^-1 K^-1]
cheat_air      = 1010.   # heat capacity of air [J kg^-1 K^-1]
latentheat_vaporization = 2.45e6 # latent heat of vaporization of water [J kg^-1]

# Isotope
R13VPDB  = 0.0112372     # 13C isotope ratio of VPDB
R18VSMOW = 2005.2e-6     # 18O isotope ratio of VSMOW
R2VSMOW  = 155.76e-6     # 2H  isotope ratio of VSMOW

# Computational
eps  = np.finfo(np.float).eps
huge = np.finfo(np.float).max
tiny = np.finfo(np.float).tiny


if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.NORMALIZE_WHITESPACE)
