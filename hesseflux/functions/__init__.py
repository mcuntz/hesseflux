"""
Purpose
=======

functions provides a variety of special functions, including
common test functions for parameter estimations such as Rosenbrock and Griewank, 
test functions for parameter sensitivity analysis such as the Ishigami and Homma function,
several forms of the logistic function and its first and second derivatives, and
a variety of functions together with robust and square cost functions to use with
scipy.optimize package.

The module is part of the JAMS Python package
https://github.com/mcuntz/jams_python
It will be synchronised with the JAMS package irregularily if used in other packages.

:copyright: Copyright 2014-2020 Matthias Cuntz, see AUTHORS.md for details.
:license: MIT License, see LICENSE for details.

Subpackages
===========
.. autosummary::
   general_functions
   fit_functions
   logistic_function
   opti_test_functions
   sa_test_functions
"""
from .general_functions   import *
from .fit_functions       import *
from .logistic_function   import *
from .opti_test_functions import *
from .sa_test_functions   import *
