"""
Purpose
=======

logtools is a Python port of the Control File Functions of
Logtools, the Logger Tools Software of Olaf Kolle, MPI-BGC Jena, (c) 2012.

From the Logtools manual:
"The functions range from simple mathematic operations to more complex
and special procedures including functions for checking data. Most of
the functions have the following appearance: `a = f(b,p1,p2,...,pn)`
where `a` is the variable in which the result of the function `f` is
stored, `b` is the input variable of the function and `p1` to `pn` are
parameters (numbers) of the function. An output variable (result of a
function) may be the same as an input variable. Some functions need
more than one input variable, some functions do not need any parameter
and some functions (`mean`, `mini`, `maxi`) may have a variable number
of input variables."

The module is part of the JAMS Python package
https://github.com/mcuntz/jams_python
It will be synchronised with the JAMS package irregularily if used in other packages.

:copyright: Copyright 2014-2020 Matthias Cuntz, see AUTHORS.md for details.
:license: MIT License, see LICENSE for details.

Subpackages
===========
.. autosummary::
   logtools
"""
from .logtools import *
