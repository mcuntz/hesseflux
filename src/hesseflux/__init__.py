"""
hesseflux provides functions used in the processing and
post-processing of the Eddy covariance flux data

It was developed for the ICOS ecosystem site FR-Hes.

The post-processing functionality for Eddy flux data is similar
to the R-package REddyProc and includes basically the steps
described in Papale et al. (Biogeosciences, 2006) plus some
extensions such as the daytime method of flux partitioning
(Lasslop et al., Global Change Biology 2010).

:copyright: Copyright 2009-2022 Matthias Cuntz, see AUTHORS.rst for details.
:license: MIT License, see LICENSE for details.

Subpackages
===========
.. autosummary::
   fgui
   gapfill
   madspikes
   nee2gpp
   ustarfilter
   logtools

History
    * Written 2017 by Matthias Cuntz (mc (at) macu (dot) de)
    * v2.0, format and docu useable with PyPI, Apr 2020, Matthias Cuntz
    * v2.0.1, more requirements for readthedocs, coveralls, etc.,
      May 2020, Matthias Cuntz
    * v2.0.2, finished setup with all dependencies, setting passwords, etc.,
      May 2020, Matthias Cuntz
    * v2.1, use en instead of eng and bugfix for Excel dates in date routines,
      Jul 2020, Matthias Cuntz
    * v2.2, seasonout in ustarfilter, flag gross fluxes if partitioning was
      not possible, Aug 2020, Matthias Cuntz
    * v3.0, include subpackages in const, functions, logtools in setup.py,
      Sep 2020, Matthias Cuntz
    * v3.1, support cftime>v1.3.0 and added eddy2nc,
      Feb 2021, Matthias Cuntz
    * v3.2, return also mean of values for error estimates,
      Jun 2021, Matthias Cuntz
    * v3.2.1, code refactoring,
      Feb 2022, Matthias Cuntz
    * v4.0, moved to pyproject.toml structure and github pages for
      documentation, Jul 2022, Matthias Cuntz
    * v5.0, use pyjams, remove old modules, Aug 2022, Matthias Cuntz

"""
# version, author
try:  # pragma: no cover
    from ._version import __version__
except ImportError:  # pragma: no cover
    # package is not installed
    __version__ = "0.0.0.dev0"
__author__  = "Matthias Cuntz"

from .fgui import directory_from_gui, directories_from_gui
from .fgui import file_from_gui, files_from_gui
from .gapfill import gapfill
from .madspikes import madspikes
from .nee2gpp import nee2gpp
from .ustarfilter import ustarfilter

from . import logtools
