"""
Purpose
=======

hesseflux provides functions used in the processing and
post-processing of the Eddy covariance flux data of the
ICOS ecosystem site FR-Hes.

The package uses several functions of the JAMS Python package
https://github.com/mcuntz/jams_python
The JAMS package and hesseflux are synchronised irregularly.

hesseflux includes a Python port of Logtools, the Logger Tools
Software of Olaf Kolle, MPI-BGC Jena, (c) 2012.

The post-processing functionality for Eddy flux data is similar
to the R-package REddyProc and includes basically the steps
described in Papale et al. (Biogeosciences, 2006) plus some
extensions such as the daytime method of flux partitioning
(Lasslop et al., Global Change Biology 2010).

:copyright: Copyright 2009-2020 Matthias Cuntz, see AUTHORS.md for details.
:license: MIT License, see LICENSE for details.

Subpackages
===========
.. autosummary::
   const
   functions
   argsort
   ascii2ascii
   date2dec
   dec2date
   division
   esat
   fgui
   fread
   fsread
   gapfill
   mad
   madspikes
   nee2gpp
   sread
   ustarfilter
   logtools
"""
from . import const
from . import functions

from .argsort     import argmax, argmin, argsort
from .ascii2ascii import ascii2ascii, ascii2en, ascii2fr, ascii2us, ascii2eng
from .ascii2ascii import en2ascii, fr2ascii, us2ascii, eng2ascii
from .date2dec    import date2dec
from .dec2date    import dec2date
from .division    import division, div
from .esat        import esat
from .fgui        import directory_from_gui, directories_from_gui, file_from_gui, files_from_gui
from .fread       import fread
from .fsread      import fsread
from .gapfill     import gapfill
from .mad         import mad
from .madspikes   import madspikes
from .nee2gpp     import nee2gpp
from .sread       import sread
from .ustarfilter import ustarfilter

from . import logtools

from .version import __version__, __author__
