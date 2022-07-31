Changelog
---------

All notable changes after its initial release in May 2020 (v2.0) are documented
in this file.

v5.0 (??? 2022)
    * Use all tools from ``pyjams`` rather than local copies.

v4.0 (Jul 2022)
    * Move to new pip structure using `pyproject.toml`.
    * Move to Github actions.
    * Change Markdown to reStructured text.

v3.2 (Jun 2021)
    * Preserve trailing whitespace in float and string read functions.
    * Set error estimates to undef by default.
    * Optionally return mean of values for error estimates to test
      goodness of filling algorithm.

v3.1 (Feb 2021)
    * Support cftime > v1.3.0 by checking Gregorian calendar with existing
      rather than skipped dates.
    * `eddypro2nc.py` script converting EddyPro output to netCDF file.
    * Use build instead of cibuildwheel to make pure Python wheels.

v3.0 (Sep 2020)
    * Included subpackages `const`, `functions`, `logtools` in automatic
      packaging.

v2.2 (Aug 2020)
    * No bootstrap by default in `ustarfilter`.
    * Optionally return u* thresholds per season in `ustarfilter`.
    * Bugfix if no threshold found in `ustarfilter`.
    * Bugfix for muli-year flags of `ustarfilter`.
    * Flag gross fluxes if partitioning was not possible due to missing meteo.
    * Make `gapfill` compatible with flake8.

v2.1 (Jul 2020)
    * Refined documentation.
    * Use python3 for Linux/macOS and python for Windows on TravisCI.
    * Use en instead of eng in date routines.
    * Bugfix for Excel dates in date routines.

v2.0.2 (May 2020)
    * Finished setup in repositories with all dependencies, TWINE password, etc.

v2.0.1 (May 2020)
    * Add more requirements for readthedocs, coverall, etc.

v2.0 (May 2020)
    * Initial public release.
