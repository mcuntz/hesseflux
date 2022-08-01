hesseflux
=========
..
  pandoc -f rst -o README.html -t html README.rst

``hesseflux`` provides functions used in the processing and post-processing of
the Eddy covariance flux data of the ICOS_ ecosystem site FR-Hes_.

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3831488.svg
   :target: https://doi.org/10.5281/zenodo.3831488
   :alt: Zenodo DOI

.. image:: https://badge.fury.io/py/hesseflux.svg
   :target: https://badge.fury.io/py/hesseflux
   :alt: PyPI version

..
   .. image:: https://img.shields.io/conda/vn/conda-forge/pyjams.svg
      :target: https://anaconda.org/conda-forge/pyjams
      :alt: Conda version

.. image:: https://img.shields.io/badge/license-MIT-blue.svg?style=flat
   :target: https://github.com/mcuntz/hesseflux/blob/master/LICENSE
   :alt: License

.. image:: https://github.com/mcuntz/hesseflux/workflows/Continuous%20Integration/badge.svg?branch=master
   :target: https://github.com/mcuntz/hesseflux/actions
   :alt: Build Status

..
   .. image:: https://coveralls.io/repos/github/mcuntz/hesseflux/badge.svg
      :target: https://coveralls.io/github/mcuntz/hesseflux
      :alt: Coverage Status


About hesseflux
---------------

``hesseflux`` collects functions used for processing Eddy covariance data of the
ICOS_ ecosystem site FR-Hes_.

The post-processing functionality for Eddy flux data is similar to the R-package
REddyProc_ and includes basically the steps described in `Papale et al.
(Biogeosciences, 2006)`_ plus some extensions such as the daytime method of flux
partitioning (`Lasslop et al., Global Change Biology 2010`_).


Documentation
-------------

The complete documentation for ``hesseflux`` is available at:

   https://mcuntz.github.io/hesseflux/


Quick usage guide
-----------------

Post-processing Eddy covariance data that is in europe-fluxdata.eu format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

An example script that makes all the steps described in `Papale et al.
(Biogeosciences, 2006)`_ is given in the example directory. It is simply called:

.. code-block:: bash

   python postproc_europe-fluxdata.py hesseflux_example.cfg

The script is governed by a configuration file in Python's standard
*configparser* format. The example configuration file *hesseflux_example.cfg* is
highly commented. See the `User Guide`_ for a step by step guide through the
script and the configuration file.


Installation
------------

The easiest way to install is via `pip`:

.. code-block:: bash

   pip install hesseflux

..
   or via `conda`:

   .. code-block:: bash

      conda install -c conda-forge hesseflux

Requirements
    * numpy_
    * scipy_
    * pandas_

    ..
       * pyjams_


License
-------

``hesseflux`` is distributed under the MIT License. See the LICENSE_ file for
details.

Copyright (c) 2009-2022 Matthias Cuntz

The project structure ``hesseflux`` has borrowed heavily from welltestpy_
by `Sebastian Müller`_.

.. _ICOS: https://www.icos-cp.eu/
.. _FR-Hes: https://www.icos-france.fr/en/static3/the-network
.. _REddyProc: https://cran.r-project.org/web/packages/REddyProc/index.html
.. _Papale et al. (Biogeosciences, 2006): https://doi.org/10.5194/bg-3-571-2006
.. _Lasslop et al., Global Change Biology 2010: https://doi.org/10.1111/j.1365-2486.2009.02041.x
.. _User Guide: https://mcuntz.github.io/pyjams/html/userguide.html
.. _numpy: https://numpy.org/
.. _scipy: https://scipy.org/
.. _pandas: https://pandas.pydata.org/
.. _pyjams: https://github.com/mcuntz/pyjams/
.. _welltestpy: https://github.com/GeoStat-Framework/welltestpy/
.. _Sebastian Müller: https://github.com/MuellerSeb
