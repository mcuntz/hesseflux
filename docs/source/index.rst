==========
Quickstart
==========

``hesseflux`` provides functions used in the processing and
post-processing of the Eddy covariance flux data of the `ICOS
<https://www.icos-cp.eu/>`_ ecosystem site `FR-Hes
<https://www.icos-france.fr/en/static3/the-network>`_.

.. toctree::
   :maxdepth: 3
   :caption: Contents:


About
=====

``hesseflux`` collects functions used for processing Eddy covariance
data of the `ICOS <https://www.icos-cp.eu/>`_ ecosystem site `FR-Hes
<https://www.icos-france.fr/en/static3/the-network>`_.

The package uses several functions of the JAMS Python package

   https://github.com/mcuntz/jams_python

The JAMS package and hesseflux are synchronised irregularly.

``hesseflux`` includes a Python port of Logtools, the Logger Tools
Software of Olaf Kolle, MPI-BGC Jena, (c) 2012.

The post-processing functionality for Eddy flux data is similar to the
R-package
`REddyProc <https://cran.r-project.org/web/packages/REddyProc/index.html>`_
and includes basically the steps described in
`Papale et al. (Biogeosciences, 2006) <https://doi.org/10.5194/bg-3-571-2006>`_
plus some extensions such as the daytime method of flux partitioning
(`Lasslop et al., Global Change Biology 2010 <https://doi.org/10.1111/j.1365-2486.2009.02041.x>`_).

The complete documentation for ``hesseflux`` is available from Read The Docs.

   http://hesseflux.readthedocs.org/en/latest/

   
Quick usage guide
=================

An example script that makes all the steps described in `Papale et
al. (Biogeosciences, 2006) <https://doi.org/10.5194/bg-3-571-2006>`_
is given in the example directory. It is simply called:

.. code-block:: bash
		
    python postproc_europe-fluxdata.py hesseflux_example.cfg

The script is governed by a configuration file in Python's standard
:any:`configparser` format. The example configuration file
`hesseflux_example.cfg` is highly commented. See the `User Guide <userguide.html>`_
for a step by step guide through the script and the configuration
file.


Installation
============

The easiest way to install is via `pip`:

.. code-block:: bash

    pip install hesseflux

See the `installation instructions <install.html>`_ for more information.


License
=======

``hesseflux`` is distributed under the MIT License.  
See the `LICENSE <https://github.com/mcuntz/hesseflux/LICENSE>`_ file for details.

Copyright (c) 2009-2020 Matthias Cuntz

The project structure is based on a `template
<https://github.com/MuellerSeb/template>`_ provided by `Sebastian
Müller <https://github.com/MuellerSeb>`_ .


Contributing to hesseflux
=========================

Users are welcome to submit bug reports, feature requests, and code
contributions to this project through GitHub.

More information is available in the
`Contributing <contributing.html>`_
guidelines.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
