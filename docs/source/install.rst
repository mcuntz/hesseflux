============
Installation
============

The easiest way to install ``hesseflux`` is via ``pip``:

.. code-block:: bash

    pip install hesseflux


Manual install
--------------

The latest version of ``hesseflux`` can be installed from source:

.. code-block:: bash

    git clone https://github.com/mcuntz/hesseflux.git
    cd hesseflux
    pip install .


Local install
-------------

Users without proper privileges can append the `--user` flag to
``pip`` either while installing from the Python Package Index (PyPI):

.. code-block:: bash

    pip install hesseflux --user

or from the top ``hesseflux`` directory:

.. code-block:: bash

    git clone https://github.com/mcuntz/hesseflux.git
    cd hesseflux
    pip install . --user

If ``pip`` is not available, then `setup.py` can still be used:

.. code-block:: bash

    python setup.py install --user

When using `setup.py` locally, it might be that one needs to append `--prefix=`
to the command:

.. code-block:: bash

    python setup.py install --user --prefix=

    
Dependencies
------------

``hesseflux`` uses the packages :mod:`numpy`, :mod:`scipy` and :mod:`pandas`.
They are all available in PyPI and ``pip`` should install them automatically.
Installations via `setup.py` might need to install the three dependencies first.
