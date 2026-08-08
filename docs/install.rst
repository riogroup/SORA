.. _Sec:install:

.. image:: images/SORA_logo.png
  :width: 500
  :align: center
  :alt: SORA: Stellar Occultation Reduction and Analysis

|
|


Installation
============


Python and runtime dependencies
-------------------------------

SORA requires Python 3.11 or later within the Python 3 series. The runtime
dependencies and their minimum versions are declared in ``pyproject.toml`` and
installed automatically with SORA:

- `NumPy <https://numpy.org/>`_ 2.2 or later;
- `PyERFA <https://pyerfa.readthedocs.io/en/latest/>`_ 2.0 or later;
- `Astropy <https://www.astropy.org/>`_ 7.0 or later;
- `jplephem <https://pypi.org/project/jplephem/>`_ 2.13 or later;
- `Astroquery <https://astroquery.readthedocs.io/en/latest/>`_ 0.4.9 or later;
- `SpiceyPy <https://spiceypy.readthedocs.io/en/main/>`_ 6.0.0 or later;
- `Matplotlib <https://matplotlib.org/>`_ 3.10.0 or later;
- `SciPy <https://scipy.org/>`_ 1.15 or later;
- `Requests <https://requests.readthedocs.io/>`_;
- `tqdm <https://tqdm.github.io/>`_ 4.66 or later;
- `Shapely <https://shapely.readthedocs.io/en/stable/>`_ 2.0.7 or later;
- `Cartopy <https://cartopy.readthedocs.io/stable/>`_ 0.24 or later.

Cartopy 0.22 and later provides binary wheels for the major operating systems,
so it is normally installed directly by ``pip`` with the other dependencies.
If a compatible wheel is unavailable for your platform, consult the
`Cartopy installation guide <https://cartopy.readthedocs.io/stable/installing.html>`_.

Installing SORA
---------------

Using a dedicated virtual environment is recommended to avoid dependency
conflicts. Install the latest SORA release and all its runtime dependencies
from PyPI with:

.. code-block:: console

   $ python -m pip install sora-astro

To update an existing installation:

.. code-block:: console

   $ python -m pip install --upgrade sora-astro

To install the current source from GitHub:

.. code-block:: console

   $ git clone https://github.com/riogroup/SORA.git
   $ cd SORA
   $ python -m pip install .

Development and documentation dependencies
------------------------------------------

The optional dependency groups defined in ``pyproject.toml`` can be installed
with the source checkout. For an editable development installation with the
future test suite and documentation tools, use:

.. code-block:: console

   $ python -m pip install -e ".[test,docs]"

The ``docs`` extra installs the Python packages used by Sphinx, including
``nbsphinx``, but the documentation build also requires the
`Pandoc <https://pandoc.org/installing.html>`_ executable. After installing
Pandoc and tox, build the documentation with:

.. code-block:: console

   $ python -m pip install tox
   $ tox run -e build_docs

Functionalities
---------------

With SORA, among other more advanced tasks, the user can easely:

#. Predict stellar occultations and obtain predictions maps;
#. Check when a stellar occultation will happen for a given observer;
#. Analyse occultation light curves and determine the immersion and 
   emersion times for the event;
#. Plot and check the chords in the skyplane;
#. Fit a circle for events with less than 3 chords or an ellipse for 
   events with more chords;
#. Determine the astrometric position of the occulting object, its 
   apparent size and projected shape.

**All these steps can be found in our Jupyter-Notebooks Tutorials.**
