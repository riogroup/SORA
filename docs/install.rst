.. _Sec:install:

.. image:: images/SORA_logo.png
  :width: 500
  :align: center
  :alt: SORA: Stellar Occultation Reduction and Analysis

|
|


Installation
============


Python package requirements
---------------------------

Several SORA functionalities use other Python-based libraries. Below are 
listed the library dependencies and their minimal version needed to use SORA. 
Most of those packages are installed on the fly using the `pip install` 
method, except for Cartopy.


-  `Astropy <https://www.astropy.org/>`_ (4.0): For astronomical related functions, 
   mainly coordinates and time.

-  `Astroquery <https://astroquery.readthedocs.io/en/latest/>`_ (0.4.1): To query 
   astronomical database as JPL and Vizier.

-  `Matplotlib <https://matplotlib.org/>`_ (3.1.1): For easy and beautiful plots.

-  `NumPy <https://numpy.org/>`_ (1.18.1): Otimized mathematical functions.

-  `SciPy <https://www.scipy.org/>`_ (1.4.1): Otimized functions for mathematics, science, and
   engineering.

-  `SpiceyPy <https://spiceypy.readthedocs.io/en/main/>`_ (3.0.2): SPICE/NAIF functions in python.

-  `PyERFA <https://pyerfa.readthedocs.io/en/latest/>`_ (2.0): Python wrapper for the ERFA library based on the SOFA library.   

-  `Cartopy <https://scitools.org.uk/cartopy/docs/latest/>`_ (0.17): Geospatial data processing to produce maps.




Installing SORA
---------------

If you are new to Python or not familiar with Python virtual environments, we 
recommend starting by installing the Anaconda Distribution.  This works on all 
platforms (Linux, macOS, Windows) and installs a full-featured scientific Python 
in a user directory without requiring root permissions. For a better experience 
with SORA, we recommend the use of Jupyter. The creation of a dedicated Conda 
environment for SORA is suggested to avoid requirement issues.

The user can install SORA and most of its requirements using **pip**, only
Cartopy should be installed by hand afterwards.

>>> user@path$> pip install sora-astro
>>> user@path$> conda install -c conda-forge cartopy

If you are a |GitHub| user, you can also use:

>>> user@path$> git clone https://github.com/riogroup/SORA/sora.git
>>> user@path$> cd sora
>>> user@path$> pip install .
>>> user@path$> conda install -c conda-forge cartopy

When new versions are available, the user can update it downloading the
last release from the SORA package in the riogroup organisation on
|GitHubRio|. If you want to be notified just follow the package.

.. |GitHubRio| raw:: html

   <a href="https://github.com/riogroup/SORA" target="_blank"> GitHub</a>

.. |GitHub| raw:: html

   <a href="https://github.com/" target="_blank"> GitHub</a>

Functionalities
---------------

With SORA (v0.2), among other more advanced tasks, the user can easely:

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

