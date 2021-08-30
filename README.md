SORA
====

PACKAGE DESCRIPTION
-------------------

SORA is the acronym for *Stellar Occultation Reduction and Analysis*.
It is a library developed in Python3 with the tools to analyze stellar
occultation data. It is based on Astropy functions and Classes.
**Full documentation at https://sora.readthedocs.io/**

A stellar occultation occurs when a solar system object passes in front
of a star for an observer on Earth, and its shadow causes a temporary
drop in the observed flux of the star. This technique allows the
determination of sizes and shapes with kilometre precision and to obtain
characteristics of the object, such as its albedo, the presence of an
atmosphere, rings, jets, or other structures around the body or even
the detection of topographic features (Sicardy et al. 2011, 2016
Braga-Ribas et al. 2013, 2014, 2019, Dias-Oliveira et al., 2015, 2017,
Benedetti-Rossi et al., 2016, 2019, Ortiz et al., 2015, 2017, 2020,
Leiva et al., 2017, Bérard et al., 2017, Morgado et al., 2019, Gomes-Júnior et al., 2020,
Souami et al., 2020, Santos-Sanz et al., 2021).

SORA is a Python-based, object-oriented library for optimal analysis of
stellar occultation data. The user can use this package to build pipelines
to analyse their stellar occultation’s data. It includes processes starting
on the prediction of such events to the resulting size, shape and position of
the Solar System object. The main modules available at version 0.2
are: **star**, **body**, **observer**, **lightcurve** and
**occultation**. It is important to note that new modules and other
improvements and implementations can be available in future versions.

AUTHORS
-------

Altair R. Gomes-Júnior (1, 2),
Bruno E. Morgado (3, 2, 4),
Gustavo Benedetti-Rossi (1, 3, 2),
Rodrigo C. Boufleur (4, 2),
Flavia L. Rommel (4, 2),
Martin B. Huarca (2, 4)

(1) UNESP - São Paulo State University, Grupo de Dinâmica Orbital e Planetologia, CEP 12516-410, Guaratinguetá, SP 12516-410, Brazil</br>
(2) Laboratório Interinstitucional de e-Astronomia - LIneA and INCT do e-Universo, Rua Gal. José Cristino 77, Rio de Janeiro, RJ 20921-400, Brazil</br>
(3) LESIA, Observatoire de Paris, Université PSL, CNRS, Sorbonne Université, Univ. Paris Diderot, Sorbonne Paris Cité, 5 place Jules Janssen, 92195 Meudon, France</br>
(4) Observatório Nacional/MCTIC, R. General José Cristino 77, Rio de Janeiro, RJ 20.921-400, Brazil</br>

CITATION
--------

If you use SORA in a scientific publication, we would appreciate that you add at your acknowledgement the following statement:

    This research made use of SORA, a python package for stellar occultations reduction and analysis, developed with the support of ERC Lucky Star and LIneA/Brazil.

SYSTEM REQUIREMENTS AND INSTALLATION
------------------------------------

SORA was developed in Python 3.7 and requires the following packages:

-  Astropy (4.3): For astronomical related functions, mainly coordinates and time.

-  Astroquery (0.4.3): To query astronomical database as JPL and Vizier.

-  Matplotlib (3.4.3): For easy and beautiful plots.

-  NumPy (1.21): Otimized mathematical functions.

-  SciPy (1.7.1): Otimized functions for mathematics, science, and engineering.

-  SpiceyPy (4.0.2): SPICE/NAIF functions in python.

-  PyERFA (2.0): Python wrapper for the ERFA library based on the SOFA library.

-  Cartopy (0.19): Geospatial data processing to produce maps.

The user can install SORA and most of its requirements using **pip**, only
Cartopy should be installed from conda afterwards.

```shell
pip install sora-astro
conda install -c conda-forge cartopy
```

If you are a GitHub user, you can also use:

```shell
git clone https://github.com/riogroup/SORA.git
cd SORA
pip install .
conda install -c conda-forge cartopy
```

For a better experience with SORA, we recommend the use of [Jupyter]. The creation of a dedicated Conda environment for SORA is suggested to avoid requirement issues.

Acknowledgements
----------------

The SORA package is hosted on a GitHub repository. It was developed with support
of the LuckyStar, that agglomerates the efforts of the Paris, Granada, and Rio
teams. The LuckyStar is funded by the ERC (European Research Council)
under the European Community’s H2020 (2014-2020/ERC Grant Agreement No. 669416). Also,
this project is supported by LIneA (Laboratório Interinstitucional de e-Astronomia),
INCT do e-Universo (CNPQ grants 465376/2014-2), by FAPESP (proc. 2018/11239-8), by CNPQ
(proc. 300472/2020-0, 150612/2020-6), and by CAPES-PRINT/UNESP (88887.571156/2020-00)
in Brazil.

The Paris, Granada, and Rio teams are professionals astronomers affiliated mainly in the following
institutions:

* LESIA - Observatoire de Paris, France;
* Institut Polytechnique des Sciences Avancées, France;
* IMCCE - Observatoire de Paris, France;
* Instituto de Astrofísica de Andalucía, Spain;
* Laboratório Interinstitucional de e-Astronomia, Brazil;
* INCT do e-Universo, Brazil;
* Observatório Nacional/MCTI, Brazil;
* Federal University of Technology - Paraná, Brazil;
* UNESP - São Paulo State University, Brazil;
* Universidade Federal do Rio de Janeiro - Observatório do Valongo, Brazil;
