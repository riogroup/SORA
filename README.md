SORA
====

PACKAGE DESCRIPTION
-------------------

SORA is the acronym for *Stellar Occultation Reduction and Analysis*.
It is a library developed in Python3 with the tools to analyze stellar
occultation data. It is based on Astropy functions and Classes.
**Discover the full documentation at https://sora.readthedocs.io/**

A stellar occultation occurs when a solar system object passes in front
of a star for an observer on Earth, and its shadow causes a temporary
drop in the observed flux of the star. This technique allows the
determination of sizes and shapes with kilometre precision and to obtain
characteristics of the object, such as its albedo, the presence of an
atmosphere, rings, jets, or other structures around the body or even
the detection of topographic features (Sicardy et al. 2011, 2016,
Braga-Ribas et al. 2013, 2014, 2019, 2023 Dias-Oliveira et al., 2015, 2017,
Benedetti-Rossi et al., 2016, 2019, Ortiz et al., 2015, 2017, 2020, 2023,
Leiva et al., 2017, Bérard et al., 2017, Morgado et al., 2019, 2022a, 2022b, 2023, 
Gomes-Júnior et al., 2020, Souami et al., 2020, Santos-Sanz et al., 2021, 2023,
Rommel et al., 2020, 2023, Pereira et al., 2023).

SORA is a Python-based, object-oriented library for optimal analysis of
stellar occultation data. The user can use this package to build pipelines
to analyse their stellar occultation’s data. It includes processes starting
on the prediction of such events to the resulting size, shape and position of
the Solar System object. The main modules available at version 0.2
are: **star**, **body**, **observer**, **spacecraft**, **lightcurve**, **chord** 
and **occultation**. It is important to note that new modules and other
improvements and implementations can be available in future versions.

AUTHORS
-------

Altair R. Gomes-Júnior (1, 2, 3),
Bruno E. Morgado (4, 3, 5),
Rodrigo C. Boufleur (5, 3),
Gustavo Benedetti-Rossi (2, 3),
Flavia L. Rommel (6, 5, 3),
Martin B. Huarca (3, 5)
Chrystian L. Pereira (5, 3)
Giuliano Margoti (6, 3)

(1) UFU - Federal University of Uberlândia, Physics Institute, Av. João Naves de Ávila 2121, Uberlândia, MG 38408-100, Brazil</br>
(2) UNESP - São Paulo State University, Grupo de Dinâmica Orbital e Planetologia, Guaratinguetá, SP 12516-410, Brazil</br>
(3) Laboratório Interinstitucional de e-Astronomia - LIneA and INCT do e-Universo, Rua Gal. José Cristino 77, Rio de Janeiro, RJ 20921-400, Brazil</br>
(4) Universidade Federal do Rio de Janeiro - Observatório do Valongo, Ladeira Pedro Antônio 43, CEP 20.080-090, Rio de Janeiro, Brazil</br>
(5) Observatório Nacional/MCTIC, R. General José Cristino 77, Rio de Janeiro, RJ 20.921-400, Brazil</br>
(6) Federal University of Technology - Paraná (UTFPR/DAFIS), Rua Sete de Setembro, 3165, CEP 80230-901, Curitiba, PR, Brazil</br>

CITATION
--------

SORA is a free software, and we kindly ask users that make use of it (or of part of it), especially in projects that results in scientific publications, to include the following citation:

::

    Gomes-Júnior et al. (2022). SORA: Stellar occultation reduction and analysis.
    MNRAS, Volume 511, Issue 1, March 2022, Pages 1167–1181
    DOI: 10.1093/mnras/stac032


The scientific documentation of SORA is described in the paper available on Monthly Notices
of the Royal Astronomical Society.

SYSTEM REQUIREMENTS AND INSTALLATION
------------------------------------

SORA requires Python 3.11 or later within the Python 3 series. Install SORA
from PyPI with:

```shell
python -m pip install sora-astro
```

This command installs all runtime dependencies declared in `pyproject.toml`,
including Cartopy. Cartopy provides binary wheels for the major operating
systems, so a separate Conda installation is not normally required. See the
[Cartopy installation guide](https://cartopy.readthedocs.io/stable/installing.html)
if a compatible wheel is unavailable for your platform.

To install the current source from GitHub:

```shell
git clone https://github.com/riogroup/SORA.git
cd SORA
python -m pip install .
```

For development, install the project in editable mode with the optional test
and documentation dependencies:

```shell
python -m pip install -e ".[test,docs]"
```

Building the documentation also requires the
[Pandoc](https://pandoc.org/installing.html) executable. Complete installation
instructions are available in the
[SORA documentation](https://sora.readthedocs.io/en/latest/install.html).

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
* Federal University of Uberlândia, Brazil;
