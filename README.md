# SORA

# ABOUT
SORA is the acronym for _Stellar Occultation Reduction and Analysis_. It is a package developed in Python3 with the tools to analyze stellar occultation data. It is based on Astropy functions and Objects. **Full documentation at (<font color=blue>link</font>)**

# AUTHORS

Altair R. Gomes-Júnior (1, 2) - altair.gomes@unesp.br or altair.gomes@linea.gov.br
Bruno E. Morgado (3, 2) - Bruno.Morgado@obspm.fr or Morgado.fis@gmail.com
Gustavo Benedetti-Rossi (3, 2) - Gustavo.Benedetti-Rossi@obspm.fr or gugabrossi@gmail.com
Rodrigo C. Boufleur (4, 2) - rodrigo.boufleur@linea.gov.br or rcboufleur@gmail.com

(1) UNESP - São Paulo State University, Grupo de Dinâmica Orbital e Planetologia, CEP 12516-410, Guaratinguetá, SP 12516-410, Brazil
(2) Laboratório Interinstitucional de e-Astronomia - LIneA and INCT do e-Universo, Rua Gal. José Cristino 77, Rio de Janeiro, RJ 20921-400, Brazil
(3) LESIA, Observatoire de Paris, Université PSL, CNRS, Sorbonne Université, Univ. Paris Diderot, Sorbonne Paris Cité, 5 place Jules Janssen, 92195 Meudon, France
(4) Observatório Nacional/MCTIC, R. General José Cristino 77, Rio de Janeiro, RJ 20.921-400, Brazil

# PACKAGE DESCRIPTION

A stellar occultation occurs when a solar system object passes in front of a star for an observer on Earth, and its shadow causes a temporary drop in the observed flux of the star. This technique allows the determination of sizes and shapes with kilometre precision and to obtain characteristics of the object, such as its albedo, the presence of an atmosphere, rings, jets, or other structures around the body (<font color=blue>Sicardy et al. 2011, 2016; Braga-Ribas et al. 2013, 2014, 2019; Dias-Oliveira et al. 2015; Benedetti-Rossi et al. 2016, 2019; Ortiz et al. 2015, 2017; Leiva et al. 2017; Bérard et al. 2017, Morgado et al. 2019, Gomes-Júnior et al., 2016, 2020</font>), or even the detection of topographic features (<font color=blue>Dias-Oliveira et al. 2017</font>).

SORA is a Python-based, object-oriented package for optimal analysis of stellar occultation data. It includes processes starting on the prediction of such events to the resulting size, shape and position of the Solar System object. The main object classes created are: **Star**, **Ephem**, **Observer**, **LightCurve** and **Occultation**. 




A stellar occultation is defined by the occulting body (**Ephem**), the occulted star (**Star**), and the time of the occultation. On the other hand, each observational station (**Observer**) will be associated with their light curve (**LightCurve**). SORA has tasks that allow the user to determine the immersion and emersion times and project them to the tangent sky plane, using the information within the **Observer**, **Ephem** and **Star** Objects. That projection will lead to chords that will be used to obtain the object's apparent size, shape and position at the moment of the occultation. Automatic processes were developed for optimizing the reduction of typical events. However, users have full control over the parameters and methods and can make changes in every step of the processes.

The Object Classes developed within SORA were built as an integrated system, and they are controlled by the **Occultation** Object Class. Combining these Object Classes and the functions therein, the user can perform the complete data reduction, pre- and post-occultation. Jupyter-Notebooks with the features of each Class can be found in the example folder.

This package was developed with support of the ERC Lucky Star, that agglomerates the efforts of the Paris, Granada, and Rio teams. The Lucky Star is funded by the ERC (European Research Council) under the European Community’s H2020 (2014-2020/ERC Grant Agreement No. 669416). Also, this project is supported by LIneA (_Laboratório Interinstitucional de e-Astronomia_), INCT do e-Universo (CNPQ grants 465376/2014-2), by FAPESP (proc. 2018/11239-8) and by CNPQ (proc. 300472/2020-0), Brazil.

# CITATION AND ACKNOWLEDGEMENT
If you use SORA in a scientific publication, we would appreciate that you add at your acknowledgement the following statement:
“This research made use of SORA, a python package for stellar occultations reduction and analysis, developed with the support of ERC Lucky Star and LIneA/Brazil.”

# SYSTEM REQUIREMENTS AND INSTALLATION

SORA was developed in Python 3.7 and requires the following packages:

  - astropy 4.0: For astronomical related functions, mainly coordinates and time.
  - astroquery 0.4.1: To query astronomical databases and services as JPL/Horizon and VizieR.
  - cartopy 0.17: Geospatial data processing to produce maps.
  - matplotlib 3.1.1: For easy and beautiful plots.
  - numpy 1.18.1: Optimized mathematical functions.
  - scipy 1.4.1: Optimized functions for mathematics, science, and engineering.
  - spiceypy 3.0.2: SPICE/NAIF functions in python.

If you have a git account

  - $ git clone git@github.com:riogroup/SORA.git
  - $ cd SORA
  - $ pip install .

If you do not have a git account

  - $ git clone https://github.com/riogroup/SORA.git
  - $ cd SORA
  - $ pip install .

# CURIOSITIES

- SORA means “sky” in Japanese and “image” in Arabic.
- Most of the v0.1 code was developed during the lockdown due to the COVID19 pandemic.

