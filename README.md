# Read-me <=> SORA

# ABOUT
SORA is the acronym for Stellar Occultation Reduction and Analysis.

It is a progam developed in Python to analyse stellar occultation data. It works as an update, automation and improvement of lots of the routines and other programs in current use.

# AUTHORS

Altair R. Gomes-Júnior - LIneA; UNESP Guaratinguetá
Bruno E. Morgado - LIneA; Observatoire de Paris / Meudon
Gustavo Benedetti-Rossi - LIneA; Observatoire de Paris / Meudon
Rodrigo C. Boufleur - LIneA; Observatório Nacional

# PROGRAM DESCRIPTION

SORA is a Python based, object oriented program to update, automate and improve Bruno Sicardy's fortran routines used for data reduction obtained in stellar occultations observations. It also include programs developed by other collaborators. It is important to stress that SORA is NOT a "translation" from fortran to Python, but it is a reprogramming of all models and scientific cases used on those programs.

SORA includes the following programs:
  - polyfit.f
  - ephem.f90 (NIMA)
  - ephem_planete.f
  - fit_d2_ksi_eta.f
  - diam.f
  - gstar.f
  - bar.f (bar2.f)
  - gene_bar.f
  - positionv.f
  - ellipse_fit.f
  - gene_ellipse.f
  - StelitoCalculus.c
  
SORA, being developed in Python, also facilitates the current method for plotting the results from the programs.

SORA is developed in modules  (see VERSION AND UPDATES section) to perform the complete data reduction, pre- and post-occultation. Each of them have their own set of input and output values (I/O). It is important to note that those I/O can be the same for different modules or even some outputs from a module will work as input for other, meaning that they are integrated and has an inside communication.

This program was also developed in order to make easy adaptations for an integration with TNO portal at LIneA.


# PROGRAM PIPELINE

For each observation, a light curve will be associated, as well as the information about the observer, the observational site, and the equipment used (telescope and camera). SORA will use the light curve to calculate the ingress and egress times and project them to the tangent sky plane, where other information about the object ephemeris and the occulted star will be added. The projection will lead to chords that will be used to obtain the object apparent shape ate the moment of the occultation. In every process the user has full control over the parameters and processes.

# SYSTEM REQUIREMENTS

SORA was developed in Python 3.7 and requires the following packages:

  - astropy 4.0: For astronomical related functions, mainly coordinates and time.
  - astroquery 0.4.1: To query astronomical database as JPL and Vizier.
  - cartopy 0.17} Geospatial data processing to produce maps.
  - matplotlib 3.1.1: For easy and beautiful plots.
  - numpy 1.18.1: Otimized mathematical functions.
  - scipy 1.4.1: Otimized functions for mathematics, science, and engineering. 
  - spiceypy 3.0.2: SPICE/NAIF functions in python.
  



# VERSION RELEASES AND UPDATES

The full documentation of all modules and fuctions are on docstrings, while the scientific part are presented in  **REFERENCE** .


## RELEASES

  - SORA v0.1 - 2020/May/18



## SORA v0.1

### MODULES

Module Prediction: Module to make the prediction of stellar occultation. It uses the object ephemeris (plus its updates and offsets) and a star catalogue to create prediction maps of the candidate events. The maps present the object shadow path on Earth and contain the occultation related information.

Light Curve fit (LC fit): Module to determine the ingress and egress instants from an observer data-set (the light curve). It uses the light curve (with time and normalised relative flux between the target star and one or more reference stars) to calculate the instants that the object enters in front of the star and leave (immersion and emersion instants, respectively).

Sky Projection: Module to project each of the times obtained with the module LC fit into positions (f, g) in the sky plane for a Geocentric reference frame. It uses the the times obtained in module LC fit and the occultation information from module Prediction, added to observer information, to convert the instants into positions in the sky plane

Ellipse fit: Module to fit an ellipse with given points. It uses the projected positions in the sky plane (f, g) obtained in module Sky Projection to fit an ellipse. 


#### INPUTS AND OUTPUTS v0.1

##### INPUTS
  - **Event Related (Star and Ephem)**
  
    - Object Name or provisory designation
    - Object Code (only for EphemKernel)
    - BSP file and name (only for EphemKernel)
    - DE file and name (only for EphemKernel)
    - Ephemeris offset for RA and DEC - $\Delta \alpha \cdot \cos \delta$, $\Delta \delta$ (set as 0,0)
    - Occultation date and time
    - Occulted star coordinates RA and DEC; or Gaia code
    - Star offset for RA and DEC - $\Delta \alpha \cdot \cos \delta$, $\Delta \delta$ (set as 0,0)

  - **Observer Related**
  
    - Site name and location (latitude, longitude, and height; or IAU/MPC code)
    - Light curve file and name; or array with fluxes and times; or immersion and emersion times
    - Exposure time in seconds
    - Observational bandwidth in microns (set as 0.7 $\pm$ 0.3 microns, Clear)

  - **Fitting Related**
  
    - Initial guess for light curve fitting: immersion, emersion and opacity.
    - Range to explore all three parameters
    - Initial guess for ellipse parameters: center (f,g), equatorial radius, oblateness, and position angle
    - Range to explore all five parameters


##### OUTPUTS

  - Generic Information (for a specific occultation):
  
  > Star
  
      - Star Gaia-DR2 ID
      - Star coordinates at 2015.5 and uncertainty - RA and DEC (hh mm ss.sss , +dd mm ss.sss, mas, mas)
      - Star proper motion - in RA, DEC - and uncertainties (mas/yr)
      - Star parallax and uncertainty (mas) 
      - Star coordinates propagated to event epoch and uncertainty - RA and DEC (hh mm ss.sss , +dd mm ss.sss, mas, mas)
      - Star magnitudes G, B, V, R, J, H, K (mag)
      - Star projected diameter and model (km and mas, model: GDR2, Van Belle, Kervella)
      - Star offset applied in RA and DEC (mas, mas)


  > Object and Ephmerid

      - Object Name
      - Object radius (km)
      - Object mass (kg)
      - Ephemeris kernel (version and DE)
      - Offset applied in RA/DEC (mas, mas)
      - Object’s distance (AU)
      - Object apparent magnitude for the date (mag)

  > Occultation

      - Event date and time (yyyy-mm-dd hh:mm:ss.sss)
      - Closest approach Angle - CA (arcsec)
      - Reference time (yyyy-mm-dd hh:mm:ss.sss)
      - Position Angle - PA (degree)
      - Shadow’s velocity relative to the geocenter (km/s)
      - Number of positive observations
      - Number of negative observations


  - Observer Information
 
      - Detection status (positive, negative, overcast, tech. problem, other)
      - Site Name
      - Site MPC/IAU code (if any)
      - Site coordinates - Latitude, Longitude and height  (dd mm ss.s ; dd mm ss.s ; m)
      - Light curve file name
      - Number of images (lines in LC)

  - Light curve fitting information (for each positive detection)

      - Acquisition start time (yyyy-mm-dd hh:mm:ss.sss)
      - Acquisition end time (yyyy-mm-dd hh:mm:ss.sss)
      - Exposure time (s)
      - Cycle time (s)
      - Time offset applied in LC (s)
      - Light curve calculated RMS
      - Calculated normalised flux and bottom flux (standard = 1, 0)
      - Band width and uncertainty (microns)
      - Shadow's velocity relative to the station (km/s)
      - Fresnel scale (s and km)
      - Projected stellar size scale (s and km)
      - Integration time scale (s and km)
      - Dead time scale (s and km)
      - Model resolution - size of synthetic LC point (s and km)
      - Immersion Time and uncertainty (yyyy-mm-dd hh:mm:ss.sss +/- s.sss)
      - Immersion Time and uncertainty - 1$\sigma$ and 3$\sigma$ (s)
      - Emersion Time and uncertainty (yyyy-mm-dd hh:mm:ss.sss +/- s.sss)
      - $\chi^2$ fit model
      - Emersion Time and uncertainty - 1$\sigma$ and 3$\sigma$ (s)
      - Minimum Chi-square - $\chi^2_{min}$
      - Number of fitted points for im- and emersion
      - Number of fitted parameters
      - Minimum Chi-square per degree of freedom - $\chi^2_{min-pdf}$

  - Elipse fit procedure
  
      - Fitted parameters: Equatorial radius and uncertainty (km); Center position ($f_0$, $g_0$) and 1$\sigma$ uncertainties (km, km); Oblateness and uncertainty; Position angle and uncertainty (degree)
      - Minimum Chi-square -  $\chi^2_{min}$
      - Minimum Chi-square per degree of freedom - $\chi^2_{min-pdf}$ 
      - Number points used to fit ( X points from Y chords )
      - Astrometric object center position at occ. time and uncertainty (hh mm ss.sss +dd mm ss.sss $\pm$ mas)

  - Plots and files (some are optional)

      - Prediction map (Lucky Star model)
      - Normalised light curve - for each site (x = time; y = flux)
      - Chi-square map for immersion and emersion times (x = time; y = $\chi^2$)
      - Light curve and synthetic LC- for each site (x = time; y = flux)
      - Chords projected in sky plane (x = $\xi$ (km); y = $\eta$ (km) )
      - Chi-square map for each ellipse parameter (x = time; y = $\chi^2_{param}$)
      - Chords projected in sky plane and the best ellipse fitted with 1$\sigma$ uncertainties (x = $\xi$ (km); y = $\eta$ (km) )
      - Log file with all information



# CURIOSITIES

- SORA means "sky" in Japanese
- Most of the v0.1 code was developed during the lockdown due to the COVID19 pandemic
