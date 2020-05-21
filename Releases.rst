# SORA

# VERSION RELEASES AND UPDATES


## RELEASES

- SORA v0.1 - 2020/May/20

## SORA v0.1 - Initial Release

### Object Classes

The documentation of all classes and functions are on their docstrings, while the scientific part is presented in the full documentation. Here follows a list with the main Object Classes:

**Ephem** Three Object Classes created to generate geocentric ephemeris for a given solar system object. **EphemJPL** queries the JPL Horizons service and download ephemeris information. **EphemKernel** reads the BSP files to calculate the ephemeris using the Spiceypy package. **EphemPlanet** reads an ASCII file with previously determined positions and interpolate them for a given instant.

JPL Horizons - https://ssd.jpl.nasa.gov/horizons.cgi

**Star** Object Class created to deal with the star parameters. From the Gaia-DR2 Source ID or a sky region, it queries the VizieR service and downloads the star’s information. From Gaia DR2 Catalog (<font color=blue>Gaia Collaboration 2016a, 2016b and 2018</font>) it gets the RA, DEC, parallax, proper motions, G magnitude and star radius; from the NOMAD Catalog (<font color=blue>Zacharias et al. 2004</font>) it gets the B, V, R, J, H and K magnitudes. The user can calculate the ICRS coordinate of the star at any epoch. It can be barycentric (corrected from proper motion) or geocentric (corrected from proper motion and parallax). Also, the apparent diameter of the star is calculated using Gaia DR2 information, or some models such as <font color=blue>Van Belle (1999)</font> and  <font color=blue>Kervella et al. (2004)</font>.

Gaia - Gaia Collaboration 2016a, 2016b and 2018
Mission: https://ui.adsabs.harvard.edu/abs/2016A\%26A...595A...1G/abstract
DR1: https://ui.adsabs.harvard.edu/abs/2016A\%26A...595A...2G/abstract
DR2: https://ui.adsabs.harvard.edu/abs/2018A\%26A...616A...1G/abstract
VizieR - https://vizier.u-strasbg.fr/viz-bin/VizieR
NOMAD - Zacharias et al. 2004
https://ui.adsabs.harvard.edu/abs/2004AAS...205.4815Z/abstract
Van Belle, 1999 - https://ui.adsabs.harvard.edu/abs/1999PASP..111.1515V/abstract
Kervella, 2004 - https://ui.adsabs.harvard.edu/abs/2004A%26A...426..297K/abstract

**Observer**: Object Class created to deal with the observer location. The user can also download the ground-based observatories from the Minor Planet Center (MPC) database.

MPC sites - https://minorplanetcenter.net/iau/lists/ObsCodesF.html

**Light Curve**: Object Class that receives the observational light curve (with time and the occulted star normalized photometry relative to reference stars) and some observational parameters (filter and exposure time). It has functions to determine the instants that the solar system object enters in front of the star and leaves, (immersion and emersion times, respectively). The model considers a sharp-edge occultation model (geometric) convolved with Fresnel diffraction, stellar diameter (projected at the body distance) and finite integration time (<font color=blue>Widemann et al., 2009; Sicardy et al., 2011</font>).

Widemann et al. 2009 -  https://ui.adsabs.harvard.edu/abs/2009Icar..199..458W/abstract
Sicardy et al. 2011 -  https://ui.adsabs.harvard.edu/abs/2011Natur.478..493S/abstract

**Occultation**: Main Object Class within SORA, created to analyze stellar occultations, and control all the other Object Classes within this package. Its functions allow converting the times for each observatory in the occulted body positions in the sky plane relative to the occulted star ($f$, $g$) (<font color=blue>IERS Conventions</font>). Also, to obtain the best ellipse parameters (centre position, apparent equatorial radius, oblateness and the position angle of the apparent polar radius) that fit the points. The results are the apparent size, shape and astrometrical position of the occulting body.

IERS Conventions: https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html

Some extra Objects Classes:

**PredictionTable**: Using the **prediction** function within SORA results in an Object Class that is a slight modification of an AstropyTable. The added changes allow to create the occultation map for each prediction, convert into specific formats, such as OccultWatcher and PRAIA (<font color=blue>Assafin et al. (2011)</font>).

OccultWatcher - https://www.occultwatcher.net/
Assafin et al., 2011 - https://ui.adsabs.harvard.edu/abs/2011gfun.conf...85A/abstract

**ChiSquare**: This Object Class is the result of the fitting functions within SORA, such as _LightCurve.occ_lcfit()_ and _Occultation.fit_ellipse()_. This Class has functions that allow viewing the values that minimize the $\chi^2$ tests, the uncertainties within $n-\sigma$, plotting the tests, and saving the values.   


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


  > Object and Ephemeris

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

