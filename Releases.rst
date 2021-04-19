SORA v0.2 (Unreleased)
========================

New Features
------------

sora.body
^^^^^^^^^^^

- Created new Body Class which downloads the occulting body information from online source.
  At the moment, it downloads only from the Small-Body DataBase. The Body class will be the manager
  for all the Body informations, such as Ephem, Shape, Ring, etc. [#51]

- New Class PhysicalData, which inherits from astropy.units.quantity.Quantity, is created to handle
  physical data with uncertainty, reference and notes. [#51]

- "pole_position_angle" and "apparent_magnitude" functions are now present in Body instead of Ephem.[#51]

- Created a hardcoded satellite database to complement missing data of SBDB. It must be replaced in the future. [#61]

sora.config
^^^^^^^^^^^

sora.ephem
^^^^^^^^^^

- A new EphemHorizons was created which is strictly equal to EphemJPL (EphemJPL may be removed in v1.0). [#51]

sora.extra
^^^^^^^^^^

- Allow two ChiSquare objects to be combined into one: `chi3 = chi1 + chi2` [#61]

- New function get_ellipse_points() that calculates the positions on the perimeter of an ellipse [#60]

sora.lightcurve
^^^^^^^^^^^^^^^

sora.observer
^^^^^^^^^^^^^

sora.occultation
^^^^^^^^^^^^^^^^

- A shortcut was created in Occultation where the user can pass the coordinate of the star directly to Occultation,
  the Star object will be created automaticaly. [#46]

- New Chord Class introduced to handle a chord with an Observer and a LightCurve. [#53]

- New ChordList Class introduced to handle the list of Chords in an Occultation. [#53]

- New function .get_impact_param() that calculatesthe impact parameter, minimal distance
  between the chord and the centre position, in Chord and ChordList.[#60]

- New function .get_theoretical_times(),that calcultates the theoretical times and chord size
  for a given ellipse in Chord and ChordList. [#60]

- New fucntion .check_time_shift() that calculates the offset in time to align the center of the chords
  in Occultation. [#60]

- New parameters sigma_result, that saves the result with an extended error bar, and ellipse_error, that
  adds a further systematic error to be considered, in Occultation.fit_ellipse(). [#60]

- New function fiter_negative_chord() that compares the ChiSquare from an Ellipse fitting with the chords
  and remove the solutions that would cross a negative chord [#60]

sora.prediction
^^^^^^^^^^^^^^^

- prediction() now makes use of the user input of the star to calculate faster the occultation parameters. [#48]

- prediction() now can make predictions using Gaia-EDR3. A new parameter "catalogue" was created
  for choosing between Gaia-DR2 and Gaia-EDR3.[#61]

- Fixed bug when plotting the heights in the map in a rotated projection. [#54]

sora.star
^^^^^^^^^^^^^^^

- Star() is now able to fully receive astrometric parameters from the user. [#48]

- Star() is able to download and use the distance from Bailer-Jones et al (2018). [#27]

API Changes
-----------

- Update the argument "log" to "verbose" on all modules. [#61]

sora.config
^^^^^^^^^^^

sora.ephem
^^^^^^^^^^

- "pole_position_angle" and "apparent_magnitude" is passed to Body Class. In Ephem, it will raise
  a FutureWarning. [#51]

- The Ephem classes are now passed through the Body Class which will have priority over Ephem
  attributes. Parameters such as "spkid", "radius", "H" and "G". [#51]

- All Ephem Classes now inherits from BaseEphem, which holds core functionality for all of them. [#51]

sora.extra
^^^^^^^^^^

sora.lightcurve
^^^^^^^^^^^^^^^

- Removed the necessity for LightCurve to have a unique name associated. [#53]

- Cycle time is now determined via mode instead of median. [#56]

sora.observer
^^^^^^^^^^^^^

- Removed the necessity for Observer to have a unique name associated. [#53]

sora.occultation
^^^^^^^^^^^^^^^^

- The new Body Class was implemented in Occultation. For backward compatibility, the previous
  usage is still possible if the Ephem object have a name. The Body Class is only required
  if the object is a planet or a planetary satellite. [#51]

- Deprecated some functions that were passed to ChordList. [#53]

sora.prediction
^^^^^^^^^^^^^^^

- prediction() now creates the time array inside each division to avoid memory overflow. [#48]

- prediction() now propagates the positions of the stars using only the proper motions
  before comparing the stars with the ephemeris. [#48]

- The new Body Class was implemented in prediction. For backward compatibility, the previous
  usage is still possible. [#51]

sora.star
^^^^^^^^^^^^^^^


Bug Fixes
---------

sora.config
^^^^^^^^^^^

sora.ephem
^^^^^^^^^^

sora.extra
^^^^^^^^^^

sora.lightcurve
^^^^^^^^^^^^^^^

- Corrected bug in LightCurve model where the size of the star was being interpreted
  as radius instead of diameter [#60]

sora.observer
^^^^^^^^^^^^^

sora.occultation
^^^^^^^^^^^^^^^^

sora.prediction
^^^^^^^^^^^^^^^

- Fixes issue that happenned in occ_params() when the instant of the occultation was outside the given range.
  The function now gives appropriate error messages. The automatic range search was increased to 50 min
  from central instant in a recursive search. [#45, #48]

sora.star
^^^^^^^^^^^^^^^

- Star now calculates the robust propagation of the position of the star and correspondent uncertainties. [#18]

- Fixed bug in Star().__str__() where pmDEC was printed wrong. [#43]

- A small bug fix was made in Star with the units of the star position error when coordinates are local. [#51]


SORA v0.1.1 (2020/Jul/30)
========================

New Features
------------

sora.config
^^^^^^^^^^^

- Module to verify if kwargs are allowed was created. This was included throughout the code. [#8]

sora.ephem
^^^^^^^^^^

sora.extra
^^^^^^^^^^

- Added a parameter that allows the used to plot a dot corresponding
  the center of the ellipse [#35]

sora.lightcurve
^^^^^^^^^^^^^^^

- Property LightCurve.time_mean that returns the mean time of the chord (positive) or
  the mean time of the observation (negative). [#34]

sora.observer
^^^^^^^^^^^^^

- Function Observer.altaz() that calculates the altitude and azimuth for a given target 
  and instant. [#34]

sora.occultation
^^^^^^^^^^^^^^^^

sora.prediction
^^^^^^^^^^^^^^^

- Four new parameters were added to `plot_occ_map()`: `path`: for the user to select
  a directory where to save the plots; `site_name`: If True, the name of the sites
  will be plotted; `chord_delta` and `chord_geo`: for the user to plot the path of
  a chord from distance of the center or passing by some coordinate, respectively. [#35]

- Two methods were added to `PredictionTable()` to help the user to remove bad events
  from table: `keep_from_selected_images()` and `remove_occ()`. [#35]

sora.star
^^^^^^^^^^^^^^^


API Changes
-----------

sora.config
^^^^^^^^^^^

- config module is now a directory. It now includes a module with decorators
  and another for variables. [#31,#35]

sora.ephem
^^^^^^^^^^

- In EphemKernel, `code` argument was replaced by `spkid`. When using 'code',
  a FutureWarning is raised stating `code` as deprecated and will be removed from v1.0. [#26]

sora.extra
^^^^^^^^^^

sora.lightcurve
^^^^^^^^^^^^^^^

- In LightCurve.immersion and LightCurve.emersion, an error will rise when these values were not 
  instanciated or fitted. [#34]

- Now the user has the possibility to redefine `tref`, `immersion`, `emersion`,
  `initial_time` and `end_time` after instantiated. [#35]

- `lambda_0` argument was replaced by `central_bandpass` and `delta_lambda` by `delta_bandpass`. 
  When using 'lambda_0' or `delta_lambda`, a FutureWarning is raised stating `lambda_0` or `delta_lambda`
  as deprecated and will be removed from v1.0. [#36]

sora.observer
^^^^^^^^^^^^^

sora.occultation
^^^^^^^^^^^^^^^^

- Occultation.new_astrometric_positions() now shows a warning when time is far
  by more than 1 day from the occultation closest approach. [#21]

- Occultation.to_log() and print(Occultation) added the polar radius, equivalent radius, 
  the Sun-Geocenter-Target angle and the Moon-Geocenter-Target angle, geocentric albedo,
  the altitude and azimuth of the target for each Observer. [#17]

- In `fit_ellipse()`, `pos_angle` and `dpos_angle` were deprecated in favor of
  `position_angle` and `dposition_angle`. [#35]

- Changed "GCRS" to "Geocentric" in the string representation to avoid confusion
  about the reference frame. [#35]
  
sora.prediction
^^^^^^^^^^^^^^^

- prediction() now calculates the ephemeris inside each division to avoid memory overflow. [#31]

- PredictionTable.to_ow() will now raise a warning if the radius or the error of
  the ephemeris is not present. [#35]

sora.star
^^^^^^^^^^^^^^^

- Now Star downloads all parameters from Gaia and saves them in the `meta_gaia` attribute [#35]


Bug Fixes
---------

sora.config
^^^^^^^^^^^

sora.ephem
^^^^^^^^^^

- Added function get_position() to EphemPlanete. This corrects a bug that prevented
  Occultation to run with EphemPlanete. [#41]

- Fixed bug in EphemJPL where `id_type` was redefined inside __init__(). [#41]

sora.extra
^^^^^^^^^^

sora.lightcurve
^^^^^^^^^^^^^^^

- Fixed error that appears when the fit was done separately (immersion and emersion times). 
  Now the final model agrees with the fitted values.   [#9]

- Fixed error when the file with the light curve has three columns. [#19]

- Fixed error when the exptime within the LightCurve was set as zero or negative. [#23]

- Fixed error in the automatic mode of LightCurve.normalize(). [#34]

- Fixed bug that was raised in LightCurve.log() when there were no initial or end times
  for lightcurves instantiated with immersion and emersion. [#35]

sora.observer
^^^^^^^^^^^^^

sora.occultation
^^^^^^^^^^^^^^^^

- Corrected error calculation using err = sqrt(star_err^2 + fit_err^2) [#18]

- Occultation.plot_occ_map() now uses the fitted ellipse to calculate the projected shadow radius [#22]

- Corrected bug that raised an error when calling Occultation.get_map_sites()
  and there were no observation added to Occultation. [#31]

- Corrected bug that did not save the fitted params in all occultations when
  more than one occultation was used in fit_ellipse(). [#35]

- Added `axis_labels` and `lw` (linewidth) to Occultation.plot_chords(). [#35]

sora.prediction
^^^^^^^^^^^^^^^

- Fixed error that was generated when only one prediction was found. [#16]

- Fixed error in the output format of PredictionTable.to_ow() when coordinate was positive [#35]

sora.star
^^^^^^^^^^^^^^^


SORA v0.1 - Initial Release (2020/May/20)
=========================================

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

  - Star
 
    - Star Gaia-DR2 ID
    - Star coordinates at 2015.5 and uncertainty - RA and DEC (hh mm ss.sss , +dd mm ss.sss, mas, mas)
    - Star proper motion - in RA, DEC - and uncertainties (mas/yr)
    - Star parallax and uncertainty (mas)
    - Star coordinates propagated to event epoch and uncertainty - RA and DEC (hh mm ss.sss , +dd mm ss.sss, mas, mas)
    - Star magnitudes G, B, V, R, J, H, K (mag)
    - Star projected diameter and model (km and mas, model: GDR2, Van Belle, Kervella)
    - Star offset applied in RA and DEC (mas, mas)


  - Object and Ephemeris

    - Object Name
    - Object radius (km)
    - Object mass (kg)
    - Ephemeris kernel (version and DE)
    - Offset applied in RA/DEC (mas, mas)
    - Object’s distance (AU)
    - Object apparent magnitude for the date (mag)

  - Occultation

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

