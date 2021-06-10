SORA v0.2 (2021/Jun/14)
=======================

New Features
------------

sora.body
^^^^^^^^^^^

- Created new Body Class which downloads the occulting body information from online source.
  At the moment, it downloads only from the Small-Body DataBase. The Body class will be the manager
  for all the Body information, such as Ephem, Shape, Ring, etc. [:issue:`51`]

- New Class PhysicalData, which inherits from astropy.units.quantity.Quantity, is created to handle
  physical data with uncertainty, reference and notes. [:issue:`51`]

- "pole_position_angle" and "apparent_magnitude" functions are now present in Body
  instead of Ephem.[:issue:`51`]

- Created a hardcoded satellite database to complement missing data of SBDB. It must be
  replaced in the future. [:issue:`61`]

sora.ephem
^^^^^^^^^^

- A new EphemHorizons was created which is strictly equal to EphemJPL (EphemJPL may be removed in v1.0). [:issue:`51`]

- A new function that downloads the kernel from JPL was added. [:issue:`33`]

sora.extra
^^^^^^^^^^

- Allow two ChiSquare objects to be combined into one: `chi3 = chi1 + chi2` [:issue:`61`]

- New function get_ellipse_points() that calculates the positions on the perimeter of an ellipse [:issue:`60`]

sora.observer
^^^^^^^^^^^^^

- New Spacecraft class developed to handle the geometry of a spacecraft observation.
  To use it,it is necessary a spkid and ephemeris. Ex:
  `spacecraft = Spacecraft(name='New Horizons', spkid='-98', ephem='horizons')`. [:issue:`63`]

- The Observer class was updated to have an ephemeris as well. [:issue:`63`]

- Now the observer can be passed as parameter to `Ephem*.get_position(observer=observer)`,
  `Star.get_position()`, `Body.get_pole_position_angle()` and `Body.apparent_magnitude()`. [:issue:`63`]

sora.occultation
^^^^^^^^^^^^^^^^

- A shortcut was created in Occultation where the user can pass the coordinate of the star directly to Occultation,
  the Star object will be created automatically. [:issue:`46`]

- New Chord Class introduced to handle a chord with an Observer and a LightCurve. [:issue:`53`]

- New ChordList Class introduced to handle the list of Chords in an Occultation. [:issue:`53`]

- New function .get_impact_param() that calculates the impact parameter, minimal distance
  between the chord and the centre position, in Chord and ChordList.[:issue:`60`]

- New function .get_theoretical_times(), that calculates the theoretical times and chord size
  for a given ellipse in Chord and ChordList. [:issue:`60`]

- New function .check_time_shift() that calculates the offset in time to align the center of the chords
  in Occultation. [:issue:`60`]

- New parameters sigma_result, that saves the result with an extended error bar, and ellipse_error, that
  adds a further systematic error to be considered, in Occultation.fit_ellipse(). [:issue:`60`]

- New function filter_negative_chord() that compares the ChiSquare from an Ellipse fitting with the chords
  and remove the solutions that would cross a negative chord [:issue:`60`]

- New method to calculate the "f" and "g" positions for observers without referring to the geocenter. [:issue:`63`]

sora.prediction
^^^^^^^^^^^^^^^

- prediction() now makes use of the user input of the star to calculate faster the occultation parameters. [:issue:`48`]

- prediction() now can make predictions using Gaia-EDR3. A new parameter "catalogue" was created
  for choosing between Gaia-DR2 and Gaia-EDR3.[:issue:`61`]

- Fixed bug when plotting the heights in the map in a rotated projection. [:issue:`54`]

- prediction() can now predict for any observer. Ex: `prediction(..., reference_center=observer)`. [:issue:`63`]

sora.star
^^^^^^^^^^^^^^^

- A new method get_position() was implemented in Star() that will replace geocentric()
  and barycentric() methods [:issue:`63`]

API Changes
-----------

- Update the argument "log" to "verbose" on all modules. [:issue:`61`]

sora.ephem
^^^^^^^^^^

- "pole_position_angle" and "apparent_magnitude" is passed to Body Class. In Ephem, it will raise
  a FutureWarning. [:issue:`51`]

- The Ephem classes are now passed through the Body Class which will have priority over Ephem
  attributes. Parameters such as "spkid", "radius", "H" and "G". [:issue:`51`]

- All Ephem Classes now inherits from BaseEphem, which holds core functionality for all of them. [:issue:`51`]

sora.lightcurve
^^^^^^^^^^^^^^^

- Removed the necessity for LightCurve to have a unique name associated. [:issue:`53`]

- Cycle time is now determined via mode instead of median. [:issue:`56`]

sora.observer
^^^^^^^^^^^^^

- Removed the necessity for Observer to have a unique name associated. [:issue:`53`]

sora.occultation
^^^^^^^^^^^^^^^^

- The new Body Class was implemented in Occultation. For backward compatibility, the previous
  usage is still possible if the Ephem object have a name. The Body Class is only required
  if the object is a planet or a planetary satellite. [:issue:`51`]

- Deprecated some functions that were passed to ChordList. [:issue:`53`]

sora.prediction
^^^^^^^^^^^^^^^

- prediction() now creates the time array inside each division to avoid memory overflow. [:issue:`48`]

- prediction() now propagates the positions of the stars using only the proper motions
  before comparing the stars with the ephemeris. [:issue:`48`]

- The new Body Class was implemented in prediction. For backward compatibility, the previous
  usage is still possible. [:issue:`51`]


Bug Fixes
---------

sora.lightcurve
^^^^^^^^^^^^^^^

- Corrected bug in LightCurve model where the size of the star was being interpreted
  as radius instead of diameter [:issue:`60`]

sora.prediction
^^^^^^^^^^^^^^^

- Fixes issue that happened in occ_params() when the instant of the occultation was outside the given range.
  The function now gives appropriate error messages. The automatic range search was increased to 50 min
  from central instant in a recursive search. [:issue:`45, 48`]


SORA v0.1.2 (2020/Dec/14)
=========================

New Features
------------

sora.star
^^^^^^^^^^^^^^^

- Star() is now able to fully receive astrometric parameters from the user. [:issue:`48`]

- Star() is able to download and use the distance from Bailer-Jones et al (2018). [:issue:`27`]

- Gaia-EDR3 was implemented in Star() and is now a default feature. [:issue:`52`]


API Changes
-----------

sora.star
^^^^^^^^^^^^^^^

- The star module was moved to its own directory. [:issue:`52`]


Bug Fixes
---------

sora.star
^^^^^^^^^^^^^^^

- Star now calculates the robust propagation of the position of the star and correspondent uncertainties. [:issue:`18`]

- Fixed bug in Star().__str__() where pmDEC was printed wrong. [:issue:`43`]

- A small bug fix was made in Star with the units of the star position error when coordinates are local. [:issue:`51`]


SORA v0.1.1 (2020/Jul/30)
=========================

New Features
------------

sora.config
^^^^^^^^^^^

- Module to verify if kwargs are allowed was created. This was included throughout the code. [:issue:`8`]

sora.extra
^^^^^^^^^^

- Added a parameter that allows the used to plot a dot corresponding
  the center of the ellipse [:issue:`35`]

sora.lightcurve
^^^^^^^^^^^^^^^

- Property LightCurve.time_mean that returns the mean time of the chord (positive) or
  the mean time of the observation (negative). [:issue:`34`]

sora.observer
^^^^^^^^^^^^^

- Function Observer.altaz() that calculates the altitude and azimuth for a given target 
  and instant. [:issue:`34`]

sora.prediction
^^^^^^^^^^^^^^^

- Four new parameters were added to `plot_occ_map()`: `path`: for the user to select
  a directory where to save the plots; `site_name`: If True, the name of the sites
  will be plotted; `chord_delta` and `chord_geo`: for the user to plot the path of
  a chord from distance of the center or passing by some coordinate, respectively. [:issue:`35`]

- Two methods were added to `PredictionTable()` to help the user to remove bad events
  from table: `keep_from_selected_images()` and `remove_occ()`. [:issue:`35`]


API Changes
-----------

sora.config
^^^^^^^^^^^

- config module is now a directory. It now includes a module with decorators
  and another for variables. [:issue:`31, 35`]

sora.ephem
^^^^^^^^^^

- In EphemKernel, `code` argument was replaced by `spkid`. When using 'code',
  a FutureWarning is raised stating `code` as deprecated and will be removed from v1.0. [:issue:`26`]

sora.lightcurve
^^^^^^^^^^^^^^^

- In LightCurve.immersion and LightCurve.emersion, an error will rise when these values were not 
  instanciated or fitted. [:issue:`34`]

- Now the user has the possibility to redefine `tref`, `immersion`, `emersion`,
  `initial_time` and `end_time` after instantiated. [:issue:`35`]

- `lambda_0` argument was replaced by `central_bandpass` and `delta_lambda` by `delta_bandpass`. 
  When using 'lambda_0' or `delta_lambda`, a FutureWarning is raised stating `lambda_0` or `delta_lambda`
  as deprecated and will be removed from v1.0. [:issue:`36`]

sora.occultation
^^^^^^^^^^^^^^^^

- Occultation.new_astrometric_positions() now shows a warning when time is far
  by more than 1 day from the occultation closest approach. [:issue:`21`]

- Occultation.to_log() and print(Occultation) added the polar radius, equivalent radius, 
  the Sun-Geocenter-Target angle and the Moon-Geocenter-Target angle, geocentric albedo,
  the altitude and azimuth of the target for each Observer. [:issue:`17`]

- In `fit_ellipse()`, `pos_angle` and `dpos_angle` were deprecated in favor of
  `position_angle` and `dposition_angle`. [:issue:`35`]

- Changed "GCRS" to "Geocentric" in the string representation to avoid confusion
  about the reference frame. [:issue:`35`]
  
sora.prediction
^^^^^^^^^^^^^^^

- prediction() now calculates the ephemeris inside each division to avoid memory overflow. [:issue:`31`]

- PredictionTable.to_ow() will now raise a warning if the radius or the error of
  the ephemeris is not present. [:issue:`35`]

sora.star
^^^^^^^^^^^^^^^

- Now Star downloads all parameters from Gaia and saves them in the `meta_gaia` attribute [:issue:`35`]


Bug Fixes
---------

sora.ephem
^^^^^^^^^^

- Added function get_position() to EphemPlanete. This corrects a bug that prevented
  Occultation to run with EphemPlanete. [:issue:`41`]

- Fixed bug in EphemJPL where `id_type` was redefined inside __init__(). [:issue:`41`]

sora.lightcurve
^^^^^^^^^^^^^^^

- Fixed error that appears when the fit was done separately (immersion and emersion times). 
  Now the final model agrees with the fitted values.   [:issue:`9`]

- Fixed error when the file with the light curve has three columns. [:issue:`19`]

- Fixed error when the exptime within the LightCurve was set as zero or negative. [:issue:`23`]

- Fixed error in the automatic mode of LightCurve.normalize(). [:issue:`34`]

- Fixed bug that was raised in LightCurve.log() when there were no initial or end times
  for lightcurves instantiated with immersion and emersion. [:issue:`35`]

sora.occultation
^^^^^^^^^^^^^^^^

- Corrected error calculation using err = sqrt(star_err^2 + fit_err^2) [:issue:`18`]

- Occultation.plot_occ_map() now uses the fitted ellipse to calculate the projected shadow radius [:issue:`22`]

- Corrected bug that raised an error when calling Occultation.get_map_sites()
  and there were no observation added to Occultation. [:issue:`31`]

- Corrected bug that did not save the fitted params in all occultations when
  more than one occultation was used in fit_ellipse(). [:issue:`35`]

- Added `axis_labels` and `lw` (linewidth) to Occultation.plot_chords(). [:issue:`35`]

sora.prediction
^^^^^^^^^^^^^^^

- Fixed error that was generated when only one prediction was found. [:issue:`16`]

- Fixed error in the output format of PredictionTable.to_ow() when coordinate was positive [:issue:`35`]


SORA v0.1 (2020/May/20)
=======================

Classes
-------

The documentation of all classes and functions are on their docstrings,
while the scientific part is presented in the full documentation.
Here follows a list with the main Classes:

**Ephem** Three Classes created to generate geocentric ephemeris for a given solar system object.
**EphemJPL** queries the JPL Horizons service and download ephemeris information.
**EphemKernel** reads the BSP files to calculate the ephemeris using the Spiceypy package.
**EphemPlanet** reads an ASCII file with previously determined positions and interpolate them for a given instant.

JPL Horizons - https://ssd.jpl.nasa.gov/horizons.cgi

**Star** Class created to deal with the star parameters. From the Gaia-DR2 Source ID
or a sky region, it queries the VizieR service and downloads the star’s information.
From Gaia DR2 Catalog (Gaia Collaboration 2016a, 2016b and 2018) it gets the RA, DEC,
parallax, proper motions, G magnitude and star radius; from the NOMAD Catalog
(Zacharias et al. 2004) it gets the B, V, R, J, H and K magnitudes.
The user can calculate the ICRS coordinate of the star at any epoch.
It can be barycentric (corrected from proper motion) or geocentric (corrected
from proper motion and parallax). Also, the apparent diameter of the star is calculated
using Gaia DR2 information, or some models such as Van Belle (1999) and  Kervella et al. (2004).

Gaia - Gaia Collaboration 2016a, 2016b and 2018
Mission: https://ui.adsabs.harvard.edu/abs/2016A\%26A...595A...1G/abstract
DR1: https://ui.adsabs.harvard.edu/abs/2016A\%26A...595A...2G/abstract
DR2: https://ui.adsabs.harvard.edu/abs/2018A\%26A...616A...1G/abstract
VizieR - https://vizier.u-strasbg.fr/viz-bin/VizieR
NOMAD - Zacharias et al. 2004 https://ui.adsabs.harvard.edu/abs/2004AAS...205.4815Z/abstract
Van Belle, 1999 - https://ui.adsabs.harvard.edu/abs/1999PASP..111.1515V/abstract
Kervella, 2004 - https://ui.adsabs.harvard.edu/abs/2004A%26A...426..297K/abstract

**Observer**: Object Class created to deal with the observer location. The user can
also download the ground-based observatories from the Minor Planet Center (MPC) database.

MPC sites - https://minorplanetcenter.net/iau/lists/ObsCodesF.html

**Light Curve**: Object Class that receives the observational light curve (with time
and the occulted star normalized photometry relative to reference stars) and some
observational parameters (filter and exposure time). It has functions to determine
the instants that the solar system object enters in front of the star and leaves,
(immersion and emersion times, respectively). The model considers a sharp-edge
occultation model (geometric) convolved with Fresnel diffraction, stellar diameter
(projected at the body distance) and finite integration time (Widemann et al.,
2009; Sicardy et al., 2011</font>).

Widemann et al. 2009 -  https://ui.adsabs.harvard.edu/abs/2009Icar..199..458W/abstract
Sicardy et al. 2011 -  https://ui.adsabs.harvard.edu/abs/2011Natur.478..493S/abstract

**Occultation**: Main Object Class within SORA, created to analyze stellar
occultations, and control all the other Object Classes within this package.
Its functions allow converting the times for each observatory in the occulted
body positions in the sky plane relative to the occulted star (f, g) (IERS
Conventions). Also, to obtain the best ellipse parameters (centre position,
apparent equatorial radius, oblateness and the position angle of the apparent
polar radius) that fit the points. The results are the apparent size, shape and
astrometrical position of the occulting body.

IERS Conventions: https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html

Some extra Objects Classes:

**PredictionTable**: Using the **prediction** function within SORA results in an
Object Class that is a slight modification of an AstropyTable. The added changes
allow to create the occultation map for each prediction, convert into specific
formats, such as OccultWatcher and PRAIA (Assafin et al. (2011)).

OccultWatcher - https://www.occultwatcher.net/
Assafin et al., 2011 - https://ui.adsabs.harvard.edu/abs/2011gfun.conf...85A/abstract

**ChiSquare**: This Object Class is the result of the fitting functions within
SORA, such as _LightCurve.occ_lcfit()_ and _Occultation.fit_ellipse()_.
This Class has functions that allow viewing the values that minimize the :math:`{\chi^2}`
tests, the uncertainties within :math:`{n\sigma}`, plotting the tests, and saving the values.


INPUTS AND OUTPUTS
------------------

INPUTS
^^^^^^
- **Event Related (Star and Ephem)**
 
  - Object Name or provisory designation
  - Object Code (only for EphemKernel)
  - BSP file and name (only for EphemKernel)
  - DE file and name (only for EphemKernel)
  - Ephemeris offset for RA and DEC - :math:`{\Delta \alpha \cdot \cos \delta}`, :math:`{\Delta \delta}` (set as 0,0)
  - Occultation date and time
  - Occulted star coordinates RA and DEC; or Gaia code
  - Star offset for RA and DEC - :math:`{\Delta \alpha \cdot \cos \delta}`, :math:`{\Delta \delta}` (set as 0,0)

- **Observer Related**
 
  - Site name and location (latitude, longitude, and height; or IAU/MPC code)
  - Light curve file and name; or array with fluxes and times; or immersion and emersion times
  - Exposure time in seconds
  - Observational bandwidth in microns (set as 0.7 :math:`{\pm}` 0.3 microns, Clear)

- **Fitting Related**

  - Initial guess for light curve fitting: immersion, emersion and opacity.
  - Range to explore all three parameters
  - Initial guess for ellipse parameters: center (f,g), equatorial radius, oblateness, and position angle
  - Range to explore all five parameters


OUTPUTS
^^^^^^^

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
  - Immersion Time and uncertainty - :math:`{1\sigma}` and :math:`{3\sigma}` (s)
  - Emersion Time and uncertainty (yyyy-mm-dd hh:mm:ss.sss +/- s.sss)
  - :math:`{\chi^2}` fit model
  - Emersion Time and uncertainty - :math:`{1\sigma}` and :math:`{3\sigma}` (s)
  - Minimum Chi-square - :math:`{\chi^2_{min}}`
  - Number of fitted points for im- and emersion
  - Number of fitted parameters
  - Minimum Chi-square per degree of freedom - :math:`{\chi^2_{min-pdf}}`

- Elipse fit procedure

  - Fitted parameters: Equatorial radius and uncertainty (km); Center position (:math:`{f_0}`, :math:`{g_0}`) and :math:`{1\sigma}` uncertainties (km, km); Oblateness and uncertainty; Position angle and uncertainty (degree)
  - Minimum Chi-square -  :math:`{\chi_{min}^2}`
  - Minimum Chi-square per degree of freedom - :math:`{\chi_{min-pdf}^2}`
  - Number points used to fit ( X points from Y chords )
  - Astrometric object center position at occ. time and uncertainty (hh mm ss.sss +dd mm ss.sss :math:`{\pm}` mas)

- Plots and files (some are optional)

  - Prediction map (Lucky Star model)
  - Normalised light curve - for each site (x = time; y = flux)
  - Chi-square map for immersion and emersion times (x = time; y = :math:`{\chi^2}`)
  - Light curve and synthetic LC- for each site (x = time; y = flux)
  - Chords projected in sky plane (x = :math:`{\xi}` (km); y = :math:`{\eta}` (km) )
  - Chi-square map for each ellipse parameter (x = time; y = :math:`{\chi^2_{param}}`)
  - Chords projected in sky plane and the best ellipse fitted with :math:`{1\sigma}` uncertainties (x = :math:`{\xi}` (km); y = :math:`{\eta}` (km))
  - Log file with all information

