import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from sora.body.ring.utils import calc_coef_projecao, project_to_ring_plane
from astropy.time import Time


__all__ = ['compute_local_properties']

# Substituir em cada função o polo por "position_angle" ficará melhor


def sky_distance(f, g, center_f=0, center_g=0):
    """Compute the radial distance of a feature in the sky plane relative to a center.

    Parameters
    ----------
    f, g : float or array-like
        Coordinates in the sky plane [km].

    center_f : float
        f-coordinate of the ring center [km].

    center_g : float
        g-coordinate of the ring center [km].

    Returns
    -------
    distance : float
        Mean radial distance, in km.

    ddistance : float
        Uncertainty of the radial distance, in km.
    """

    r = np.sqrt((np.array(f) - center_f)**2 + (np.array(g) - center_g)**2)

    if np.ndim(r) > 0 and len(np.atleast_1d(r)) > 1:
        distance = (r.max() + r.min()) / 2
        ddistance = (r.max() - r.min()) / 2

    else:
        distance = float(r)
        ddistance = 0.0

    return distance, ddistance



def sky_width(fi, gi, fe, ge, center_f=0, center_g=0, dfi=None, dgi=None, dfe=None, dge=None):
    """Compute the width of a feature in the sky plane between immersion and emersion points.

    Parameters
    ----------
    fi, gi : float or array-like
        Immersion point coordinates in the sky plane [km].

    fe, ge : float or array-like
        Emersion point coordinates in the sky plane [km].

    center_f : float
        f-coordinate of the ring center [km].

    center_g : float
        g-coordinate of the ring center [km].

    dfi, dgi, dfe, dge : float, optional
        Uncertainties in fi, gi, fe, ge [km].
        Only used if coordinates are scalars.

    Returns
    -------
    width : float
        Mean width of the feature in km.

    dwidth : float
        Half-width (uncertainty) of the width in km.
    """
    ri = np.sqrt((np.array(fi) - center_f)**2 + (np.array(gi) - center_g)**2)
    re = np.sqrt((np.array(fe) - center_f)**2 + (np.array(ge) - center_g)**2)

    w = np.abs(ri - re)

    if np.ndim(w) > 0 and len(np.atleast_1d(w)) > 1:
        width = (w.max() + w.min()) / 2
        dwidth = (w.max() - w.min()) / 2

    else:
        width = float(w)
        dwidth = 0.0

    return width, dwidth


def equatorial_distance(f, g, pole_orientation, ephem, tref, center_f=0, center_g=0):
    """Compute the radial distance projected onto the ring's equatorial plane.

    Parameters
    ----------
    f, g : float or array-like
        Coordinates in the sky plane [km].

    pole_orientation : astropy.coordinates.SkyCoord
        Ring pole orientation in ICRS.

    ephem : object
        Ephemeris object to get body position.

    center_f : float
        f-coordinate of the ring center [km].

    center_g : float
        g-coordinate of the ring center [km].

    tref : astropy.time.Time
        Reference time of the occultation.

    Returns
    -------
    equatorial_distance : float
        Mean radial distance in the ring plane [km].

    dequatorial_distance : float
        Half-width (uncertainty) of the radial distance in the ring plane [km].
    """
    pos = ephem.get_position(tref)
    position_angle = pos.position_angle(pole_orientation.icrs).rad * u.rad
    opening_angle = np.arcsin(-(np.sin(pole_orientation.icrs.dec) * np.sin(pos.dec) +
                               np.cos(pole_orientation.icrs.dec) * np.cos(pos.dec) *
                               np.cos(pole_orientation.icrs.ra - pos.ra)))

    earth_pole = SkyCoord('12h00m00s', '+90d00m00s', frame='icrs')
    coef, coef_polo = calc_coef_projecao(pos, pole_orientation, opening_angle.value,
                                        position_angle.value, earth_pole)

    x, y = project_to_ring_plane(f, g, coef, coef_polo, ksi_0=center_f, eta_0=center_g)

    r = np.sqrt(x**2 + y**2)
    equatorial_distance = (r.max() + r.min()) / 2
    dequatorial_distance = (r.max() - r.min()) / 2
    return equatorial_distance, dequatorial_distance


def equatorial_width(fi, gi, fe, ge, pole_orientation, ephem, tref, center_f=0, center_g=0):
    """Compute the width projected onto the ring's equatorial plane between immersion and emersion points.

    Parameters
    ----------
    fi, gi : float or array-like
        Immersion point coordinates in the sky plane [km].

    fe, ge : float or array-like
        Emersion point coordinates in the sky plane [km].

    pole_orientation : astropy.coordinates.SkyCoord
        Ring pole orientation in ICRS.

    ephem : object
        Ephemeris object to get body position.

    center_f : float
        f-coordinate of the ring center [km].

    center_g : float
        g-coordinate of the ring center [km].

    tref : astropy.time.Time
        Reference time of the occultation.

    Returns
    -------
    equatorial_width : float
        Mean width in the ring plane [km].

    dequatorial_width : float
        Half-width (uncertainty) of the width in the ring plane [km].
    """
    pos = ephem.get_position(tref)
    position_angle = pos.position_angle(pole_orientation.icrs).rad * u.rad
    opening_angle = np.arcsin(-(np.sin(pole_orientation.icrs.dec) * np.sin(pos.dec) +
                               np.cos(pole_orientation.icrs.dec) * np.cos(pos.dec) *
                               np.cos(pole_orientation.icrs.ra - pos.ra)))

    earth_pole = SkyCoord('12h00m00s', '+90d00m00s', frame='icrs')
    coef, coef_polo = calc_coef_projecao(pos, pole_orientation, opening_angle.value,
                                        position_angle.value, earth_pole)

    xi, yi = project_to_ring_plane(fi, gi, coef, coef_polo, ksi_0=center_f, eta_0=center_g)
    xe, ye = project_to_ring_plane(fe, ge, coef, coef_polo, ksi_0=center_f, eta_0=center_g)

    wri = np.sqrt(xi**2 + yi**2)
    wre = np.sqrt(xe**2 + ye**2)
    w = abs(wri - wre)
    equatorial_width = (w.max() + w.min()) / 2
    dequatorial_width = (w.max() - w.min()) / 2

    return equatorial_width, dequatorial_width


def corrected_normal_opacity(opacity, pole_orientation, ephem, tref):
    """Compute the diffraction-corrected normal opacity of the ring.

    Parameters
    ----------
    opacity : float or array-like
        Apparent opacity from the light curve.

    pole_orientation : astropy.coordinates.SkyCoord
        Ring pole orientation in ICRS.

    ephem : object
        Ephemeris object to get body position.

    tref : astropy.time.Time
        Reference time of the occultation.

    Returns
    -------
    cnormal_opacity : float
        Mean corrected normal opacity.

    dcnormal_opacity : float
        Uncertainty in the corrected normal opacity.
    """
    pos = ephem.get_position(tref)
    opening_angle = np.arcsin(-(np.sin(pole_orientation.icrs.dec) * np.sin(pos.dec) +
                               np.cos(pole_orientation.icrs.dec) * np.cos(pos.dec) *
                               np.cos(pole_orientation.icrs.ra - pos.ra)))
    cn_opacity = opacity * np.abs(np.sin(opening_angle).value)
    cnormal_opacity = (cn_opacity.max() + cn_opacity.min()) / 2
    dcnormal_opacity = (cn_opacity.max() - cn_opacity.min()) / 2
    return cnormal_opacity, dcnormal_opacity


def normal_opacity(opacity, pole_orientation, ephem, tref):
    """Compute the geometric normal opacity without diffraction correction.

    Parameters
    ----------
    opacity : float or array-like
        Apparent opacity from the light curve.

    pole_orientation : astropy.coordinates.SkyCoord
        Ring pole orientation in ICRS.

    ephem : object
        Ephemeris object to get body position.

    tref : astropy.time.Time
        Reference time of the occultation.

    Returns
    -------
    normal_opacity : float
        Mean geometric normal opacity.

    dnormal_opacity : float
        Uncertainty in the geometric normal opacity.
    """
    pos = ephem.get_position(tref)
    opening_angle = np.arcsin(-(np.sin(pole_orientation.icrs.dec) * np.sin(pos.dec) +
                               np.cos(pole_orientation.icrs.dec) * np.cos(pos.dec) *
                               np.cos(pole_orientation.icrs.ra - pos.ra)))
    n_opacity = (1 - (1 - opacity)**2) * np.abs(np.sin(opening_angle).value)
    normal_opacity = (n_opacity.max() + n_opacity.min()) / 2
    dnormal_opacity = (n_opacity.max() - n_opacity.min()) / 2
    return normal_opacity, dnormal_opacity


def corrected_normal_optical_depth(opacity, pole_orientation, ephem, tref):
    """Compute the normal optical depth corrected for diffraction by particles.

    Parameters
    ----------
    opacity : float or array-like
        Apparent opacity from the light curve.

    pole_orientation : astropy.coordinates.SkyCoord
        Ring pole orientation in ICRS.

    ephem : object
        Ephemeris object to get body position.

    tref : astropy.time.Time
        Reference time of the occultation.

    Returns
    -------
    normal_opt_depth : float
        Mean corrected normal optical depth.

    dnormal_opt_depth : float
        Uncertainty in the corrected normal optical depth.
    """
    pos = ephem.get_position(tref)
    opening_angle = np.arcsin(-(np.sin(pole_orientation.icrs.dec) * np.sin(pos.dec) +
                               np.cos(pole_orientation.icrs.dec) * np.cos(pos.dec) *
                               np.cos(pole_orientation.icrs.ra - pos.ra)))

    app_tau = -np.log((1 - opacity)**2)
    normal_tau = (app_tau / 2) * np.abs(np.sin(opening_angle).value)
    normal_opt_depth = (normal_tau.max() + normal_tau.min()) / 2
    dnormal_opt_depth = (normal_tau.max() - normal_tau.min()) / 2
    return normal_opt_depth, dnormal_opt_depth


def normal_optical_depth(opacity, pole_orientation, ephem, tref):
    """Compute the normal optical depth assuming no diffraction.

    Parameters
    ----------
    opacity : float or array-like
        Apparent opacity from the light curve.

    pole_orientation : astropy.coordinates.SkyCoord
        Ring pole orientation in ICRS.

    ephem : object
        Ephemeris object to get body position.

    tref : astropy.time.Time
        Reference time of the occultation.

    Returns
    -------
    normal_opt_depth : float
        Mean normal optical depth (uncorrected).

    dnormal_opt_depth : float
        Uncertainty in the normal optical depth.
    """
    pos = ephem.get_position(tref)
    opening_angle = np.arcsin(-(np.sin(pole_orientation.icrs.dec) * np.sin(pos.dec) +
                               np.cos(pole_orientation.icrs.dec) * np.cos(pos.dec) *
                               np.cos(pole_orientation.icrs.ra - pos.ra)))

    app_tau = -np.log((1 - opacity)**2)
    normal_tau = app_tau * np.abs(np.sin(opening_angle).value)
    normal_opt_depth = (normal_tau.max() + normal_tau.min()) / 2
    dnormal_opt_depth = (normal_tau.max() - normal_tau.min()) / 2
    return normal_opt_depth, dnormal_opt_depth 


def compute_local_properties(immersion, emersion, opacity, chord, center_f=0, center_g=0, pole_orientation=None):
    """Compute physical properties of a feature detected in an occultation light curve.

    Parameters
    ----------
    immersion : float or array-like
        Immersion time(s) in seconds relative to the chord reference time.

    emersion : float or array-like
        Emersion time(s) in seconds relative to the chord reference time.

    opacity : float or array-like
        Apparent opacity from the light curve model.

    chord : sora.occultation.chords.Chord
        Chord object containing light curve and geometry data.

    center_f : float, optional
        f-coordinate of the ring center [km]. Default is 0.

    center_g : float, optional
        g-coordinate of the ring center [km]. Default is 0.

    pole_orientation : astropy.coordinates.SkyCoord, optional
        Ring pole orientation in ICRS. If provided, computes equatorial-plane properties.

    Returns
    -------
    properties : dict
        Dictionary of computed properties including:
        - central_time, dcentral_time
        - sky_distance, dsky_distance
        - sky_width, dsky_width
        - transmittance, dtransmittance
        - equatorial_distance, dequatorial_distance (if pole_orientation is given)
        - equatorial_width, dequatorial_width (if pole_orientation is given)
        - normal_opacity, dnormal_opacity (if pole_orientation is given)
        - normal_optical_depth, dnormal_optical_depth (if pole_orientation is given)
    """
    properties = {}

    immersion = np.array(immersion, ndmin=1)
    emersion = np.array(emersion, ndmin=1)
    opacity = np.array(opacity, ndmin=1)

    airy_scale = (chord.lightcurve.central_bandpass * u.micrometer.to('km') /
                  (2 * 0.001)) * chord.lightcurve.dist * u.au.to('km')  # assumed particle size = 0.001 km

    tref = Time(chord.lightcurve.tref)

    ti = tref + immersion * u.s
    te = tref + emersion * u.s
    tc = ti + (te - ti) * 0.5
    tc_sec = (tc - tref).sec

    fgs = [chord.get_fg(time=t) for t in [ti, tc, te]]
    (fi, gi), (f, g), (fe, ge) = fgs

    ephem = chord._shared_with['chordlist']['body'].ephem

    properties['central_time'] = tc_sec.mean()
    properties['dcentral_time'] = 0.0 if np.isnan(tc_sec.std(ddof=1)) else tc_sec.std(ddof=1)

    distance, ddistance = sky_distance(f, g, center_f=center_f, center_g=center_g)
    properties['sky_distance'], properties['dsky_distance'] = distance, ddistance

    skywidth, dskywidth = sky_width(fi=fi, gi=gi, fe=fe, ge=ge,
                                    center_f=center_f, center_g=center_g)
    properties['sky_width'], properties['dsky_width'] = skywidth, dskywidth

    properties['transmittance'] = (1 - opacity).mean()
    properties['dtransmittance'] = (1 - opacity).std(ddof=1)

    if pole_orientation:
        pole_orientation = SkyCoord(pole_orientation)
        equat_distance, dequat_distance = equatorial_distance(
            f=f, g=g,
            pole_orientation=pole_orientation.icrs,
            ephem=ephem,
            center_f=center_f,
            center_g=center_g,
            tref=tref)

        equat_width, dequat_width = equatorial_width(
            fi=fi, gi=gi, fe=fe, ge=ge,
            pole_orientation=pole_orientation.icrs,
            ephem=ephem,
            center_f=center_f,
            center_g=center_g,
            tref=tref)

        properties['equatorial_distance'], properties['dequatorial_distance'] = equat_distance, dequat_distance
        properties['equatorial_width'], properties['dequatorial_width'] = equat_width, dequat_width

        if airy_scale > skywidth:
            n_opacity, dn_opacity = corrected_normal_opacity(
                opacity=opacity,
                pole_orientation=pole_orientation.icrs,
                ephem=ephem,
                tref=tref)

            normal_opt_depth, dnormal_opt_depth = corrected_normal_optical_depth(
                opacity=opacity,
                pole_orientation=pole_orientation.icrs,
                ephem=ephem,
                tref=tref)

            properties['normal_opacity'], properties['dnormal_opacity'] = n_opacity, dn_opacity
            properties['normal_optical_depth'], properties['dnormal_optical_depth'] = normal_opt_depth, dnormal_opt_depth

        else:
            n_opacity, dn_opacity = normal_opacity(
                opacity=opacity,
                pole_orientation=pole_orientation.icrs,
                ephem=ephem,
                tref=tref)

            normal_opt_depth, dnormal_opt_depth = normal_optical_depth(
                opacity=opacity,
                pole_orientation=pole_orientation.icrs,
                ephem=ephem,
                tref=tref)

            properties['normal_opacity'], properties['dnormal_opacity'] = n_opacity, dn_opacity
            properties['normal_optical_depth'], properties['dnormal_optical_depth'] = normal_opt_depth, dnormal_opt_depth

    return properties
