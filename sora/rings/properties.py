import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from .utils import calc_coef_projecao, project_to_ring_plane
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



def sky_width(fi, gi, fe, ge, center_f=0, center_g=0):
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


def squarewell_properties(model, chord, sigma=1, center_f=0, center_g=0, ring=None, pole_orientation=None, particle_size=0.01):
    # inserir kwargs immersion, emersion, opacity, outros?
    # o que acontece se imm, eme e opa não forem arrays?
    
    properties = {}
    
    label = model.fit_meta['label']    
    chisquare = model.lightcurve.chi2_maps[label]
    samples = chisquare.get_values(sigma=sigma)
    
    imm = np.array(samples['immersion'], ndmin=1)
    eme = np.array(samples['emersion'], ndmin=1)
    opa = np.array(samples['opacity'], ndmin=1)       # p = 1 - sqrt(T)
    
    transmittance = (1 - opa)**2
    apparent_opacity = 1 - transmittance
    
    airy_scale = (chord.lightcurve.central_bandpass * u.micrometer.to('km') / (2 * particle_size)) * chord.lightcurve.dist * u.au.to('km')  # assumed particle size = 0.1 km
    properties['airy_scale'] = airy_scale
    
    tref = Time(chord.lightcurve.tref)
    ti = tref + imm * u.s
    te = tref + eme * u.s
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

    properties['transmittance'] = np.median(transmittance)
    properties['dtransmittance'] = (transmittance.max() - transmittance.min()) / 2
    properties['apparent_opacity'] = np.median(apparent_opacity)
    properties['dapparent_opacity'] = (apparent_opacity.max() - apparent_opacity.min()) / 2
    
    
    if airy_scale > skywidth:
        opacity = opa
    else:
        opacity = apparent_opacity
    
    properties['opacity'] = np.median(opacity)
    properties['dopacity'] = (opacity.max() - opacity.min()) / 2 
    
    pole = pole_orientation if pole_orientation is not None else (getattr(ring, "pole_orientation", None) if ring is not None else None)
    
    if pole is not None:
        pole_orientation = SkyCoord(pole).icrs
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
                opacity=opa,
                pole_orientation=pole_orientation.icrs,
                ephem=ephem,
                tref=tref)

            normal_opt_depth, dnormal_opt_depth = corrected_normal_optical_depth(
                opacity=opa,
                pole_orientation=pole_orientation.icrs,
                ephem=ephem,
                tref=tref)

            properties['normal_opacity'], properties['dnormal_opacity'] = n_opacity, dn_opacity
            properties['normal_optical_depth'], properties['dnormal_optical_depth'] = normal_opt_depth, dnormal_opt_depth

        else:
            n_opacity, dn_opacity = normal_opacity(
                opacity=opa,
                pole_orientation=pole_orientation.icrs,
                ephem=ephem,
                tref=tref)

            normal_opt_depth, dnormal_opt_depth = normal_optical_depth(
                opacity=opa,
                pole_orientation=pole_orientation.icrs,
                ephem=ephem,
                tref=tref)

            properties['normal_opacity'], properties['dnormal_opacity'] = n_opacity, dn_opacity
            properties['normal_optical_depth'], properties['dnormal_optical_depth'] = normal_opt_depth, dnormal_opt_depth

    return properties

def _signed_radius(r):
    r = np.asarray(r, dtype=float)
    i0 = int(np.argmin(r))
    sign = np.ones_like(r)
    sign[:i0] = -1.0
    return sign * r, sign

def _center_fwhm_km(r, y, sign=None):
    r = np.asarray(r, dtype=float)
    y = np.asarray(y, dtype=float)
    s = r if sign is None else np.asarray(sign, dtype=float) * r

    idx = np.argsort(s)
    s = s[idx]
    y = y[idx]

    if len(s) < 3 or not np.isfinite(y).any():
        return np.nan, np.nan

    i0 = int(np.nanargmax(y))
    y0 = y[i0]
    if not np.isfinite(y0) or y0 <= 0:
        return np.nan, np.nan

    half = 0.5 * y0

    left = np.where(y[:i0] <= half)[0]
    right = np.where(y[i0:] <= half)[0] + i0
    if len(left) == 0 or len(right) == 0:
        return s[i0], np.nan

    iL = left[-1]
    iR = right[0]

    def interp_x(i_a, i_b):
        x1, x2 = s[i_a], s[i_b]
        y1, y2 = y[i_a], y[i_b]
        if y2 == y1:
            return 0.5 * (x1 + x2)
        return x1 + (half - y1) * (x2 - x1) / (y2 - y1)

    xL = interp_x(iL, iL + 1)
    xR = interp_x(iR - 1, iR)

    return s[i0], (xR - xL)

def lorentzian_properties(*, model, chord, ring=None, center_f=0, center_g=0):
    if chord is None:
        raise ValueError("Chord must be provided.")

    lc = chord.lightcurve
    tref = Time(lc.tref)
    t_sec = np.asarray(lc.time, dtype=float)

    if t_sec.size < 3:
        raise ValueError("Not enough samples in the selected time window.")

    t_time = tref + t_sec * u.s

    F = np.asarray(model.compute(time=t_sec), dtype=float)
    tau_app = -np.log(F)
    ew_app_integrand = 1.0 - F

    f, g = chord.get_fg(time=t_time)  
    r_sky = np.sqrt((f - center_f)**2 + (g - center_g)**2)

    rsky_signed, sign_sky = _signed_radius(r_sky)

    idx = np.argsort(rsky_signed)
    x = rsky_signed[idx]
    tau_x = tau_app[idx]
    ew_x = ew_app_integrand[idx]

    out = {}
    out["apparent_equivalent_width_model"] = float(np.trapz(ew_x, x))
    out["apparent_optical_depth_integral_model"] = float(np.trapz(tau_x, x))

    i_peak = int(np.argmax(tau_x))
    out["peak_distance_sky_model"] = float(x[i_peak])
    c_sky, w_sky = _center_fwhm_km(np.abs(x), tau_x, sign=np.sign(x))
    out["fwhm_sky_km_model"] = float(w_sky)

    if ring is not None:
        time_evt = chord._shared_with["chordlist"]["time"]
        _, B = ring.get_ring_orientation(time=time_evt)
        sinB = abs(np.sin(B).value)

        x_ring, y_ring = ring.to_ring_plane(
            f=f, g=g,
            time=time_evt,
            center_f=center_f,
            center_g=center_g
        )
        r_ring = np.sqrt(x_ring**2 + y_ring**2)
        rring_signed, sign_ring = _signed_radius(r_ring)

        idx = np.argsort(rring_signed)
        xr = rring_signed[idx]
        tauN = tau_app[idx] * 0.5 * sinB
        pN = ew_app_integrand[idx] * 0.5 * sinB

        out["normal_equivalent_width_model"] = float(np.trapz(pN, xr))
        out["normal_optical_depth_integral_model"] = float(np.trapz(tauN, xr))

        i_peak_r = int(np.argmax(tauN))
        out["peak_distance_ring_model"] = float(xr[i_peak_r])
        c_ring, w_ring = _center_fwhm_km(np.abs(xr), tauN, sign=np.sign(xr))
        out["fwhm_ring_km_model"] = float(w_ring)

    return out