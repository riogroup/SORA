import astropy.units as u
import numpy as np
from astropy.time import Time

from sora.config.decorators import deprecated_alias

__all__ = ['fit_ellipse']


@deprecated_alias(pos_angle='position_angle', dpos_angle='dposition_angle', log='verbose')  # remove this line for v1.0
def fit_ellipse(*args, equatorial_radius, dequatorial_radius=0, center_f=0, dcenter_f=0, center_g=0,
                dcenter_g=0, oblateness=0, doblateness=0, position_angle=0, dposition_angle=0,
                loop=10000000, number_chi=10000, dchi_min=None, verbose=False, ellipse_error=0, sigma_result=1):
    """Fits an ellipse to given occultation using given parameters.

    Parameters
    ----------
    center_f : `int`, `float`, default=0
        The coordinate in f of the ellipse center.

    center_g : `int`, `float`, default=0
        The coordinate in g of the ellipse center.

    equatorial_radius : `int`, `float`
        The Equatorial radius (semi-major axis) of the ellipse.

    oblateness : `int`, `float`, default=0
        The oblateness of the ellipse.

    position_angle : `int`, `float`, default=0
        The pole position angle of the ellipse in degrees.
        Zero is in the North direction ('g-positive'). Positive clockwise.

    dcenter_f : `int`, `float`
        Interval for coordinate f of the ellipse center.

    dcenter_g : `int`, `float`
        Interval for coordinate g of the ellipse center.

    dequatorial_radius `int`, `float`
        Interval for the Equatorial radius (semi-major axis) of the ellipse.

    doblateness : `int`, `float`
        Interval for the oblateness of the ellipse

    dposition_angle : `int`, `float`
        Interval for the pole position angle of the ellipse in degrees.

    loop : `int`, default=10000000
        The number of ellipses to attempt fitting.

    dchi_min : `int`, `float`
        If given, it will only save ellipsis which chi square are smaller than
        chi_min + dchi_min.

    number_chi : `int`, default=10000
        If dchi_min is given, the procedure is repeated until number_chi is
        reached.

    verbose : `bool`, default=False
        If True, it prints information while fitting.

    ellipse_error : `int`, `float`
        Model uncertainty to be considered in the fit, in km.

    sigma_result : `int`, `float`
        Sigma value to be considered as result.

    Returns
    -------
    chisquare : `sora.ChiSquare`
        A ChiSquare object with all parameters.

    Important
    ---------
    Each occultation is added as the first argument(s) directly.

    Mandatory input parameters: 'center_f', 'center_g', 'equatorial_radius',
    'oblateness', and 'position_angle'.

    Parameters fitting interval: 'dcenter_f', 'dcenter_g', 'dequatorial_radius',
    'doblateness', and 'dposition_angle'. Default values are set to zero.
    Search done between (value - dvalue) and (value + dvalue).


    Examples
    --------
    To fit the ellipse to the chords of occ1 Occultation object:

    >>> fit_ellipse(occ1, **kwargs)

    To fit the ellipse to the chords of occ1 and occ2 Occultation objects together:

    >>> fit_ellipse(occ1, occ2, **kwargs)
    """
    from sora.extra import ChiSquare
    from sora.config.visuals import progressbar_show
    from astropy.coordinates import Angle
    from .core import Occultation

    v = {'dcenter_f': dcenter_f, 'dcenter_g': dcenter_g, 'doblateness': doblateness, 'dposition_angle': dposition_angle,
         'dequatorial_radius': dequatorial_radius, 'ellipse_error': ellipse_error, 'sigma_result': sigma_result,
         'dchi_min': dchi_min}
    for key, item in v.items():
        if item is not None and item < 0:
            raise ValueError("{} must be a positive number.".format(key))

    values = []
    chord_name = []
    if len(args) == 0:
        raise ValueError('No occultation have been given as input.')
    for occ in args:
        if not isinstance(occ, Occultation):
            raise TypeError('Given argument must be an Occultation object.')
        for name, chord in occ.chords.items():
            if chord.status() == 'positive':
                if chord.is_able['immersion']:
                    f, g, vf, vg = chord.get_fg(time='immersion', vel=True)
                    err = np.linalg.norm([vf, vg])*chord.lightcurve.immersion_err
                    values.append([f, g, err])
                    chord_name.append(name + '_immersion')
                if chord.is_able['emersion']:
                    f, g, vf, vg = chord.get_fg(time='emersion', vel=True)
                    err = np.linalg.norm([vf, vg])*chord.lightcurve.emersion_err
                    values.append([f, g, err])
                    chord_name.append(name + '_emersion')

    controle_f0 = Time.now()
    f0_chi = np.array([])
    g0_chi = np.array([])
    a_chi = np.array([])
    obla_chi = np.array([])
    posang_chi = np.array([])
    chi2_best = np.array([])

    while len(f0_chi) < number_chi:
        progressbar_show(len(f0_chi), number_chi, prefix='Ellipse fit:')
        chi2 = np.zeros(loop)
        f0 = center_f + dcenter_f*(2*np.random.random(loop) - 1)
        g0 = center_g + dcenter_g*(2*np.random.random(loop) - 1)
        a = equatorial_radius + dequatorial_radius*(2*np.random.random(loop) - 1)
        obla = oblateness + doblateness*(2*np.random.random(loop) - 1)
        obla[obla < 0], obla[obla > 1] = 0, 1
        phi_deg = position_angle + dposition_angle*(2*np.random.random(loop) - 1)
        controle_f1 = Time.now()

        for fi, gi, si in values:
            b = a - a*obla
            phi = phi_deg*(np.pi/180.0)
            dfi = fi-f0
            dgi = gi-g0
            theta = np.arctan2(dgi, dfi)
            ang = theta+phi
            r_model = (a*b)/np.sqrt((a*np.sin(ang))**2 + (b*np.cos(ang))**2)
            f_model = f0 + r_model*np.cos(theta)
            g_model = g0 + r_model*np.sin(theta)
            chi2 += ((fi - f_model)**2 + (gi - g_model)**2)/(si**2 + ellipse_error**2)

        controle_f2 = Time.now()
        if dchi_min is not None:
            region = np.where(chi2 < chi2.min() + dchi_min)[0]
        else:
            region = np.arange(len(chi2))
        chi2_best = np.append(chi2_best, chi2[region])
        if verbose:
            print('Elapsed time: {:.3f} seconds.'.format((controle_f2 - controle_f1).sec))
            print(len(chi2[region]), len(chi2_best))
        f0_chi = np.append(f0_chi, f0[region])
        g0_chi = np.append(g0_chi, g0[region])
        a_chi = np.append(a_chi, a[region])
        obla_chi = np.append(obla_chi, obla[region])
        posang_chi = np.append(posang_chi, phi_deg[region])

    progressbar_show(number_chi, number_chi, prefix='Ellipse fit:')
    chisquare = ChiSquare(chi2_best, len(values), center_f=f0_chi, center_g=g0_chi, equatorial_radius=a_chi,
                          oblateness=obla_chi, position_angle=posang_chi)
    controle_f4 = Time.now()
    if verbose:
        print('Total elapsed time: {:.3f} seconds.'.format((controle_f4 - controle_f0).sec))

    result_sigma = chisquare.get_nsigma(sigma=sigma_result)
    a = result_sigma['equatorial_radius'][0]
    f0 = result_sigma['center_f'][0]
    g0 = result_sigma['center_g'][0]
    obla = result_sigma['oblateness'][0]
    phi_deg = result_sigma['position_angle'][0]
    radial_dispersion = np.array([])
    error_bar = np.array([])
    position_angle_point = np.array([])

    for fi, gi, si in values:
        b = a - a*obla
        phi = phi_deg*(np.pi/180.0)
        dfi = fi-f0
        dgi = gi-g0
        r = np.sqrt(dfi**2 + dgi**2)
        theta = np.arctan2(dgi, dfi)
        ang = theta+phi
        r_model = (a*b)/np.sqrt((a*np.sin(ang))**2 + (b*np.cos(ang))**2)
        radial_dispersion = np.append(radial_dispersion, r - r_model)
        error_bar = np.append(error_bar, si)
        position_angle_point = np.append(position_angle_point, Angle(90*u.deg - theta*u.rad).wrap_at(360 * u.deg).degree)

    for occ in args:
        if isinstance(occ, Occultation):
            occ.fitted_params = {i: result_sigma[i] for i in ['equatorial_radius', 'center_f', 'center_g',
                                                              'oblateness', 'position_angle']}
            occ.chi2_params = {'chord_name': chord_name, 'radial_dispersion': radial_dispersion,
                               'position_angle': position_angle_point, 'radial_error': error_bar,
                               'chi2_min': chisquare.get_nsigma(sigma=sigma_result)['chi2_min'],
                               'nparam': chisquare.nparam, 'npts': chisquare.npts}
    return chisquare
