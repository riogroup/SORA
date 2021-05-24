import astropy.units as u
import numpy as np

__all__ = ['filter_negative_chord', 'positionv']


def positionv(star, ephem, observer, time):
    """Calculates the position and velocity of the occultation shadow relative
    to the observer.

    Parameters
    ----------
    star : `sora.Star`
        The coordinate of the star in the same reference frame as the ephemeris.
        It must be a Star object.

    ephem : `sora.Ephem`
        The object ephemeris. It must be an Ephemeris object.

    observer : `sora.Observer`
        The Observer information. It must be an Observer object.

    time : `astropy.time.Time`
        Reference instant to calculate position and velocity. It can be a string
        in the ISO format (yyyy-mm-dd hh:mm:ss.s) or an astropy Time object.

    Returns
    -------
    f, g, vf, vg : `list`
        The orthographic projection of the shadow relative to the observer.
        ``'f'`` is in the x-axis (East-West direction; East positive).
        ``'g'`` is in the y-axis (North-South direction; North positive).
    """
    from sora.ephem import EphemPlanete, EphemJPL, EphemKernel, EphemHorizons
    from sora.observer import Observer
    from sora.star import Star
    from astropy.time import Time

    if type(star) != Star:
        raise ValueError('star must be a Star object')
    if type(ephem) not in [EphemPlanete, EphemJPL, EphemKernel, EphemHorizons]:
        raise ValueError('ephem must be an Ephemeris object')
    if type(observer) != Observer:
        raise ValueError('observer must be an Observer object')
    time = Time(time)

    coord = star.geocentric(time)
    dt = 0.1*u.s

    if type(ephem) == EphemPlanete:
        ephem.fit_d2_ksi_eta(coord, verbose=False)
    ksio1, etao1 = observer.get_ksi_eta(time=time, star=coord)
    ksie1, etae1 = ephem.get_ksi_eta(time=time, star=coord)

    f = ksio1-ksie1
    g = etao1-etae1

    ksio2, etao2 = observer.get_ksi_eta(time=time+dt, star=coord)
    ksie2, etae2 = ephem.get_ksi_eta(time=time+dt, star=coord)

    nf = ksio2-ksie2
    ng = etao2-etae2

    vf = (nf-f)/0.1
    vg = (ng-g)/0.1

    return f, g, vf, vg


def filter_negative_chord(chord, chisquare, step=1, sigma=0):
    """Get points for the ellipse with the given input parameters.

    Parameters
    ----------
    chord : `sora.observer.Chord`
        Chord object, must be associated to an Occultation to work.

    chisquare : `sora.extra.ChiSquare`
        Resulted ChiSquare object of fit_ellipse.

    sigma : `int`, `float`
        Uncertainty of the expected ellipse, in km.

    step : `int`, `float`, `str`
        If a number, it corresponds to the step, in seconds, for each point of
        the chord path. The step can also be equal to ``'exposure'``. In this
        case, the chord path will consider the lightcurve individual times and
        exptime.
    """
    from sora.config.visuals import progressbar
    from sora.extra import ChiSquare
    from sora.extra.utils import get_ellipse_points

    keep = []
    if step == 'exposure':
        try:
            step = chord.lightcurve.exptime/10.
        except AttributeError:
            raise AttributeError('Chord.lightcurve does not have "exptime"')
        time_all = np.arange(chord.lightcurve.time.min(), chord.lightcurve.time.max(), step)
        time_exposure = np.array([])
        for i in range(len(chord.lightcurve.time)):
            event_model = (time_all > chord.lightcurve.time[i] - chord.lightcurve.exptime/2.) & (
                        time_all < chord.lightcurve.time[i] + chord.lightcurve.exptime/2.)
            time_exposure = np.append(time_exposure, time_all[event_model])
            f_all, g_all = chord.get_fg(time=time_exposure*u.s + chord.lightcurve.tref)
    else:
        f_all, g_all = chord.path(segment='full', step=step)
    for i in progressbar(range(len(chisquare.data['chi2'])), 'Filter chord: {}'.format(chord.name)):
        df_all = (f_all - chisquare.data['center_f'][i])
        dg_all = (g_all - chisquare.data['center_g'][i])

        r_all = np.sqrt(df_all**2 + dg_all**2)
        cut = r_all < 1.5*chisquare.data['equatorial_radius'][i]
        df_path = df_all[cut]
        dg_path = dg_all[cut]
        r_path = r_all[cut]
        theta_path = np.arctan2(dg_path, df_path)

        r_ellipse = get_ellipse_points(theta_path,
                                       equatorial_radius=chisquare.data['equatorial_radius'][i],
                                       oblateness=chisquare.data['oblateness'][i],
                                       center_f=chisquare.data['center_f'][i],
                                       center_g=chisquare.data['center_g'][i],
                                       position_angle=chisquare.data['position_angle'][i])[2]
        keep.append(np.all(r_path - r_ellipse + sigma > 0))

    filtered_chisquare = ChiSquare(chisquare.data['chi2'][keep], chisquare.npts,
                                   center_f=chisquare.data['center_f'][keep],
                                   center_g=chisquare.data['center_g'][keep],
                                   equatorial_radius=chisquare.data['equatorial_radius'][keep],
                                   oblateness=chisquare.data['oblateness'][keep],
                                   position_angle=chisquare.data['position_angle'][keep])
    return filtered_chisquare
