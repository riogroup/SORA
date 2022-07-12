import astropy.units as u
import numpy as np
from astropy.coordinates import Angle
from astropy.time import Time
from tqdm import tqdm

from sora.body.shape.limb import limb_radial_residual
from sora.config.decorators import deprecated_alias
from sora.config.visuals import progressbar_show
from sora.extra import ChiSquare

__all__ = ['fit_ellipse', 'fit_to_limb', 'fit_shape']


@deprecated_alias(pos_angle='position_angle', dpos_angle='dposition_angle', log='verbose')  # remove this line for v1.0
def fit_ellipse(*args, equatorial_radius, dequatorial_radius=0, center_f=0, dcenter_f=0, center_g=0,
                dcenter_g=0, oblateness=0, doblateness=0, position_angle=0, dposition_angle=0,
                loop=10000, number_chi=10000, dchi_min=None, verbose=False, ellipse_error=0, sigma_result=1,
                method='least_squares', threads=1):
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
        chi_min + dchi_min. By default `None` when used with `method='chisqr`, and
        `3` for other methods.

    number_chi : `int`, default=10000
        In the `chisqr` method if `dchi_min` is given, the procedure is repeated until
        `number_chi` is reached.
        In other methods it is the number of values (simulations) that should lie within
        the provided `sigma_result`.

    verbose : `bool`, default=False
        If True, it prints information while fitting.

    ellipse_error : `int`, `float`
        Model uncertainty to be considered in the fit, in km.

    sigma_result : `int`, `float`
        Sigma value to be considered as result.

    method : `str`, default=`least_squares`
        Method used to perform the fit. Available methods are:
        `chisqr` : monte carlo computation method used in versions of SORA <= 0.2.1.
        `fastchi` : monte carlo computation method, allows multithreading .
        `least_squares` or `ls`: best fit done used levenberg marquardt convergence algorithm.
        `differential_evolution` or ``de`: best fit done using genetic algorithms.
        All methods return a Chisquare object.

    threads : `int`
        Number of threads/workers used to perform parallel computations of the chi square
        object. It works with all methods except `chisqr`, by default 1.

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
    from astropy.coordinates import Angle
    from .core import Occultation
    from sora.stats import Parameters, least_squares, differential_evolution
    from warnings import warn
    from multiprocessing import Pool

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

    method = method.lower()

    if method not in ['chisqr', 'least_squares', 'ls', 'fastchi', 'differential_evolution', 'de']:
        warn(f'Invalid method `{method}` provided. Setting to default.')
        method = 'chisqr'

    set_bestchi = False # variable used with convergence algorithms and fastchi

    if (method == 'chisqr'):

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
            f0_chi = np.append(f0_chi, f0[region])
            g0_chi = np.append(g0_chi, g0[region])
            a_chi = np.append(a_chi, a[region])
            obla_chi = np.append(obla_chi, obla[region])
            posang_chi = np.append(posang_chi, phi_deg[region])

        progressbar_show(number_chi, number_chi, prefix='Ellipse fit:')
        chisquare = ChiSquare(chi2_best, len(values), center_f=f0_chi, center_g=g0_chi, equatorial_radius=a_chi,
                            oblateness=obla_chi, position_angle=posang_chi)


    else:
        # generate parameters object
        initial = Parameters()
        eqradius_vary = False if dequatorial_radius == 0 else True
        initial.add(name='equatorial_radius', value=equatorial_radius,
                    minval=-np.inf if not eqradius_vary else (equatorial_radius - dequatorial_radius),
                    maxval=np.inf if not eqradius_vary else (equatorial_radius + dequatorial_radius),
                    free=eqradius_vary)

        center_f_vary = False if dcenter_f == 0 else True
        initial.add(name='center_f', value=center_f,
                    minval=-np.inf if not center_f_vary else (center_f-dcenter_f),
                    maxval=np.inf if not center_f_vary else (center_f+dcenter_f), free=center_f_vary)

        center_g_vary = False if dcenter_g == 0 else True
        initial.add(name='center_g', value=center_g,
                    minval=-np.inf if not center_g_vary else (center_g-dcenter_g),
                    maxval=np.inf if not center_g_vary else (center_g+dcenter_g), free=center_g_vary)

        oblateness_vary = False if doblateness == 0 else True
        ob_minval = 0 if ((oblateness - doblateness) < 0) else (oblateness - doblateness)
        ob_maxval = 1 if ((oblateness + doblateness) >1) else (oblateness + doblateness)
        if doblateness >= 0.5:
            ob_minval, ob_maxval = 0, 1
        initial.add(name='oblateness', value=oblateness,
                    minval=-np.inf if not oblateness_vary else ob_minval,
                    maxval=np.inf if not oblateness_vary else ob_maxval, free=oblateness_vary)

        pangle_vary = False if dposition_angle == 0 else True
        initial.add(name='position_angle', value=position_angle,
                    minval=-np.inf if not pangle_vary else (position_angle - dposition_angle),
                    maxval=np.inf if not pangle_vary else (position_angle + dposition_angle),
                    free=pangle_vary)

        fi, gi, si = np.array(values).T


        #run least_squares
        if (method == 'least_squares') or (method == 'ls'):
            result = least_squares(ellipseError, initial, args=(fi, gi, si, ellipse_error), algorithm='trf', sigma=sigma_result)
            # pass the best solution as input to the fastchi method
            equatorial_radius = result.params['equatorial_radius'].value
            center_f = result.params['center_f'].value
            center_g = result.params['center_g'].value
            oblateness = result.params['oblateness'].value
            position_angle = result.params['position_angle'].value
            bestchi, set_bestchi = result.chisqr, True
            method = 'fastchi'

        # run differential_evolution
        if (method == 'differential_evolution') or (method == 'de'):
            result = differential_evolution(ellipseError, initial, args=(fi, gi, si, ellipse_error), sigma=sigma_result)
            # pass the best solution as input to the fastchi method
            equatorial_radius = result.params['equatorial_radius'].value
            center_f = result.params['center_f'].value
            center_g = result.params['center_g'].value
            oblateness = result.params['oblateness'].value
            position_angle = result.params['position_angle'].value
            bestchi, set_bestchi = result.chisqr, True
            method = 'fastchi'

        # run fastchi
        if (method == 'fastchi'):

            if not set_bestchi:
                bestchi = None

            if threads is None:
                threads = 1

            argsloop = [values, bestchi, equatorial_radius, dequatorial_radius, center_f, dcenter_f,
                    center_g, dcenter_g, oblateness, doblateness, position_angle, dposition_angle,
                    loop, np.ceil(number_chi/threads), dchi_min, ellipse_error, False]

            argsloop_verbose = [values, bestchi, equatorial_radius, dequatorial_radius, center_f, dcenter_f,
                            center_g, dcenter_g, oblateness, doblateness, position_angle, dposition_angle,
                            loop, np.ceil(number_chi/threads), dchi_min, ellipse_error, True]

            pool_args = []
            if verbose:
                pool_args.append(argsloop_verbose)
                for i in range(threads-1):
                    pool_args.append(argsloop)
            else:
                pool_args = [argsloop for t in range(threads)]


            with Pool(processes=threads) as pool:
                pool_result = pool.starmap(__fit_ellipse_parallel, pool_args)

            result = [[],[],[],[],[],[]]

            for j in range(6):
                for i in range(threads):
                    for k in pool_result[i][j]:
                        result[j].append(k)

        chisquare = ChiSquare(np.array(result[0]), len(values), center_f=np.array(result[1]), center_g=np.array(result[2]),
                              equatorial_radius=np.array(result[3]), oblateness=np.array(result[4]), position_angle=np.array(result[5]))

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


def ellipse(parameters, x_values, y_values):
    '''
    Returns an ellipse give a set of parameters.

    Parameters
    ----------
    parameters : `stats.Parameters`
        Parameters object that describe the ellipse
    x_values : `np.ndarray`, `list`
        X axis array of values
    y_values : `np.ndarray`, `list`
        Y axis array of values

    Returns
    -------
    [x_model, y_model] : `list`
        Array containing the ellipse data points.
    '''
    p = parameters.valuesdict()

    b = p['equatorial_radius'] - p['equatorial_radius']*p['oblateness']
    phi = p['position_angle'] * np.pi/180.0
    theta = np.arctan2( y_values - p['center_g'], x_values - p['center_f'])
    angle = theta + phi
    radial_model = ( p['equatorial_radius'] * b )/np.sqrt( ( p['equatorial_radius'] * np.sin( angle ) )**2 + ( b * np.cos( angle ) )**2 )
    x_model = p['center_f'] + radial_model*np.cos( theta )
    y_model = p['center_g'] + radial_model*np.sin( theta )

    return [x_model, y_model]


def ellipseError(parameters, f, g, uncertainty, ellipse_error=0):
    '''
    Returns an array of residuals of an ellipse. Depends on the
    ellipse() fuction.

    Parameters
    ----------
    parameters : `stats.Parameters`
        Object containing parameters information
    f : `float`
        f measurements.
    g : `float`
        g measurements.
    uncertainty : `float`
        Uncertainty of the measurements.
    ellipse_error : int, optional
        Ellipse additional uncertainty, by default 0.

    Returns
    -------
    [array] : float
        Array containing the residuals
    '''
    f_model, g_model = ellipse(parameters, f, g)
    return (( (f - f_model)**2 + (g - g_model)**2 )/(uncertainty**2 + ellipse_error**2) )


def __fit_ellipse_parallel(values, bestchi, equatorial_radius, dequatorial_radius, center_f, dcenter_f,
                 center_g, dcenter_g, oblateness, doblateness, position_angle, dposition_angle,
                loop, number, dchi_min, ellipse_error, verbose):

    """Private function that fits an ellipse to given occultation using given parameters.

    Parameters
    ----------
    values : `float`
        Array containing the data to be fitted [f, g, uncertainty]

    bestchi : `bool` or None
        Variable used to allow passing bestfit restults found with
        convergence methods to the chisqr maps.

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

    loop : `int`
        The number of ellipses to attempt fitting.

    dchi_min : `int`, `float`
        If given, it will only save ellipsis which chi square are smaller than
        chi_min + dchi_min. By default `None` when used with `method='chisqr`, and
        `3` for other methods.

    number_chi : `int`, default=10000
        In the `chisqr` method if `dchi_min` is given, the procedure is repeated until
        `number_chi` is reached.
        In other methods it is the number of values (simulations) that should lie within
        the provided `sigma_result`.

    verbose : `bool`, default=False
        If True, it prints information while fitting.

    ellipse_error : `int`, `float`
        Model uncertainty to be considered in the fit, in km.

    sigma_result : `int`, `float`
        Sigma value to be considered as result.
    """

    f0_chi = np.array([]) if bestchi is None else np.array([center_f])
    g0_chi = np.array([]) if bestchi is None else np.array([center_g])
    a_chi = np.array([]) if bestchi is None else np.array([equatorial_radius])
    obla_chi = np.array([]) if bestchi is None else np.array([oblateness])
    posang_chi = np.array([]) if bestchi is None else np.array([position_angle])
    chi2_best = np.array([]) if bestchi is None else np.array([bestchi])

    while len(f0_chi) < number:
        if verbose:
            progressbar_show(len(f0_chi), number, prefix='Ellipse fit:')

        chi2 = np.zeros(loop)
        f0 = center_f + dcenter_f*(2*np.random.RandomState().random(loop) - 1)
        g0 = center_g + dcenter_g*(2*np.random.RandomState().random(loop) - 1)
        a = equatorial_radius + dequatorial_radius*(2*np.random.RandomState().random(loop) - 1)
        obla = oblateness + doblateness*(2*np.random.RandomState().random(loop) - 1)
        obla[obla < 0], obla[obla > 1] = 0, 1
        phi_deg = position_angle + dposition_angle*(2*np.random.RandomState().random(loop) - 1)
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
        # if verbose:
        #     print('Elapsed time: {:.3f} seconds.'.format((controle_f2 - controle_f1).sec), end='\r')
        #     print(len(chi2[region]), len(chi2_best), end='\r')
        f0_chi = np.append(f0_chi, f0[region])
        g0_chi = np.append(g0_chi, g0[region])
        a_chi = np.append(a_chi, a[region])
        obla_chi = np.append(obla_chi, obla[region])
        posang_chi = np.append(posang_chi, phi_deg[region])

    if verbose:
        progressbar_show(number, number, prefix='Ellipse fit:')

    return [chi2_best, f0_chi, g0_chi, a_chi,  obla_chi, posang_chi]


def fit_to_limb(limb, fg, error, center_f=0, dcenter_f=0, center_g=0, dcenter_g=0, scale=1, dscale=0, loop=150000):
    """

    Parameters
    ----------
    limb : `sora.body.shape.Limb`
        Generic limb to fit.

    fg : `numpy.array`
        Matrix nx2 with the `xy` coordinates of each of the `n` points
        to fit the limb. See example below.

    error : `numpy.array`
        Array with n values corresponding to the error of each point.
        See `fg` parameter. It may be the same format as `fg`, thus
        being a vector error.

    center_f : `int`, `float`, default=0
        The coordinate in f of the ellipse center.

    center_g : `int`, `float`, default=0
        The coordinate in g of the ellipse center.

    dcenter_f : `int`, `float`
        Interval for coordinate f of the ellipse center.

    dcenter_g : `int`, `float`
        Interval for coordinate g of the ellipse center.

    scale : `number`
        Scale factor of the limb

    dscale : `number`
        Interval for scale

    loop : `int`, default=150000
        The number of centers to attempt fitting.

    Returns
    -------

    chisquare : `sora.ChiSquare`
        A ChiSquare object with all parameters.

    Examples
    ________

    fg = np.array([[-107.3,   57.8],
                   [ 103.7,   53.2],
                   [ -20.9,  172.4],
                   [   1.9,  171.9]])
    """
    if scale - np.absolute(dscale) <= 0:
        raise ValueError('Scale can not be 0 or negative. Please provide proper scale and dscale')
    if fg.ndim != 2 and fg.shape[1] != 2:
        raise ValueError("'fg' parameter must be a nx2 matrix")
    if len(error) != len(fg):
        raise ValueError("'error' and 'fg' must have the same number of points")
    if error.ndim != 1:
        error = np.linalg.norm(error, axis=-1)
    x0 = center_f + dcenter_f * (2 * np.random.random(loop) - 1)
    y0 = center_g + dcenter_g * (2 * np.random.random(loop) - 1)
    s0 = scale + dscale * (2 * np.random.random(loop) - 1)
    err2 = np.square(error)

    def calc_dist(item):
        residual = limb_radial_residual(limb, fg, center_f=x0[item], center_g=y0[item], scale=s0[item])
        return np.sum(residual ** 2 / err2)

    chi2 = np.array(list(map(calc_dist, tqdm(range(loop)))))

    return ChiSquare(chi2, center_f=x0, center_g=y0, scale=s0, npts=len(fg))


def fit_shape(occ, center_f=0, dcenter_f=0, center_g=0, dcenter_g=0, scale=1, dscale=0, loop=150000, sigma_result=1):
    """
    Parameters
    ----------
    occ : `sora.Occultation`
        The Occultation object with information to fit
    center_f : `int`, `float`, default=0
        The coordinate in f of the ellipse center.
    center_g : `int`, `float`, default=0
        The coordinate in g of the ellipse center.
    dcenter_f : `int`, `float`
        Interval for coordinate f of the ellipse center.
    dcenter_g : `int`, `float`
        Interval for coordinate g of the ellipse center.
    scale : `number`
        Scale factor of the limb
    dscale : `number`
        Interval for scale
    loop : `int`, default=150000
        The number of centers to attempt fitting.
    sigma_result : `int`, `float`
        Sigma value to be considered as result.

    Returns
    -------

    chisquare : `sora.ChiSquare`
        A ChiSquare object with all parameters.

    """
    orientation = occ.body.get_orientation(time=occ.tca)
    limb = occ.body.shape.get_limb(**orientation)
    chord_names, fg, error = occ.chords.get_limb_points()
    chisquare = fit_to_limb(limb, fg, error, center_f=center_f, dcenter_f=dcenter_f, center_g=center_g, dcenter_g=dcenter_g,
                            scale=scale, dscale=dscale, loop=loop)

    result_sigma = chisquare.get_nsigma(sigma=sigma_result)
    occ.fitted_params = {i: result_sigma[i] for i in ['center_f', 'center_g', 'scale']}
    radial_dispersion = limb_radial_residual(limb, fg, **{i: val[0] for i, val in occ.fitted_params.items()})
    xy = fg - np.array([[occ.fitted_params['center_f'][0]], [occ.fitted_params['center_g'][0]]]).T
    theta = np.arctan2(xy.T[1], xy.T[0])
    position_angle_point = Angle(90 * u.deg - theta * u.rad).wrap_at(360 * u.deg).degree
    occ.chi2_params = {'chord_name': chord_names, 'radial_dispersion': radial_dispersion,
                       'position_angle': position_angle_point, 'radial_error': np.linalg.norm(error, axis=-1),
                       'chi2_min': result_sigma['chi2_min'],
                       'nparam': chisquare.nparam, 'npts': chisquare.npts}
    return chisquare
