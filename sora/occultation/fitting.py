import astropy.units as u
import numpy as np
from astropy.time import Time

from sora.config.decorators import deprecated_alias

__all__ = ['fit_ellipse']


@deprecated_alias(pos_angle='position_angle', dpos_angle='dposition_angle', log='verbose')  # remove this line for v1.0
def fit_ellipse(*args, equatorial_radius, dequatorial_radius=0, center_f=0, dcenter_f=0, center_g=0,
                dcenter_g=0, oblateness=0, doblateness=0, position_angle=0, dposition_angle=0,
                loop=10000000, number_chi=10000, dchi_min=None, verbose=False, ellipse_error=0, sigma_result=1,
                sigma_samples=None, method='chisqr', threads=1, marching_grid=False):
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
    
    sigma_samples : `int`
        The minimum number of samples that lie within the defined `sigma` interval,
        by default None.
    
    method : `str`, default=`chisqr`
        Method used to perform the fit. Available methods are:
        `chisqr` : monte carlo computation method used in versions of SORA <= 0.2.1.
        `fastchi` : monte carlo computation method, allows multithreading and can be
         enhanced by marching grid algorithm.
        `least_squares` or `ls`: best fit done used levenberg marquardt convergence algorithm.
        `differential_evolution` : best fit done using genetic algorithms.
        All methods return a Chisquare object.
        
    threads : `int`
        Number of threads/workers used to perform parallel computations of the chi square 
        object. It works with all methods except `chisqr`, by default 1.

    marching_grid : `bool`
        When using in the methods:`least_squares`, `ls`, `differential_evolution`,
        or `fastchi`, ff `True`, it will use the marching grid approach to march 
        faster towards the best solution. When `False` sampling is done uniformly 
        within the boundings provided, by default True. By default False.

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
    from sora.stats import Parameters, least_squares, differential_evolution, fastchi, ellipse, ellipseError
    from warnings import warn

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

    if method not in ['chisqr', 'least_squares', 'ls', 'fastchi', 'differential_evolution']:
        warn(f'Invalid method `{method}` provided. Setting to default.')
        method = 'chisqr'

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
            if verbose:
                print('Elapsed time: {:.3f} seconds.'.format((controle_f2 - controle_f1).sec), end='\r')
                print(len(chi2[region]), len(chi2_best), end='\r')
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
                    minval=-np.inf if not eqradius_vary else (equatorial_radius - dequatorial_radius/2), 
                    maxval=np.inf if not eqradius_vary else (equatorial_radius + dequatorial_radius/2), 
                    free=eqradius_vary)
        
        center_f_vary = False if dcenter_f == 0 else True
        initial.add(name='center_f', value=center_f, 
                    minval=-np.inf if not center_f_vary else (center_f-dcenter_f/2), 
                    maxval=np.inf if not center_f_vary else (center_f+dcenter_f/2), free=center_f_vary)

        center_g_vary = False if dcenter_g == 0 else True
        initial.add(name='center_g', value=center_g, 
                    minval=-np.inf if not center_g_vary else (center_g-dcenter_g/2), 
                    maxval=np.inf if not center_g_vary else (center_g+dcenter_g/2), free=center_g_vary)

        oblateness_vary = False if doblateness == 0 else True
        ob_minval = 0 if ((oblateness - doblateness/2) < 0) else (oblateness - doblateness/2)
        ob_maxval = 1 if ((oblateness + doblateness/2) >1) else (oblateness + doblateness/2)
        if doblateness >= 0.5:
            ob_minval, ob_maxval = 0, 1
        initial.add(name='oblateness', value=oblateness, 
                    minval=-np.inf if not oblateness_vary else ob_minval, 
                    maxval=np.inf if not oblateness_vary else ob_maxval, free=oblateness_vary)

        pangle_vary = False if dposition_angle == 0 else True
        initial.add(name='position_angle', value=position_angle, 
                    minval=-np.inf if not pangle_vary else (position_angle - dposition_angle/2), 
                    maxval=np.inf if not pangle_vary else (position_angle + dposition_angle/2), 
                    free=pangle_vary)
        
        fi, gi, si = np.array(values).T
        if dchi_min is None:
            dchi_min = 3

        # run fastchi
        if (method == 'fastchi'):
            chi_result = fastchi(ellipseError, initial, args=(fi, gi, si, ellipse_error), samples=number_chi, sigma_range=dchi_min+sigma_result, sigma=sigma_result, 
                                 sigma_samples=sigma_samples, marching_grid=marching_grid, threads=threads, show_progress=(True if verbose else False),
                                 run_size = 10000)

        #run least_squares
        if (method == 'least_squares') or (method == 'ls'):
            result = least_squares(ellipseError, initial, args=(fi, gi, si, ellipse_error), algorithm='trf', sigma=sigma_result)
            chi_result = fastchi(ellipseError, result.params, args=(fi, gi, si, ellipse_error), samples=number_chi, sigma_range=dchi_min+sigma_result, sigma=sigma_result, 
                             sigma_samples=sigma_samples, threads=threads, show_progress=(True if verbose else False), marching_grid=marching_grid,
                             run_size = 10000)
        
        # run differential_evolution  
        if (method == 'differential_evolution'):
            result = differential_evolution(ellipseError, initial, args=(fi, gi, si, ellipse_error), sigma=sigma_result)
            chi_result = fastchi(ellipseError, result.params, args=(fi, gi, si, ellipse_error), samples=number_chi, sigma_range=dchi_min+sigma_result, sigma=sigma_result, 
                             sigma_samples=sigma_samples, threads=threads, show_progress=(True if verbose else False), marching_grid=marching_grid,
                             run_size = 10000)



        names, chi_len, chi_values = chi_result.var_names, len(chi_result.chi), chi_result.chi

        # condition to update chi values and samples with
        conv_result = (method == 'least_squares') or (method == 'ls') or (method == 'differential_evolution')
        if conv_result:
            chi_values = np.append(chi_values, result.chisqr)
            chi_len += 1

        if 'equatorial_radius' in names:
            a_chi = chi_result.samples[:][names.index('equatorial_radius')]
            if conv_result:
                a_chi = np.append(a_chi, result.params['equatorial_radius'].value)
        else:
            a_chi = np.repeat(initial['equatorial_radius'].value, chi_len)

        if 'center_f' in names:
            f0_chi = chi_result.samples[:][names.index('center_f')]
            if conv_result:
                f0_chi = np.append(f0_chi, result.params['center_f'].value)
        else:
            f0_chi = np.repeat(initial['center_f'].value, chi_len)

        if 'center_g' in names:
            g0_chi = chi_result.samples[:][names.index('center_g')]
            if conv_result:
                g0_chi = np.append(g0_chi, result.params['center_g'].value)
        else:
            g0_chi = np.repeat(initial['center_g'].value, chi_len)

        if 'oblateness' in names:
            obla_chi = chi_result.samples[:][names.index('oblateness')]
            if conv_result:
                obla_chi = np.append(obla_chi, result.params['oblateness'].value)
        else:
            obla_chi = np.repeat(initial['oblateness'].value, chi_len)
        
        if 'position_angle' in names:
            posang_chi = chi_result.samples[:][names.index('position_angle')]
            if conv_result:
                posang_chi = np.append(posang_chi, result.params['position_angle'].value)
        else:
            posang_chi = np.repeat(initial['position_angle'].value, chi_len)

        chisquare = ChiSquare(chi_values, len(values), center_f=f0_chi, center_g=g0_chi, equatorial_radius=a_chi,
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